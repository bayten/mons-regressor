/* Copyright 2017 Baytekov Nikita */
#ifndef INCLUDE_SAMPLESET_H_
#define INCLUDE_SAMPLESET_H_

#include <algorithm>
#include "../include/default_types.h"

template<class T>
class ClassSamples {
 private:
    Mat<T> objs;
    int class_tag;

 public:
    ClassSamples() {}
    ClassSamples(Mat<T> init_objs, int init_tag) :
                objs(init_objs), class_tag(init_tag) {}
    ~ClassSamples() {}

    void tearup(int tear_num, Mat<T>* torn_x, Vec<int>* torn_y);
    int get_size() const { return objs.get_sx(); }
    int get_tag() const { return class_tag; }

    Vec<T>& operator[](int index) { return objs[index]; }
    const Vec<T>& operator[](int index) const { return objs[index]; }
};


template<class T>
class SampleSet {
 private:
    Vec< ClassSamples<T> > samples;
    int total_size;

 public:
    explicit SampleSet(int size = 0) : samples(size), total_size(size) {}  // default constructor
    explicit SampleSet(Vec< ClassSamples<T> > init_samples);
    explicit SampleSet(const SampleSet& copy_obj) :  // copy constructor
                samples(copy_obj.samples),
                total_size(copy_obj.get_total_size()) {}

    explicit SampleSet(SampleSet&& move_obj) :  // move constructor
                samples(std::move(move_obj.get_samples())),
                total_size(move_obj.get_total_size()) {}
    ~SampleSet() {}

    void append(Mat<T> X, Vec<int> y);
    void shuffle();
    int get_class_num() const { return samples.get_size(); }
    int get_total_size() const { return total_size; }
    Vec< ClassSamples <T> >& get_samples() { return samples; }
    const Vec< ClassSamples <T> >& get_samples() const { return samples; }
    ClassSamples<T>& get_class(int index_tag);
    SampleSet<T> get_anticlass(int index_tag);
    bool delete_class(int index_tag);

    SampleSet<T>& operator= (const SampleSet<T>& copy_obj);  // copy assignment
    Vec<T>& operator[](int abs_index);
    const Vec<T>& operator[](int abs_index) const;
};


template<typename T>
void ClassSamples<T>::tearup(int tear_num, Mat<T>* torn_x, Vec<int>* torn_y) {
    std::random_shuffle(&(objs[0]), &(objs[-1]));
    *torn_x = objs.cutslice(0, tear_num);
    *torn_y = Vec<int>(tear_num);
    for (int j = 0; j < tear_num; j++)
        (*torn_y)[j] = class_tag;
}

template<typename T>
SampleSet<T>::SampleSet(Vec< ClassSamples<T> > init_samples) :
        samples(init_samples), total_size(0) {
    int class_num = samples.get_size();
    for (int i = 0; i < class_num; i++)
        total_size += samples[i].get_size();
}

template<typename T>
void SampleSet<T>::append(Mat<T> X, Vec<int> y) {
    int obj_num = y.get_size();
    int samples_size = samples.get_size();
    int i = 0;

    if (!samples_size) {  // if appending for the first time...
        samples.append(ClassSamples<T>(X[0], y[0]));
        i++;
    }

    for (; i < obj_num; i++) {
        bool was_found = 0;

        for (int j = 0; j < samples_size; j++) {  // searching within existing classes
            if (samples[j].class_tag == y[i]) {
                samples[j].objs.append(X[i]);
                was_found = 1;
                break;
            }
        }
        if (!was_found) {  // otherwise adding new class container =)
            samples.append(ClassSamples<T>(X[i], y[i]));
        }
    }
}

template<typename T>
void SampleSet<T>::shuffle() {
    int samples_size = samples.get_size();
    for (int i = 0; i < samples_size; i++) {
        std::random_shuffle(&(samples[i].obj[0]), &(samples[i].obj[-1]));
    }
}

template<typename T>
ClassSamples<T>& SampleSet<T>::get_class(int index_tag) {
    int samples_size = samples.get_size();
    for (int i = 0; i < samples_size; i++)
        if (samples[i].get_tag() == index_tag)
            return samples[i];

    return samples[0];
}

template<typename T>
SampleSet<T> SampleSet<T>::get_anticlass(int index_tag) {
    SampleSet<T> new_sample_set = *this;
    new_sample_set.delete_class(index_tag);
    return new_sample_set;
}

template<typename T>
bool SampleSet<T>::delete_class(int index_tag) {
    int samples_size = samples.get_size();
    for (int i = 0; i < samples_size; i++) {
        if (samples[i].get_tag() == index_tag) {
            samples.erase(i);
            return;
        }
    }
    return;
}

template<typename T>
Vec<T>& SampleSet<T>::operator[](int abs_index) {
    int samples_num = samples.get_size();
    for (int i = 0; i < samples_num; i++) {
        if (abs_index < samples[i].get_size())
            return samples[i][abs_index];

        abs_index -= samples[i].get_size();
    }

    return samples[-1][-1];
}

template<typename T>
const Vec<T>& SampleSet<T>::operator[](int abs_index) const {
    int samples_num = samples.get_size();
    for (int i = 0; i < samples_num; i++) {
        if (abs_index < samples[i].get_size())
            return samples[i][abs_index];

        abs_index -= samples[i].get_size();
    }

    return samples[-1][-1];
}

template<typename T>
SampleSet<T>& SampleSet<T>::operator= (const SampleSet<T>& copy_obj) {
    samples = copy_obj.get_samples();
    total_size = copy_obj.get_total_size();
    return (*this);
}

#endif  // INCLUDE_SAMPLESET_H_
