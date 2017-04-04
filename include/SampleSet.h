/* Copyright 2017 Baytekov Nikita */
#ifndef INCLUDE_SAMPLESET_H_
#define INCLUDE_SAMPLESET_H_

#include <algorithm>
#include "./include/default_types.h"

template<class T>
class ClassSamples {
 private:
    Mat<T> objs;
    int class_tag;

 public:
    ClassSamples(Mat<T> init_objs, int init_tag) :
                objs(init_objs), class_tag(init_tag) {}
    ~ClassSamples() {}

    void tearup(int tear_num, Mat<T>* torn_x, Vec<int>* torn_y);
    int get_size() const { return objs.get_sx(); }
    int get_tag() const { return class_tag; }

    Vec<T>& operator[](int index) { return objs[index]; }
};

template<class T>
class SampleSet {
 private:
    Vec< ClassSamples<T> > samples;
    int total_size;

 public:
    explicit SampleSet(int size = 0) : samples(size), total_size(size) {}  // default constructor
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
    Vec< ClassSamples<T> >& get_samples() { return samples; }
    ClassSamples<T>& get_class(int index_tag);

    SampleSet<T>& operator= (const SampleSet<T>& copy_obj);  // copy assignment
    SampleSet<T>& operator= (SampleSet<T>&& copy_obj);  // move assignment
    Vec<T>& operator[](int abs_index);
};

#endif  // INCLUDE_SAMPLESET_H_
