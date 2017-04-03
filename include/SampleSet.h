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
    ClassSamples(Mat<T> init_objs, int init_tag);
    ~ClassSamples();

    void tearup(int tear_num, Mat<T>* torn_x, Vec<int>* torn_y);
    int get_size() {
        return objs.get_sx();
    }
    int get_tag() {
        return class_tag;
    }
};

template<class T>
class SampleSet {
 private:
    Vec< ClassSamples<T> > samples;

 public:
    explicit SampleSet(int size = 0);  // default constructor
    explicit SampleSet(const SampleSet& copy_obj);  // copy constructor
    explicit SampleSet(SampleSet&& move_obj);  // move constructor
    ~SampleSet();

    void append(Mat<T> X, Vec<int> y);
    void shuffle();
    int get_size() {
        return samples.get_size();
    }
    ClassSamples<T>& get_sample(int vec_index) {
        return samples[vec_index];
    }

    SampleSet<T>& operator= (const SampleSet<T>& copy_obj);  // copy assignment
    SampleSet<T>& operator= (SampleSet<T>&& copy_obj);  // move assignment
    ClassSamples<T>& operator[](int index_tag);
};

#endif  // INCLUDE_SAMPLESET_H_
