/* Copyright 2017 Baytekov Nikita */
#include "../include/default_types.h"
#include "../include/CollFamily.h"
#include "../include/SampleSet.h"
#ifndef INCLUDE_SAMPLEHANDLER_H_
#define INCLUDE_SAMPLEHANDLER_H_

template<class T>
class SampleHandler {
 public:
    SampleHandler();
    ~SampleHandler();

    SampleSet<T> make_samples(const Mat<T>& X, const Vec<int>& y, bool mix = true);
    int make_train_and_valid(const SampleSet<T>& sample_set, \
                             SampleSet<T>* train, \
                             SampleSet<T>* valid);
};

#endif  // INCLUDE_SAMPLEHANDLER_H_
