/* Copyright 2017 Baytekov Nikita */
#include "default_types.h"
#include "log_wrapper.h"

#ifndef INCLUDE_ML_METRICS_H_
#define INCLUDE_ML_METRICS_H_

namespace ml_metrics {

template<typename T>
float ml_accuracy(Vec<T> fst, Vec<T> sec) {
    if (fst.get_size() != sec.get_size()) {
        LOG_(error) << "Cannot compare target vectors: different sizes!";
        return -1.0;
    }

    int sz = fst.get_size();
    float res = 0.0;
    for (int i = 0; i < sz; i++)
        res += fst[i] == sec[i];
    return res/sz;
}

template<typename T>
float ml_mse(Vec<T> fst, Vec<T> sec) {
    int sz = fst.get_size();
    float res = 0.0;
    for (int i = 0; i < sz; i++)
        res += (fst[i] - sec[i])*(fst[i] - sec[i]);
    return res/float(sz);
}

}

#endif  // INCLUDE_ML_METRICS_H_
