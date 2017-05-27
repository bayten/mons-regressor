/* Copyright 2017 Baytekov Nikita */

#include "default_types.h"
#include "MONSClassifier.h"

#ifndef INCLUDE_MONSREGRESSOR_H_
#define INCLUDE_MONSREGRESSOR_H_

enum RegressionCastType {
    kPhenoType  = 0,
    kGenoType   = 1,
    kGenoPhenoType = 2
    // kCustomType = 3  - TODO: make custom selection of features by bool vec
};
enum MetricType {  // TODO: Add more metrics!
    // kManhattanType = 0,
    kEuclideType   = 1,
    // kCosineType    = 2
}
template<typename S, typename T, typename U>
class MONSRegressor {
    MONSClassifier<S, T> mons_instance;
    RegressionCastType reg_cast;
    MetricType metric;

    SampleSet<S, U> original_sset;
    SampleSet<S, int> casted_sset;

 public:
    MONSRegressor();
    ~MONSRegressor();

    fit();
    predict();

    SampleSet<S, int> cast_sset(SampleSet<S, U> in_sset);
};



#endif  // INCLUDE_MONSREGRESSOR_H_
