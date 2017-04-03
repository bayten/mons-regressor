/* Copyright 2017 Baytekov Nikita */
#include "../include/default_types.h"
#include "../include/containers.h"
#include "../include/SampleSet.h"

#ifndef INCLUDE_LBBUILDER_H_
#define INCLUDE_LBBUILDER_H_

// Abstract class for future realisations of BOOSTLO and TREELO

template<class T>
class LBBuilder {
 public:
    LBBuilder();
    virtual ~LBBuilder() = 0;
    virtual ElColl<T>& build_lb(const SampleSet<T>& train, int class_tag) = 0;
};


template<class T>
class RandomLBBuilder : public LBBuilder<T> {
 public:
    RandomLBBuilder();
    ~RandomLBBuilder();

    ElColl<T>& build_lb(const SampleSet<T>& train, int class_tag);
};
#endif  // INCLUDE_LBBUILDER_H_
