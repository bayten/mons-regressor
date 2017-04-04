/* Copyright 2017 Baytekov Nikita */
#include "../include/default_types.h"
#include "../include/CollFamily.h"
#include "../include/SampleHandler.h"

#ifndef INCLUDE_GENETICDUALIZER_H_
#define INCLUDE_GENETICDUALIZER_H_

class Individual {
};

class Population {
};

template<class T>
class GeneticInitiator {
 public:
    GeneticInitiator();
    virtual ~GeneticInitiator();
    virtual Population get_init_population(SampleSet<T> sample_set);
};

class GeneticCrossover {
};

class GeneticMutator {
};

template<class T>
class GeneticDualizer {
 private:
    Population population;
    GeneticInitiator<T> initiator;
    GeneticCrossover cross_op;
    GeneticMutator mutator;

 public:
    GeneticDualizer();
    ~GeneticDualizer();

    Vec< CollFamily<T> > execute_ga(const Mat<T> & Train, int class_tag, \
                                    const ElColl<T> & LB);
};

#endif  // INCLUDE_GENETICDUALIZER_H_
