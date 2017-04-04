/* Copyright 2017 Baytekov Nikita */
#include "../include/default_types.h"
#include "../include/CollFamily.h"
#include "../include/SampleHandler.h"
#include "../include/Population.h"

#ifndef INCLUDE_GENETICDUALIZER_H_
#define INCLUDE_GENETICDUALIZER_H_

template<class T>
class GeneticInitiator {
 public:
    GeneticInitiator();
    virtual ~GeneticInitiator() = 0;
    virtual Population<T> get_init_population(SampleSet<T> sample_set) = 0;
};


enum CoveringFormType {
    kBinForm = 0,
    kIntForm = 1
};

// TODO(Baytekov): How to initiate GA correctly?
template<class T>
class DefaultGeneticInitiator : public GeneticInitiator<T> {
    CoveringFormType covering_form;

 public:
    explicit DefaultGeneticInitiator(CoveringFormType init_form = kIntForm):
            covering_form(init_form) {}
    ~DefaultGeneticInitiator() {}

    Population<T> get_init_population(SampleSet<T> sample_set);

    // private:
    // build_covering();
    // make_covering_deadend();
    // get_covering_form();
};

template<class T>
class GeneticBreeder {
 public:
    GeneticBreeder();
    virtual ~GeneticBreeder() = 0;
    virtual Population<T> breed_new_population(Population<T> in_population) = 0;
};

template<class T>
class GeneticMutator {
 public:
    GeneticMutator();
    virtual ~GeneticMutator() = 0;
    virtual Population<T> mutate_population(Population<T> in_population) = 0;
};

// TODO(Baytekov): Tournament and Roulette Selectors
template<class T>
class GeneticSelector {  // Abstract class for further implementations
};

template<class T>
class GeneticDualizer {
 private:
    Population<T> population;
    SampleSet<T> basic;
    SampleSet<T> valid;
    int target_tag;

    GeneticInitiator<T> initiator;
    GeneticBreeder<T> breeder;
    GeneticMutator<T> mutator;
    GeneticSelector<T> selector;

 public:
    GeneticDualizer();
    ~GeneticDualizer();

    Vec< CollFamily<T> > execute_ga(const SampleSet<T>& Train, int init_tag, \
                                    const ElColl<T> & LB);
};

#endif  // INCLUDE_GENETICDUALIZER_H_
