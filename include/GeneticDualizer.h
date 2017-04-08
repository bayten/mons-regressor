/* Copyright 2017 Baytekov Nikita */
#include "../include/default_types.h"
#include "../include/CollFamily.h"
#include "../include/SampleHandler.h"
#include "../include/GeneticAlgorithm.h"

#ifndef INCLUDE_GENETICDUALIZER_H_
#define INCLUDE_GENETICDUALIZER_H_

enum CoveringFormType {
    kBinForm = 0,
    kIntForm = 1
};
template<typename T>
class GenDualInitiator : public GeneticInitiator<T> {
    CoveringFormType covering_form;

 public:
    explicit GenDualInitiator(int init_num = 0, CoveringFormType init_form = kIntForm);
    ~GenDualInitiator() {}

    Population<T> get_init_population(SampleSet<T> sample_set);

    // private:
    // build_covering();
    // make_covering_deadend();
    // get_covering_form();
};

template<typename T>
class GenDualMutator : public GeneticMutator<T> {
 public:
    GenDualMutator();
    virtual ~GenDualMutator() = 0;
    virtual Population<T> mutate_population(Population<T> in_popul);
};

template<class T>
class GeneticDualizer : public GeneticAlgorithm<T> {
 private:
    SampleSet<T> basic;
    SampleSet<T> valid;
    int target_tag;

 public:
    GeneticDualizer(TerminationCriterionType init_term_crit = kPopulConverged,
                    float init_term_crit_val = 1.0);
    ~GeneticDualizer();

    virtual void set_sample_set(SampleSet<T> init_set);
    void set_target_tag(int init_target_tag) { target_tag = init_target_tag; }
};

#endif  // INCLUDE_GENETICDUALIZER_H_
