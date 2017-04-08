/* Copyright 2017 Baytekov Nikita */
#include "../include/default_types.h"
#include "../include/CollFamily.h"
#include "../include/SampleHandler.h"
#include "../include/GeneticAlgorithm.h"

#ifndef INCLUDE_GENETICDUALIZER_H_
#define INCLUDE_GENETICDUALIZER_H_

template<typename T>
class CoveringInd : public Individual < ElColl<T> > {
 public:
    explicit CoveringInd(ElColl<T> init_coll) : Individual< ElColl<T> >(init_coll) {}
    ~CoveringInd() {}

    float get_quality(const SampleSet<T>& basic, \
                      const SampleSet<T>& valid, \
                      int target_tag);

    virtual Individual< ElColl<T> > operator*(const Individual<ElColl<T> >& cross_obj);
    virtual float operator%(const Individual<ElColl<T> >& geno_obj);
    virtual bool operator== (const Individual<T>& comp_obj);
};


template<typename T>
class CoveringPopul : public Population<CoveringInd, T> {
 public:
    explicit CoveringPopul(Vec< CoveringInd<T> > init_vec);
    CoveringPopul() {}
    ~CoveringPopul() {}

    virtual bool update_costs();
    virtual Population<CoveringInd, T>& operator+ (const Population<CoveringInd, T>& merge_obj);
};


enum CoveringFormType {
    kBinForm = 0,
    kIntForm = 1
};
template<typename T>
class GenDualInitiator : public GeneticInitiator<CoveringInd, T> {
    CoveringFormType covering_form;

 public:
    explicit GenDualInitiator(int init_num = 0, CoveringFormType init_form = kIntForm);
    ~GenDualInitiator() {}

    Population<CoveringInd, T> get_init_population(SampleSet<T> sample_set);

    // private:
    // build_covering();
    // make_covering_deadend();
    // get_covering_form();
};

template<typename T>
class GenDualMutator : public GeneticMutator<CoveringInd, T> {
 public:
    GenDualMutator();
    virtual ~GenDualMutator() = 0;
    virtual Population<CoveringInd, T> mutate_population(Population<CoveringInd, T> in_popul);
};

template<class T>
class GeneticDualizer : public GeneticAlgorithm<CoveringInd, T> {
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
