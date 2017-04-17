/* Copyright 2017 Baytekov Nikita */
#include "../include/default_types.h"
#include "../include/CollFamily.h"
#include "../include/SampleHandler.h"
#include "../include/CoveringHandler.h"
#include "../include/GeneticAlgorithm.h"

#ifndef INCLUDE_GENETICDUALIZER_H_
#define INCLUDE_GENETICDUALIZER_H_

template<typename S, typename T>
class GeneticDualizer;

template<typename S, typename T>
class GenDualInitiator : public GeneticInitiator<S, T> {
    GeneticDualizer<S, T>* my_parent;
    CoveringHandler<T> covering_handler;
 public:
    explicit GenDualInitiator(GeneticDualizer<S, T>* init_parent_ptr,
                              int init_num = 0);
    ~GenDualInitiator() {}

    Population<T> get_init_population(const SampleSet<S>& sample_set);
};

template<typename S, typename T>
class GenDualMutator : public GeneticMutator<T> {
    GeneticDualizer<S, T>* my_parent;
    CoveringHandler<T> covering_handler;
    float mutation_rate;
    int64_t curr_iter;
 public:
    explicit GenDualMutator(GeneticDualizer<S, T>* init_parent_ptr,
                            float init_pfrac = 0.0,
                            float init_mrate = 0.2);
    ~GenDualMutator() {}
    virtual Population<T> mutate_population(const Population<T>& in_popul);
 private:
    Population<T> recover_admissibility(const Population<T>& in_popul);
};


template<typename S, typename T>
class GeneticDualizer : public GeneticAlgorithm<S, T> {
    SampleHandler<S> sample_handler;
    ElColl<S> local_basis;
    SampleSet<S> basic;
    SampleSet<S> valid;
    int target_tag;

    Mat<bool> gen_matrix;

 public:
    explicit GeneticDualizer(int popul_size, float sel_frac, float max_mut,
                             float mut_frac, TerminationCriterionType init_tcrit = kPopulConverged,
                             float init_tcrit_val = 1.0);
    ~GeneticDualizer() {}

    void set_init_data(const SampleSet<S>& init_set,
                       const ElColl<S>& init_lb,
                       int init_target_tag);
    Vec<ElColl<S> > decode_collections(Vec< Vec<T> > in_code);
    virtual void update_costs(Population<T>* in_popul);

 private:
    float get_quality(const Chromosome<T>& chromo);
    const Mat<bool>& get_gen_matrix() const { return gen_matrix; }
};

#endif  // INCLUDE_GENETICDUALIZER_H_
