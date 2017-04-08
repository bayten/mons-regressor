/* Copyright 2017 Baytekov Nikita */
#ifndef INCLUDE_GENETICALGORITHM_H_
#define INCLUDE_GENETICALGORITHM_H_

#include "../include/default_types.h"
#include "../include/SampleHandler.h"

template<typename T>
class Individual {
 protected:
    T data;
    float score;

 public:
    explicit Individual(T init_data, float init_score = 0.0) : data(init_data), score(init_score) {}
    ~Individual() {}

    float get_score() const { return score; }
    T get_data() const { return data; }
    void set_score(float score_val) { score = score_val; }

    virtual Individual<T> operator*(const Individual<T>& cross_obj);  // crossover op
    virtual float operator%(const Individual<T>& geno_obj);  // genotypes comparation op
    virtual bool operator== (const Individual<T>& comp_obj);
};

template<template <typename IT> class I, typename T>
class Population {
    Vec< I<T> > ind_vec;

 public:
    explicit Population(Vec< I<T> > init_vec) : ind_vec(init_vec) {}
    Population() {}
    ~Population() {}

    virtual bool update_costs() { return 0; }

    bool add_individual(I<T> add_obj);
    bool del_individual(int idx);

    const Vec<T>& get_popul_data() const { return ind_vec; }
    int get_popul_size() const { return ind_vec.get_size(); }

    float get_avg_score() const;
    const I<T>& get_best_ind() const;

    I<T>& operator[] (int idx) { return ind_vec[idx]; }
    virtual Population<I, T>& operator+ (const Population<I, T>& merge_obj) { return (*this); }
};


template<template <typename IT> class I, typename T>
class GeneticInitiator {
    int popul_num;

 public:
    explicit GeneticInitiator(int init_num = 0) : popul_num(init_num) {}
    virtual ~GeneticInitiator() = 0;
    virtual Population<I, T> get_init_population(SampleSet<T> sample_set) = 0;
};


template<template <typename IT> class I, typename T>
class GeneticSelector {
    float popul_frac;
 public:
    explicit GeneticSelector(int init_pfrac = 0.5) : popul_frac(init_pfrac) {}
    virtual ~GeneticSelector() = 0;
    virtual Population<I, T> select_population(const Population<I, T>& in_popul) = 0;
};

template<template <typename IT> class I, typename T>
class TournamentSelector : public GeneticSelector<I, T> {
 public:
    explicit TournamentSelector(int init_pfrac = 0.5) : GeneticSelector<I, T>(init_pfrac) {}
    virtual ~TournamentSelector();
    virtual Population<I, T> select_population(const Population<I, T>& in_popul);
};

template<template <typename IT> class I, typename T>
class RouletteSelector : public GeneticSelector<I, T> {
 public:
    explicit RouletteSelector(int init_pfrac = 0.5) : GeneticSelector<I, T>(init_pfrac) {}
    virtual ~RouletteSelector();
    virtual Population<I, T> select_population(const Population<I, T>& in_popul);
};

template<template <typename IT> class I, typename T>
class RankingSelector : public GeneticSelector<I, T> {
    int uniform_thresh;
 public:
    explicit RankingSelector(int init_pfrac = 0.5, int init_uniform = 0);
    virtual ~RankingSelector();
    virtual Population<I, T> select_population(const Population<I, T>& in_popul);
};

template<template <typename IT> class I, typename T>
class SigmaTruncSelector : public GeneticSelector<I, T> {
 public:
    explicit SigmaTruncSelector(int init_pfrac = 0.5) : GeneticSelector<I, T>(init_pfrac) {}
    virtual ~SigmaTruncSelector();
    virtual Population<I, T> select_population(const Population<I, T>& in_popul);
};


template<template <typename IT> class I, typename T>
class GeneticBreeder {
    int popul_num;
 public:
    explicit GeneticBreeder(int init_pnum = 0) : popul_num(init_pnum) {}
    virtual ~GeneticBreeder() = 0;
    virtual Population<I, T> breed_new_population(const Population<I, T>& in_popul) = 0;
};

template<template <typename IT> class I, typename T>
class PanmixiaBreeder : public GeneticBreeder<I, T> {
 public:
    explicit PanmixiaBreeder(int init_pnum = 0) : GeneticBreeder<I, T>(init_pnum) {}
    ~PanmixiaBreeder();
    virtual Population<I, T> breed_new_population(const Population<I, T>& in_popul);
};

enum InOutBreederType {
    kInPhenoType  = 0,
    kInGenoType   = 1,
    kOutPhenoType = 2,
    kOutGenoType  = 3
};
template<template <typename IT> class I, typename T>
class InOutBreeder : public GeneticBreeder<I, T> {
    InOutBreederType breeder_type;

 public:
    explicit InOutBreeder(InOutBreederType init_type = kInPhenoType, int init_pnum = 0);
    ~InOutBreeder();
    virtual Population<I, T> breed_new_population(const Population<I, T>& in_popul);

 private:
    virtual const I<T>& find_match(const Population<I, T>& in_popul, const I<T>& suitor);
};


template<template <typename IT> class I, typename T>
class GeneticMutator {
    float popul_frac;
 public:
    explicit GeneticMutator(int init_pfrac = 0.0) : popul_frac(init_pfrac) {}
    virtual ~GeneticMutator() = 0;
    virtual Population<I, T> mutate_population(Population<I, T> in_popul) = 0;
};


enum TerminationCriterionType {
    kPopulConverged = 0,  // term_crit_val = min delta between two populations
    kBestConverged  = 1,  // term_crit_val = min delta between two fitnesses of best individual
    kMaxPopulNum    = 2   // term_crit_val = amount of populations
    // kMaxCompTime    = 3   // term_crit_val = amount of ticks for computation
};
template<template <typename IT> class I, typename T>
class GeneticAlgorithm {
 protected:
    SampleSet<T> sample_set;

    GeneticInitiator<I, T> initiator;
    GeneticSelector<I, T> selector;
    GeneticBreeder<I, T> breeder;
    GeneticMutator<I, T> mutator;

    TerminationCriterionType term_crit;
    float term_crit_val;

 public:
    GeneticAlgorithm(GeneticInitiator<I, T> init_initiator, \
                     GeneticSelector<I, T>  init_selector, \
                     GeneticBreeder<I, T>   init_breeder, \
                     GeneticMutator<I, T>   init_mutator, \
                     TerminationCriterionType init_term_crit = kPopulConverged,
                     float init_term_crit_val = 1.0);
    virtual ~GeneticAlgorithm();

    virtual void set_sample_set(SampleSet<T> init_set) { sample_set = init_set; }
    virtual Vec<T> execute_ga();
    virtual bool check_term_crit(Population<I, T> curr_popul, \
                                 Population<I, T> new_popul, \
                                 int64_t iter_cnt);
};

#endif  // INCLUDE_GENETICALGORITHM_H_
