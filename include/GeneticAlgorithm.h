/* Copyright 2017 Baytekov Nikita */
#ifndef INCLUDE_GENETICALGORITHM_H_
#define INCLUDE_GENETICALGORITHM_H_

#include "../include/genetic_types.h"
#include "../include/GeneticSelector.h"
#include "../include/GeneticBreeder.h"

// TODO(Baytekov): Add weights!

template<typename S, typename T>
class GeneticInitiator {
 protected:
    int popul_num;

 public:
    explicit GeneticInitiator(int init_num = 0) : popul_num(init_num) {}
    GeneticInitiator(const GeneticInitiator<S, T>& gi_obj) : popul_num(gi_obj.popul_num) {}
    virtual ~GeneticInitiator() {}

    virtual Population<T> get_init_population(const SampleSet<S>& sample_set) = 0;
};


template<typename T>
class GeneticMutator {
 protected:
    float popul_frac;

 public:
    explicit GeneticMutator(int init_pfrac = 0.0) : popul_frac(init_pfrac) {}
    GeneticMutator(const GeneticMutator<T>& gm_obj) : popul_frac(gm_obj.popul_frac) {}
    virtual ~GeneticMutator() {}

    virtual Population<T> mutate_population(const Population<T>& in_popul) = 0;
};


enum TerminationCriterionType {
    kPopulConverged = 0,  // term_crit_val = min delta between two populations
    kBestConverged  = 1,  // term_crit_val = min delta between two fitnesses of best individual
    kMaxPopulNum    = 2   // term_crit_val = amount of populations
};
template<typename S, typename T>
class GeneticAlgorithm {
 protected:
    SampleSet<S> sample_set;

    GeneticInitiator<S, T>* initiator;
    GeneticSelector<T>* selector;
    GeneticBreeder<T>* breeder;
    GeneticMutator<T>* mutator;

    TerminationCriterionType term_crit;
    float term_crit_val;

 public:
    GeneticAlgorithm(GeneticInitiator<S, T>* init_initiator, GeneticSelector<T>* init_selector,
                     GeneticBreeder<T>* init_breeder, GeneticMutator<T>* init_mutator,
                     TerminationCriterionType init_tcrit = kPopulConverged,
                     float init_tcrit_val = 1.0);
    virtual ~GeneticAlgorithm() {}

    virtual Vec< Vec<T> > execute_ga();
    virtual bool check_term_crit(Population<T> curr_popul, Population<T> new_popul,
                                 int64_t iter_cnt);

     virtual void set_sample_set(SampleSet<S> init_set) { sample_set = init_set; }
     virtual void update_costs(Population<T>* in_popul) = 0;
};



template<typename S, typename T>
GeneticAlgorithm<S, T>::GeneticAlgorithm(GeneticInitiator<S, T>* init_initiator, \
                                         GeneticSelector<T>*  init_selector, \
                                         GeneticBreeder<T>*   init_breeder, \
                                         GeneticMutator<T>*   init_mutator, \
                                         TerminationCriterionType init_term_crit,
                                         float init_term_crit_val) :
                initiator(init_initiator),
                selector(init_selector),
                breeder(init_breeder),
                mutator(init_mutator),
                term_crit(init_term_crit),
                term_crit_val(init_term_crit_val) {
    LOG_(debug) << "Created GeneticAlgorithm instance (by usual constructor).";
}

template<typename S, typename T>
Vec< Vec<T> > GeneticAlgorithm<S, T>::execute_ga() {
    Population<T> curr_population = initiator->get_init_population(sample_set);
    int64_t iter_count = 0L;

    while (true) {
        update_costs(&curr_population);

        Population<T> breed_population = selector->select_population(curr_population);
        Population<T> new_population = breeder->breed_new_population(breed_population);

        new_population = mutator->mutate_population(new_population);
        update_costs(&new_population);

        if (check_term_crit(curr_population, new_population, iter_count))
            return new_population.get_popul_data();
    }
}

template<typename S, typename T>
bool GeneticAlgorithm<S, T>::check_term_crit(Population<T> curr_popul, \
                                             Population<T> new_popul, \
                                             int64_t iter_cnt) {
    switch (term_crit) {
        case kPopulConverged:
        {
            if (new_popul.get_avg_score() - curr_popul.get_avg_score() < term_crit_val)
                return 1;
            return 0;
        }

        case kBestConverged:
        {
            float new_best_score = new_popul.get_best_ind().get_score();
            float curr_best_score = curr_popul.get_best_ind().get_score();
            if (new_best_score - curr_best_score < term_crit_val)
                return 1;
            return 0;
        }

        case kMaxPopulNum:
        {
            if (iter_cnt < static_cast<int64_t>(term_crit_val))
                return 1;
            return 0;
        }

        default:
            return 0;
    }
    return 0;
}

#endif  // INCLUDE_GENETICALGORITHM_H_
