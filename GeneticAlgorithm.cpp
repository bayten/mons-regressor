/* Copyright 2017 Baytekov Nikita */
#include "./include/GeneticAlgorithm.h"

template<typename T>
bool Individual<T>::operator== (const Individual<T>& comp_obj) {
    if (data == comp_obj.get_data() && score == comp_obj.get_score())
        return 1;
    return 0;
}

template<typename T>
Individual<T> Individual<T>::operator*(const Individual<T>& cross_obj) {
    Individual<T> new_obj = *this;
    return new_obj;
}

template<typename T>
float Individual<T>::operator%(const Individual<T>& geno_obj) {
    return 0.0;
}


template<template <typename IT> class I, typename T>
bool Population<I, T>::add_individual(I<T> add_obj) {
    int ind_num = ind_vec.get_size();
    T add_obj_data = add_obj.get_data();

    for (int i = 0; i < ind_num; i++)
        if (add_obj_data == ind_vec[i].get_data())
            return 0;
    ind_vec.append(add_obj);

    return 1;
}

template<template <typename IT> class I, typename T>
bool Population<I, T>::del_individual(int idx) {
    return ind_vec.erase(idx);
}

template<template <typename IT> class I, typename T>
float Population<I, T>::get_avg_score() const {
    int ind_num = ind_vec.get_size();
    float total_score = 0.0;
    for (int i = 0; i < ind_num; i++)
        total_score += ind_vec[i].get_score();

    return total_score;
}

template<template <typename IT> class I, typename T>
const I<T>& Population<I, T>::get_best_ind() const {
    int ind_num = ind_vec.get_size();
    int idx = 0;
    float max_score = ind_vec[idx];
    for (int i = 1; i < ind_num; i++)
        if (ind_vec[i].get_score() > max_score) {
            idx = i;
            max_score = ind_vec[idx];
        }

    return ind_vec[idx];
}


template<template<typename IT> class I, typename T>
Population<I, T> TournamentSelector<I, T>::select_population(const Population<I, T>& in_popul) {
    Population<I, T> out_popul;
    int popul_size = in_popul.get_popul_size();
    int needed_size = in_popul * this->popul_frac;

    for (int i = 0; i < needed_size; i++) {
        I<T> fst_candidate = in_popul[rand_r(time(0)) % popul_size];
        I<T> sec_candidate = in_popul[rand_r(time(0)) % popul_size];

        if (fst_candidate == sec_candidate) {
            i--;
            continue;
        } else if (fst_candidate.get_score() >= sec_candidate.get_score()) {
            if (!out_popul.add_individual(fst_candidate)) {
                i--;
                continue;
            }
        } else if (!out_popul.add_individual(sec_candidate)) {
            i--;
            continue;
        }
    }
    return out_popul;
}


template<template<typename IT> class I, typename T>
Population<I, T> RouletteSelector<I, T>::select_population(const Population<I, T>& in_popul) {
    Population<I, T> out_popul;
    int popul_size = in_popul.get_popul_size();
    int needed_size = in_popul * this->population_fraction;
    Vec<float> probabilities(popul_size);
    float total_score = 0;
    for (int i = 0; i < popul_size; i++) {
        probabilities[i] = out_popul[i].get_score();
        total_score += probabilities[i];
    }
    for (int i = 0; i < popul_size; i++)
        probabilities[i] /= total_score;

    for (int i = 0; i < needed_size; i++) {
        float rand_frac = static_cast<float>(rand_r(time(0))) / (RAND_MAX);
        int candidate_idx = 0;
        for (; candidate_idx < popul_size; candidate_idx++) {
            if (rand_frac < probabilities[candidate_idx])
                break;
            rand_frac -= probabilities[candidate_idx];
        }
        if (!out_popul.add_individual(in_popul[candidate_idx]))
            i--;
    }
    return out_popul;
}

template<template<typename IT> class I, typename T>
RankingSelector<I, T>::RankingSelector(int init_pfrac, int init_uniform) :
        GeneticSelector<I, T>(init_pfrac),
        uniform_thresh(init_uniform) {
}

template<template<typename IT> class I, typename T>
Population<I, T> RankingSelector<I, T>::select_population(const Population<I, T>& in_popul) {
    Population<I, T> out_popul;
    int popul_size = in_popul.get_popul_size();
    int needed_size = in_popul * this->popul_frac;
    Vec<float> probabilities(popul_size);

    for (int i = 0; i < popul_size; i++)
        probabilities[i] = in_popul.ind_vec[i].get_score();
    Vec<int> indices = probabilities.sort_indices();

    if (uniform_thresh > 1) {
        float val = 1.0/uniform_thresh;
        for (int i = 0; i < popul_size; i++)
            probabilities[i] = (popul_size-indices[i] <= uniform_thresh) ? val : 0;
    } else {
        float divisor = popul_size*(popul_size+1)/2.0;
        for (int i = 0; i < popul_size; i++)
            probabilities[i] = (popul_size-indices[i])/divisor;
    }

    for (int i = 0; i < needed_size; i++)
        if (!out_popul.add_individual(in_popul[rand_r(time(0)) % popul_size]))
            i--;
    return out_popul;
}


template<template<typename IT> class I, typename T>
Population<I, T> SigmaTruncSelector<I, T>::select_population(const Population<I, T>& in_popul) {
    Population<I, T> out_popul;
    int popul_size = in_popul.get_popul_size();
    int needed_size = in_popul * this->popul_frac;
    Vec<float> probabilities(popul_size);

    for (int i = 0; i < popul_size; i++)
        probabilities[i] = in_popul.ind_vec[i].get_score();
    Vec<int> indices = probabilities.sort_indices();

    for (int i = 0; i < popul_size; i++)
        if (popul_size-indices[i] <= needed_size)
            out_popul.add_individual(in_popul[i]);
    return out_popul;
}


template<template<typename IT> class I, typename T>
Population<I, T> PanmixiaBreeder<I, T>::breed_new_population(const Population<I, T>& in_popul) {
    Population<I, T> out_popul;
    int in_popul_size = in_popul.get_popul_size();
    if (in_popul_size <= 1)
        return in_popul;

    for (int i = 0; i < this->popul_num; i++) {
        I<T> fst_parent = in_popul[rand_r(time(0)) % in_popul_size];
        I<T> sec_parent = in_popul[rand_r(time(0)) % in_popul_size];
        while (sec_parent == fst_parent)
            sec_parent = in_popul[rand_r(time(0)) % in_popul_size];
        if (!out_popul.add_individual(fst_parent * sec_parent))
            i--;
    }
    return out_popul;
}

template<template <typename IT> class I, typename T>
InOutBreeder<I, T>::InOutBreeder(InOutBreederType init_type, int init_pnum) :
            GeneticBreeder<I, T>(init_pnum),
            breeder_type(init_type) {
}

template<template <typename IT> class I, typename T>
Population<I, T> InOutBreeder<I, T>::breed_new_population(const Population<I, T>& in_popul) {
    Population<I, T> out_popul;
    int in_popul_size = in_popul.get_popul_size();
    if (in_popul_size <= 1)
        return in_popul;

    for (int i = 0; i < this->popul_num; i++) {
        I<T> fst_parent = in_popul[rand_r(time(0)) % in_popul_size];
        I<T> sec_parent = find_match(in_popul, fst_parent);
        if (!out_popul.add_individual(fst_parent * sec_parent))
            i--;
    }
    return out_popul;
}

template<template <typename IT> class I, typename T>
const I<T>& InOutBreeder<I, T>::find_match(const Population<I, T>& in_popul, const I<T>& suitor) {
    int in_popul_size = in_popul.get_size();
    int partner_idx = (in_popul[0] == suitor) ? 1 : 0;
    float param_diff = 0.0, curr_param_diff = 0.0;
    float suitor_score = suitor.get_score();

    switch (breeder_type) {
        case kInPhenoType:
            param_diff = fabs(in_popul[partner_idx].get_score() - suitor_score);
            for (int i = (partner_idx+1); i < in_popul_size; i++) {
                curr_param_diff = fabs(in_popul[i].get_score() - suitor_score);
                if (curr_param_diff < param_diff && !(in_popul[i] == suitor)) {
                    partner_idx = i;
                    param_diff = curr_param_diff;
                }
            }
            break;

        case kInGenoType:
            param_diff = in_popul[partner_idx] % suitor;
            for (int i = (partner_idx+1); i < in_popul_size; i++) {
                curr_param_diff = in_popul[i] % suitor;
                if (curr_param_diff < param_diff && !(in_popul[i] == suitor)) {
                    partner_idx = i;
                    param_diff = curr_param_diff;
                }
            }
            break;

        case kOutPhenoType:
            param_diff = fabs(in_popul[partner_idx].get_score() - suitor_score);
            for (int i = (partner_idx+1); i < in_popul_size; i++) {
                curr_param_diff = fabs(in_popul[i].get_score() - suitor_score);
                if (curr_param_diff > param_diff && !(in_popul[i] == suitor)) {
                    partner_idx = i;
                    param_diff = curr_param_diff;
                }
            }
            break;

        case kOutGenoType:
            param_diff = in_popul[partner_idx] % suitor;
            for (int i = (partner_idx+1); i < in_popul_size; i++) {
                curr_param_diff = in_popul[i] % suitor;
                if (curr_param_diff > param_diff && !(in_popul[i] == suitor)) {
                    partner_idx = i;
                    param_diff = curr_param_diff;
                }
            }
            break;
    }
    return in_popul[partner_idx];
}

template<template <typename IT> class I, typename T>
GeneticAlgorithm<I, T>::GeneticAlgorithm(GeneticInitiator<I, T> init_initiator, \
                                         GeneticSelector<I, T>  init_selector, \
                                         GeneticBreeder<I, T>   init_breeder, \
                                         GeneticMutator<I, T>   init_mutator, \
                                         TerminationCriterionType init_term_crit,
                                         float init_term_crit_val) :
                initiator(init_initiator),
                selector(init_selector),
                breeder(init_breeder),
                mutator(init_mutator),
                term_crit(init_term_crit),
                term_crit_val(init_term_crit_val) {
}

template<template <typename IT> class I, typename T>
Vec<T> GeneticAlgorithm<I, T>::execute_ga() {
    Population<I, T> curr_population = initiator->get_init_population(sample_set);
    int64_t iter_count = 0L;

    while (true) {
        Population<I, T> breed_population = selector.select_population(curr_population);
        Population<I, T> new_population = breeder.breed_new_population(breed_population);
        new_population = mutator.mutate_population(new_population);

        if (check_term_crit(curr_population, new_population, iter_count))
            return new_population.get_popul_data;
    }
}

template<template <typename IT> class I, typename T>
bool GeneticAlgorithm<I, T>::check_term_crit(Population<I, T> curr_popul, \
                                             Population<I, T> new_popul, \
                                             int64_t iter_cnt) {
    switch (term_crit) {
        case kPopulConverged:
            if (new_popul.get_avg_score() - curr_popul.get_avg_score() < term_crit_val)
                return 1;
            return 0;

        case kBestConverged:
            float new_best_score = new_popul.get_best_ind().get_score();
            float curr_best_score = curr_popul.get_best_ind().get_score();
            if (new_best_score - curr_best_score < term_crit_val)
                return 1;
            return 0;

        case kMaxPopulNum:
            if (iter_cnt < static_cast<int64_t>(term_crit_val))
                return 1;
            return 0;
    }
}
