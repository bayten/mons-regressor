/* Copyright 2017 Baytekov Nikita */

#include "../include/genetic_types.h"
#include "../include/SampleSet.h"

#ifndef INCLUDE_GENETICSELECTOR_H_
#define INCLUDE_GENETICSELECTOR_H_

template<typename T>
class GeneticSelector {
 protected:
    float popul_frac;
 public:
    explicit GeneticSelector(int init_pfrac = 0.5) : popul_frac(init_pfrac) {}
    virtual ~GeneticSelector() {}
    virtual Population<T> select_population(const Population<T>& in_popul) = 0;
};


template<typename T>
class TournamentSelector : public GeneticSelector<T> {
 public:
    explicit TournamentSelector(int init_pfrac = 0.5) : GeneticSelector<T>(init_pfrac) {}
    virtual ~TournamentSelector() {}
    virtual Population<T> select_population(const Population<T>& in_popul);
};

template<typename T>
class RouletteSelector : public GeneticSelector<T> {
 public:
    explicit RouletteSelector(int init_pfrac = 0.5) : GeneticSelector<T>(init_pfrac) {}
    virtual ~RouletteSelector() {}
    virtual Population<T> select_population(const Population<T>& in_popul);
};

template<typename T>
class RankingSelector : public GeneticSelector<T> {
    int uniform_thresh;
 public:
    explicit RankingSelector(int init_pfrac = 0.5, int init_uniform = 0);
    virtual ~RankingSelector() {}
    virtual Population<T> select_population(const Population<T>& in_popul);
};

template<typename T>
class SigmaTruncSelector : public GeneticSelector<T> {
 public:
    explicit SigmaTruncSelector(int init_pfrac = 0.5) : GeneticSelector<T>(init_pfrac) {}
    virtual ~SigmaTruncSelector() {}
    virtual Population<T> select_population(const Population<T>& in_popul);
};



template<typename T>
Population<T> TournamentSelector<T>::select_population(const Population<T>& in_popul) {
    Population<T> out_popul;
    int popul_size = in_popul.get_size();
    int needed_size = in_popul * this->popul_frac;
    unsigned int rand_seed = static_cast<unsigned int>(time(0));

    for (int i = 0; i < needed_size; i++) {
        Chromosome<T> fst_chromo = in_popul[rand_r(&rand_seed) % popul_size];
        Chromosome<T> sec_chromo = in_popul[rand_r(&rand_seed) % popul_size];

        if (fst_chromo == sec_chromo) {
            i--;
            continue;
        } else if (fst_chromo.get_score() >= sec_chromo.get_score()) {
            if (!out_popul.add_chromo(fst_chromo)) {
                i--;
                continue;
            }
        } else if (!out_popul.add_chromo(sec_chromo)) {
            i--;
            continue;
        }
    }
    return out_popul;
}


template<typename T>
Population<T> RouletteSelector<T>::select_population(const Population<T>& in_popul) {
    Population<T> out_popul;
    int popul_size = in_popul.get_size();
    int needed_size = popul_size * this->popul_frac;
    unsigned int rand_seed = static_cast<unsigned int>(time(0));
    Vec<float> probs(popul_size);  // probs will be stored here
    float total_score = 0;
    for (int i = 0; i < popul_size; i++) {
        probs[i] = out_popul[i].get_score();
        total_score += probs[i];
    }
    for (int i = 0; i < popul_size; i++)
        probs[i] /= total_score;

    for (int i = 0; i < needed_size; i++) {
        float rand_frac = static_cast<float>(rand_r(&rand_seed)) / (RAND_MAX);
        int chromo_idx = 0;
        for (; chromo_idx < popul_size; chromo_idx++) {
            if (rand_frac < probs[chromo_idx])
                break;
            rand_frac -= probs[chromo_idx];
        }
        if (!out_popul.add_chromo(in_popul[chromo_idx]))
            i--;
    }
    return out_popul;
}

template<typename T>
RankingSelector<T>::RankingSelector(int init_pfrac, int init_uniform) :
        GeneticSelector<T>(init_pfrac), uniform_thresh(init_uniform) {
}

template<typename T>
Population<T> RankingSelector<T>::select_population(const Population<T>& in_popul) {
    Population<T> out_popul;
    int popul_size = in_popul.get_size();
    int needed_size = in_popul * this->popul_frac;
    Vec<float> probs(popul_size);

    for (int i = 0; i < popul_size; i++)
        probs[i] = in_popul.ind_vec[i].get_score();
    Vec<int> indices = probs.sort_indices();

    if (uniform_thresh > 1) {
        float val = 1.0/uniform_thresh;
        for (int i = 0; i < popul_size; i++)
            probs[i] = (popul_size-indices[i] <= uniform_thresh) ? val : 0;
    } else {
        float divisor = popul_size*(popul_size+1)/2.0;
        for (int i = 0; i < popul_size; i++)
            probs[i] = (popul_size-indices[i])/divisor;
    }

    for (int i = 0; i < needed_size; i++)
        if (!out_popul.add_chromo(in_popul[rand_r(time(0)) % popul_size]))
            i--;
    return out_popul;
}


template<typename T>
Population<T> SigmaTruncSelector<T>::select_population(const Population<T>& in_popul) {
    Population<T> out_popul;
    int popul_size = in_popul.get_size();
    int needed_size = in_popul * this->popul_frac;
    Vec<float> probs(popul_size);

    for (int i = 0; i < popul_size; i++)
        probs[i] = in_popul.ind_vec[i].get_score();
    Vec<int> indices = probs.sort_indices();

    for (int i = 0; i < popul_size; i++)
        if (popul_size-indices[i] <= needed_size)
            out_popul.add_chromo(in_popul[i]);
    return out_popul;
}

#endif  // INCLUDE_GENETICSELECTOR_H_
