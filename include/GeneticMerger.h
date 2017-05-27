/* Copyright 2017 Baytekov Nikita */

#include "default_types.h"
#include "genetic_types.h"

#ifndef INCLUDE_GENETICMERGER_H_
#define INCLUDE_GENETICMERGER_H_


template<typename T>
class GeneticMerger {
 protected:
    bool is_max_sf;

 public:
    GeneticMerger(bool init_max_sf = true) : is_max_sf(init_max_sf) {}
    GeneticMerger(const GeneticMerger<T>& gm_obj) {}
    virtual ~GeneticMerger() {}

    virtual Population<T> merge_populations(const Population<T>& parents,
                                            const Population<T>& children) = 0;
};

template<typename T>
class SequentialMerger : public GeneticMerger<T> {
 public:
    SequentialMerger(bool init_max_sf = true) : GeneticMerger<T>(init_max_sf) {}
    SequentialMerger(const GeneticMerger<T>& gm_obj) {}
    virtual ~SequentialMerger() {}

    virtual Population<T> merge_populations(const Population<T>& parents,
                                            const Population<T>& children);
};



template<typename T>
Population<T> SequentialMerger<T>::merge_populations(const Population<T>& parents,
                                                     const Population<T>& children) {
    Population<T> out_popul = parents;

    int parents_num = parents.get_size();
    int children_num = children.get_size();

    Vec<float> scores(parents_num);
    for (int i = 0; i < parents_num; i++)
        scores[i] = parents[i].get_score();
    for (int i = 0; i < children_num; i++) {
        Vec<int> replace_idx;
        float my_score = children[i].get_score();

        for (int j = 0; j < parents_num; j++)
            if (this->is_max_sf == (scores[j] < my_score)) {
                replace_idx.append(j);
            }

        LOG_(trace) << "Final candidates to be replaced:" << replace_idx;
        if (!replace_idx.get_size())
            continue;
        int replaced_idx = replace_idx[rand() % replace_idx.get_size()];
        out_popul[replaced_idx] = children[i];
        scores[replaced_idx] = my_score;
    }

    return out_popul;
}

#endif  // INCLUDE_GENETICMERGER_H_
