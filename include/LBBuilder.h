/* Copyright 2017 Baytekov Nikita */
#include <ctime>
#include <vector>
#include <algorithm>
#include "../include/default_types.h"
#include "../include/CollFamily.h"
#include "../include/SampleSet.h"

#ifndef INCLUDE_LBBUILDER_H_
#define INCLUDE_LBBUILDER_H_

// Abstract class for future realisations of BOOSTLO and TREELO
template<typename T>
class LBBuilder {
 public:
    LBBuilder() {}
    virtual ~LBBuilder() {}
    virtual ElColl<T> build_lb(const SampleSet<T>& train, int class_tag) { return ElColl<T>(); }
};


template<typename T>
class RandomLBBuilder : public LBBuilder<T> {
 public:
    RandomLBBuilder() {}
    ~RandomLBBuilder() {}

    ElColl<T> build_lb(const SampleSet<T>& train, int class_tag);

 private:
    ElClass<T> get_elclass(const Vec<T>& rand_obj);
};

template<typename T>
ElColl<T> RandomLBBuilder<T>::build_lb(const SampleSet<T>& train, int class_tag) {
    int total_num = train.get_total_size();
    int limit = static_cast<int>(2 * log(total_num) + 3);
    unsigned int rand_seed = static_cast<unsigned int>(time(0));

    ElColl<T> local_basis;
    for (int i = 0; i < limit; i++) {
        ElClass<T> new_elclass(get_elclass(train[rand_r(&rand_seed) % total_num]));
        if (!local_basis.add(new_elclass))
            i--;
    }

    return local_basis;
}

template<typename T>
ElClass<T> RandomLBBuilder<T>::get_elclass(const Vec<T>& rand_obj) {
    int num_of_features = rand_obj.get_size();
    unsigned int rand_seed = static_cast<unsigned int>(time(0));
    int chosen_features = rand_r(&rand_seed) % num_of_features;
    std::vector<int> features = {};
    for (int i = 0; i < num_of_features; i++)
        features.push_back(i);
    random_shuffle(features.begin(), features.end());

    int* init_cols = new int[chosen_features];
    int* init_vals = new int[chosen_features];

    for (int i = 0; i < chosen_features; i++) {
        init_cols[i] = features[i];
        init_vals[i] = rand_obj[features[i]];
    }

    ElClass<T> out_ec(chosen_features, init_cols, init_vals);
    delete [] init_cols;
    delete [] init_vals;

    return out_ec;
}

#endif  // INCLUDE_LBBUILDER_H_
