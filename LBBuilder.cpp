/* Copyright 2017 Baytekov Nikita */
#include "./include/LBBuilder.h"
#include <vector>
#include <algorithm>

template<class T>
RandomLBBuilder<T>::RandomLBBuilder() {
}

template<class T>
RandomLBBuilder<T>::~RandomLBBuilder() {
}

template<class T>
ElColl<T> RandomLBBuilder<T>::build_lb(const SampleSet<T>& train, int class_tag) {
    int total_num = train.get_total_size();
    int limit = static_cast<int>(2 * log(total_num) + 3);

    ElColl<T> local_basis;
    for (int i = 0; i < limit; i++) {
        ElClass<T> new_elclass = get_elclass(train[rand_r(time(0)) % total_num]);
        if (!local_basis.add(new_elclass))
            i--;
    }

    return local_basis;
}

template<class T>
ElClass<T> RandomLBBuilder<T>::get_elclass(const Vec<T>& rand_obj) {
    int num_of_features = rand_obj.get_size();
    int chosen_features = rand_r(time(0)) % num_of_features;
    std::vector<int> features = {};
    for (int i = 0; i < num_of_features; i++)
        features.push_back(i);
    random_shuffle(features.begin(), features.end());

    int init_cols[chosen_features] = {};
    int init_vals[chosen_features] = {};
    for (int i = 0; i < chosen_features; i++) {
        init_cols[i] = features[i];
        init_vals[i] = rand_obj[features[i]];
    }

    return ElClass<T>(chosen_features, init_cols, init_vals);
}
