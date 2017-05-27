/* Copyright 2017 Baytekov Nikita */
#include <ctime>
#include <vector>
#include <algorithm>
#include "default_types.h"
#include "CollFamily.h"
#include "SampleSet.h"

#ifndef INCLUDE_LBBUILDER_H_
#define INCLUDE_LBBUILDER_H_

// Abstract class for future realisations of BOOSTLO and TREELO
template<typename T, typename S>
class LBBuilder {
 public:
    LBBuilder() {}
    virtual ~LBBuilder() {}
    virtual ElColl<T> build_lb(const SampleSet<T, S>& train, S group_tag) = 0;
};


template<typename T>
class RandomLBBuilder : public LBBuilder<T, int> {
 public:
    RandomLBBuilder() { LOG_(trace) << "Created RandomLBBuilder instance (by usual constructor)"; }
    ~RandomLBBuilder() {}

    ElColl<T> build_lb(const SampleSet<T, int>& train, int group_tag);

 private:
    ElClass<T> get_rand_ec(const Vec<T>& rand_obj);
    ElClass<T> get_diff_ec(const Vec<T>& fst_obj, const Vec<T>& sec_obj);
};



template<typename T>
ElColl<T> RandomLBBuilder<T>::build_lb(const SampleSet<T, int>& train, int class_tag) {
    LOG_(trace) << "Building local basis...";

    ElColl<T> local_basis;
    const GroupSamples<T, int>& tag_class = train.get_group(class_tag);
    SampleSet<T, int> tag_anticlass(train.get_antigroup(class_tag));
    int tag_class_size = tag_class.get_size();
    int limit = static_cast<int>(3*sqrt(tag_class_size))+3;

    local_basis = ElColl<T>();
    std::vector<int>idx(tag_class_size);
    std::iota(std::begin(idx), std::end(idx), 0);
    std::random_shuffle(std::begin(idx), std::end(idx));

    LOG_(trace) << "Working with train set: " << train;
    LOG_(trace) << "Need to add " << limit << " random elementary classifiers.";
    for (int i = 0; local_basis.get_size() < limit; i++) {
        ElClass<T> new_elclass(get_rand_ec(tag_class[idx[i]]));
        LOG_(trace) << "New elementary classifier to add(" << new_elclass << ")...";
        if (!local_basis.add(new_elclass)) {
            LOG_(trace) << "...was rejected.";
        } else {
            LOG_(trace) << "...was added to local basis!";
        }
        if (i == tag_class_size) {
            LOG_(warning) << "Reached the end of train objects - going to the start";
            i = 0;
        }
    }

    LOG_(trace) << "Accomplishing local basis... ";
    int tclass_size = tag_class.get_size();
    int tanticlass_size = tag_anticlass.get_total_size();

    for (int i = 0; i < tclass_size; i++) {
        for (int j = 0; j < tanticlass_size; j++) {
            Vec<bool> class_vec = local_basis.apply_to_object(tag_class[i]);
            Vec<bool> aclass_vec = local_basis.apply_to_object(tag_anticlass[j]);
            bool was_found = 0;
            int lb_size = local_basis.get_size();
            for (int k = 0; k < lb_size; k++) {
                if(class_vec[k] && !aclass_vec[k]) {
                    was_found = 1;
                    break;
                }
            }

            if(!was_found) {
                LOG_(trace) << "Zero row was found at " << i*tanticlass_size+j << " entry.";
                ElClass<T> new_elclass(get_diff_ec(tag_class[i], tag_anticlass[j]));
                if (!local_basis.add(new_elclass)) {
                    LOG_(error) << "Can't add " << new_elclass << " to distinguish objects!";
                }
            }
        }
    }

    return local_basis;
}

template<typename T>
ElClass<T> RandomLBBuilder<T>::get_rand_ec(const Vec<T>& rand_obj) {
    int num_of_features = rand_obj.get_size();
    int chosen_features = rand() % (std::max(1, static_cast<int>(0.8 * num_of_features))) + 1;

    std::vector<int> features(num_of_features);
    std::iota(std::begin(features), std::end(features), 0);
    std::random_shuffle(std::begin(features), std::end(features));

    int* init_cols = new int[chosen_features];
    T* init_vals = new T[chosen_features];

    for (int i = 0; i < chosen_features; i++) {
        init_cols[i] = features[i];
        init_vals[i] = rand_obj[features[i]];
    }

    ElClass<T> out_ec(chosen_features, init_cols, init_vals);
    delete [] init_cols;
    delete [] init_vals;

    return out_ec;
}

template<typename T>
ElClass<T> RandomLBBuilder<T>::get_diff_ec(const Vec<T>& fst_obj, const Vec<T>& sec_obj) {
    LOG_(trace) << "Computing differential elementary classifier...";
    LOG_(trace) << "Target objects: " << fst_obj << " and " << sec_obj;
    int num_of_features = fst_obj.get_size();
    Vec<int> diffs;
    Vec<int> equals;
    for (int i = 0; i < num_of_features; i++)
        if (fst_obj[i] != sec_obj[i])
            diffs.append(i);
        else
            equals.append(i);

    int df_size = diffs.get_size();
    int eq_size = equals.get_size();
    std::random_shuffle(&(diffs[0]), &(diffs[-1])+1);

    int df_features = rand() % (std::max(1, df_size-1)) + 1;
    int eq_features = (eq_size) ? rand() % (std::max(1, static_cast<int>(0.4 * eq_size))) + 1 : 0;
    LOG_(trace) << "Taking " << df_features << "+" << eq_features << " features";

    int* init_cols = new int[df_features + eq_features];
    T* init_vals = new T[df_features + eq_features];

    for (int i = 0; i < df_features; i++) {
        init_cols[i] = diffs[i];
        init_vals[i] = fst_obj[diffs[i]];
    }

    if (eq_size) {
        std::random_shuffle(&(equals[0]), &(equals[-1])+1);
        for (int i = 0; i < eq_features; i++) {
            init_cols[df_features + i] = equals[i];
            init_vals[df_features + i] = fst_obj[equals[i]];
        }
    }

    ElClass<T> out_ec(df_features + eq_features, init_cols, init_vals);

    delete [] init_cols;
    delete [] init_vals;

    if(!out_ec.apply_to_object(fst_obj) || out_ec.apply_to_object(sec_obj)) {
        LOG_(error) << "Incorrect differential elementary classifier!!!";
        LOG_(error) << "Fst obj: " << out_ec.apply_to_object(fst_obj);
        LOG_(error) << "Sec obj: " << out_ec.apply_to_object(sec_obj);
    }

    return out_ec;
}

#endif  // INCLUDE_LBBUILDER_H_
