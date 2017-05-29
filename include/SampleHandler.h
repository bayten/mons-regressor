/* Copyright 2017 Baytekov Nikita */
#include <type_traits>

#include "default_types.h"
#include "CollFamily.h"
#include "SampleSet.h"
#ifndef INCLUDE_SAMPLEHANDLER_H_
#define INCLUDE_SAMPLEHANDLER_H_

template<typename S, typename T>  // S - type of features, T - type of answers
class SampleHandler {
 public:
    SampleHandler() { LOG_(debug) << "Created SampleHandler instance successfully"; }
    ~SampleHandler() {}

    SampleSet<S, T> make_samples(const Mat<S>& X, const Vec<T>& y, bool mix = true);
    SampleSet<S, T> make_samples(const Vec<S>& X, const Vec<T>& y, bool mix = true);
    int make_train_and_valid(const SampleSet<S, T>& sample_set, \
                             SampleSet<S, T>* train, \
                             SampleSet<S, T>* valid);
};


template<typename S, typename T>
SampleSet<S, T> SampleHandler<S, T>::make_samples(const Mat<S>& X, const Vec<T>& y, bool mix) {
    LOG_(trace) << "Making samples...";
    SampleSet<S, T> sample_set;

    sample_set.append(X, y);
    LOG_(trace) << "Sample Set was made.";

    if (mix) {
        sample_set.shuffle();
        LOG_(trace) << "Samples were mixed";
    }
    return sample_set;
}

template<typename S, typename T>
SampleSet<S, T> SampleHandler<S, T>::make_samples(const Vec<S>& X, const Vec<T>& y, bool mix) {
    LOG_(trace) << "Making samples...";
    SampleSet<S, T> sample_set;
    int x_size = X.get_size();
    Mat<S> transposed_X(x_size, 1);
    for (int i = 0; i < x_size; i++)
        transposed_X[i][0] = X[i];
    sample_set.append(transposed_X, y);
    LOG_(trace) << "Sample Set was made.";

    if (mix) {
        sample_set.shuffle();
        LOG_(trace) << "Samples were mixed";
    }
    return sample_set;
}

template<typename S, typename T>
int SampleHandler<S, T>::make_train_and_valid(const SampleSet<S, T>& sample_set, \
                                                SampleSet<S, T>* train, \
                                                SampleSet<S, T>* valid) {
    LOG_(trace) << "Splitting sample set on train and valid sets...";
    LOG_(trace) << "Input sample set:" << sample_set;
    (*train) = SampleSet<S, T>();
    (*valid) = SampleSet<S, T>();

    Mat<S> slice_X;
    Vec<T> slice_y;

    int group_num = sample_set.get_group_num();  // doing this to have samples of every
    Vec<T> sset_tags = sample_set.get_tags();
    LOG_(trace) << "Group amount: " << group_num;


    if (std::is_same<T, int>::value) {  // classification case
        for (int i = 0; i < group_num; i++) {  // class in validation set
            const GroupSamples<S, T>& curr_group = sample_set.get_group(sset_tags[i]);
            int curr_group_size = curr_group.get_size();
            int slice_num = static_cast<int>(2 * log(curr_group_size)+1);

            if (2*slice_num >= curr_group_size) {  // this case can occur in case of small train sets
                int new_slice = std::max(static_cast<int>(0.4 * curr_group_size), 1);
                LOG_(trace) << "Using size " << new_slice << " instead of " << slice_num;
                slice_num = new_slice;
            }

            LOG_(trace) << "Taking " << slice_num << "/" << curr_group_size << ", group:" << sset_tags[i];
            curr_group.slice_rand(slice_num, &slice_X, &slice_y);

            valid->append(slice_X, slice_y);
            Vec<int> tag_vec(curr_group_size - slice_num);
            for (int j = 0; j < curr_group_size - slice_num; j++)
                tag_vec[j] = sset_tags[i];

            train->append(curr_group.get_objs().get_rect(slice_num, 0), tag_vec);
        }
    } else {
        // int needed_size = sample_set.get_total_size()
        // TODO: Implement
        LOG_(error) << "Not implemented yet.";
    }
    return group_num;
}

#endif  // INCLUDE_SAMPLEHANDLER_H_
