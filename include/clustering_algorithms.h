/* Copyright 2017 Baytekov Nikita */
#include "default_types.h"
#include "log_wrapper.h"
#include "SampleSet.h"

#include <cmath>

#ifndef INCLUDE_CLUSTERING_ALGORITHMS_H_
#define INCLUDE_CLUSTERING_ALGORITHMS_H_

namespace cluster_algos{

enum MetricType {
    kManhattan = 0,
    kEuclidean = 1
    // kCosine = 2
};

template<typename S>
float find_dist(const Vec<S>& fst_obj, const Vec<S>& sec_obj, MetricType used_metric) {
    // vectors will always have the same length (just not to check it every time)
    float out_dist = 0.0;
    int obj_len = fst_obj.get_size();

    switch (used_metric) {
        case kManhattan: {
            for (int i = 0; i < obj_len; i++)
                out_dist += fabs(fst_obj[i] - sec_obj[i]);
            break;
        }

        case kEuclidean: {
            for (int i = 0; i < obj_len; i++) {
                S diff = fst_obj[i]-sec_obj[i];
                out_dist += diff*diff;
            }
            out_dist = sqrt(out_dist);
            break;
        }

        default: {
            LOG_(error) << "Unknown metric type(" << used_metric << ")";
            return -1.0;
        }
    }

    return out_dist;
}

template<typename S>
Vec<int> dbscan(const Mat<S>& data, MetricType used_metric = kEuclidean,
                float eps = 0.1, int min_neigh = 3) {
    int samples_num = data.get_sx();
    Mat<float> dist_mat(samples_num, samples_num);

    for (int i = 0; i < samples_num; i++) {
        dist_mat[i][i] = -1;
        for (int j = i+1; j < samples_num; j++) {
            float dist = find_dist(data[i], data[j], used_metric);
            dist = (dist <= eps) ? dist : -1.0;
            dist_mat[i][j] = dist;
            dist_mat[j][i] = dist;
        }
    }
    LOG_(trace) << "DIST_MAT:" << dist_mat;
    Vec<int> target_class(samples_num);
    Vec<bool> was_processed(samples_num);
    for (int i = 0; i < samples_num; i++) {
        target_class[i] = -1;
        was_processed[i] = 0;
        int neighboors = 0;
        for (int j = 0; j < samples_num; j++)
            if (dist_mat[i][j] >= 0)
                neighboors++;
        if (neighboors < min_neigh) {
            target_class[i] = 0;
            was_processed[i] = 1;
        }
    }
    int class_label = 1;
    for (int i = 0; i < samples_num; i++) {
        if (was_processed[i]) {
            LOG_(trace) << "Object " << i << " was processed already";
            continue;
        }
        LOG_(trace) << "Processing object " << i;

        target_class[i] = class_label;
        bool was_changed = 1;
        while (was_changed) {
            LOG_(trace) << "TARGET" << target_class;
            was_changed = 0;
            for (int j = 0; j < samples_num; j++) {
                if (target_class[j] == class_label && !was_processed[j]) {
                    LOG_(trace) << "Inner processing object " << j << "...";
                    was_processed[j] = 1;
                    for (int k = 0; k < samples_num; k++)
                        if (dist_mat[j][k] >= 0) {
                            LOG_(trace) << "Object " << k << " within eps radius!";
                            target_class[k] = class_label;
                            if(!was_processed[k])
                                was_changed = 1;
                        }
                }
            }
        }
        LOG_(trace) << "Exited processing loop";
        class_label++;
    }
    return target_class;
}

Vec<float> get_diff_deriv(const Vec<float>& data, int h = 1) {
    int vec_len = data.get_size();
    Vec<float> out_vec(vec_len);

    for (int i = h; i < vec_len-h; i++)
        out_vec[i] = (data[i+h]- data[i-h])/static_cast<float>(2*h);

    for (int i = 0; i < h; i++) {
        out_vec[i] = out_vec[h];
        out_vec[vec_len-i-1] = out_vec[vec_len-h-1];
    }
    return out_vec;
}

template<typename S>
Vec<int> dmdbscan(const Mat<S>& data,
                           MetricType used_metric = kEuclidean,
                           int min_neigh = 3) {
   int samples_num = data.get_sx();
   Mat<float> dist_mat(samples_num, samples_num);

    for (int i = 0; i < samples_num; i++) {
        dist_mat[i][i] = -1;
        for (int j = i+1; j < samples_num; j++) {
            float dist = find_dist(data[i], data[j], used_metric);
            dist_mat[i][j] = dist;
            dist_mat[j][i] = dist;
        }
    }

    Vec<int> target_class(samples_num);
    Vec<bool> was_processed(samples_num);

    Vec<float> kdist_data(samples_num);
    for (int i = 0; i < samples_num; i++) {
        target_class[i] = -1;
        was_processed[i] = 0;

        Vec<float> vals;
        Vec<int> idxs;
        for (int j = 0; j < samples_num; j++) {
            if (dist_mat[i][j] >= 0.0) {
                vals.append(dist_mat[i][j]);
                idxs.append(j);
            }
        }
        Vec<int> sort_idxs = vals.sort_indices();
        kdist_data[i] = dist_mat[i][idxs[sort_idxs[min_neigh]]];
    }

    Vec<float> kdist_fst_deriv = get_diff_deriv(kdist_data.sort());
    Vec<float> kdist_sec_deriv = get_diff_deriv(kdist_fst_deriv);

    Vec<float> epsilons;

    for (int i = 0; i < samples_num-1; i++)
        if (kdist_sec_deriv[i] > 0 && kdist_sec_deriv[i+1] < 0)
            epsilons.append((kdist_data[i]+kdist_data[i+1])/2.0);
    LOG_(trace) << "Estimated epsilons: " << epsilons;

    int epses_num = epsilons.get_size();
    Vec<float> matching_epses(samples_num);

    for (int i = 0; i < samples_num; i++) {
        if (kdist_data[i] >= epsilons[-1]) {
            target_class[i] = 0;
            was_processed[i] = 1;
        }

        for (int j = 0; j < epses_num; j++) {
            if (kdist_data[i] < epsilons[j]) {
                matching_epses[i] = epsilons[j];
                break;
            }
        }
    }

    int class_label = 1;
    for (int i = 0; i < samples_num; i++) {
        if (was_processed[i]) {
            LOG_(trace) << "Object " << i << " was processed already";
            continue;
        }
        LOG_(trace) << "Processing object " << i;

        target_class[i] = class_label;
        bool was_changed = 1;
        while (was_changed) {
            LOG_(trace) << "TARGET" << target_class;
            was_changed = 0;
            for (int j = 0; j < samples_num; j++) {
                if (target_class[j] == class_label && !was_processed[j]) {
                    LOG_(trace) << "Inner processing object " << j << "...";
                    was_processed[j] = 1;
                    for (int k = 0; k < samples_num; k++)
                        if (dist_mat[j][k] <= matching_epses[j]) {
                            LOG_(trace) << "Object " << k << " within eps radius!";
                            target_class[k] = class_label;
                            if(!was_processed[k])
                                was_changed = 1;
                        }
                }
            }
        }
        LOG_(trace) << "Exited processing loop";
        class_label++;
    }

    return target_class;
}

}
#endif  // INCLUDE_GENETIC_TYPES_H_
