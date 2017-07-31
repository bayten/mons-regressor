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

Vec<float> mean_filter(Vec<float> raw_data, int x) {  // Savitsky-Golay filter - smmothing curve
    int data_size = raw_data.get_size();
    float normaliser = 2.0*static_cast<float>(x)+1.0;
    Vec<float> smoothed(data_size);
    for (int i = 0; i < data_size; i++) {
        smoothed[i] = 0;
        for (int j = -x; j <= x; j++) {
            smoothed[i] += raw_data[std::min(std::max(0, i+j), data_size-1)];
        }
        smoothed[i] = smoothed[i] / normaliser;
    }
    return smoothed;
}

template<typename S>
Mat<S> norm_data(const Mat<S>& data) {
    int sx = data.get_sx(), sy = data.get_sy();

    Vec<float> norm_vals = data[0];
    Vec<float> min_vals = data[0];

    for (int i = 0; i < sy; i++) {
        for (int j = 1; j < sx; j++) {
            if (data[j][i] > norm_vals[i])
                norm_vals[i] = data[j][i];
            if (data[j][i] < min_vals[i])
                min_vals[i] = data[j][i];
        }
        norm_vals[i] = norm_vals[i] - min_vals[i];
    }

    Mat<S> normed_data = data;
    for (int i = 0; i < sx; i++)
        for (int j = 0; j < sy; j++)
            normed_data[i][j] = (data[i][j] - min_vals[j])/norm_vals[j];
    return normed_data;
}

template<typename S>
Vec<int> dmdbscan(const Mat<S>& data, MetricType used_metric = kEuclidean,
                  int min_neigh = 7) {
    int samples_num = data.get_sx();
    int k_param = min_neigh * 0.85;

    Mat<S> normed_data = norm_data(data);
    Mat<float> dist_mat(samples_num, samples_num);

    for (int i = 0; i < samples_num; i++) {  // computing distance matrix
        dist_mat[i][i] = -1;
        for (int j = i+1; j < samples_num; j++) {
            float dist = find_dist(normed_data[i], normed_data[j], used_metric);
            dist_mat[i][j] = dist;
            dist_mat[j][i] = dist;
        }
    }

    Vec<int> target_class(samples_num);
    Vec<int> was_processed(samples_num);  // 0 - wasn't processed, 1 - core, 2 - border/core

    Vec<float> kdist_data(samples_num);
    for (int i = 0; i < samples_num; i++) {  // filling k-dist vector
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
        kdist_data[i] = dist_mat[i][idxs[sort_idxs[k_param]]];
    }

    Vec<float> kdist_data_sorted = kdist_data.sort();
    Vec<float> kdist_smoothed = mean_filter(kdist_data_sorted, min_neigh);
    Vec<float> kdist_fst_deriv = get_diff_deriv(kdist_smoothed);
    Vec<float> kdist_sec_deriv = get_diff_deriv(kdist_fst_deriv);
    LOG_(trace) << "DATA SORTED: " << kdist_data_sorted;
    LOG_(trace) << "DATA SMOOTHED: " << kdist_smoothed;
    LOG_(trace) << "FST DERIV: " << kdist_fst_deriv;
    LOG_(trace) << "SEC DERIV: " << kdist_sec_deriv;

    Vec<float> epsilons;
    Vec<float> noise_borders;
    for (int i = 0; i < samples_num-1; i++) {  // estimating epsilons for different densities
        if (kdist_sec_deriv[i] > 0 && kdist_sec_deriv[i+1] < 0 &&
            kdist_sec_deriv[i]-kdist_sec_deriv[i+1] > 0.0001) {
            LOG_(trace) << "EPSILON SITUATION(S):" << kdist_smoothed[i] << " --- " << kdist_smoothed[i+1];
            LOG_(trace) << "EPSILON SITUATION(2D):" << kdist_sec_deriv[i] << " --- " << kdist_sec_deriv[i+1];
            LOG_(trace) << "EPSILON SITUATION(C): appending " << 0.7*kdist_smoothed[i]+0.3*kdist_smoothed[i+1];
            epsilons.append(0.7*kdist_smoothed[i]+0.3*kdist_smoothed[i+1]);
        } else if (kdist_sec_deriv[i] < 0 && kdist_sec_deriv[i+1] > 0 &&
                   kdist_sec_deriv[i+1]-kdist_sec_deriv[i] > 0.0001) {
            LOG_(trace) << "NOISE SITUATION:" <<  kdist_smoothed[i] << " --- " << kdist_smoothed[i+1];
            LOG_(trace) << "NOISE SITUATION(C): appending " << kdist_smoothed[i];
            noise_borders.append(kdist_smoothed[i]);
        }
    }
    noise_borders.append(kdist_smoothed[-1]);
    LOG_(trace) << "Estimated epsilons: " << epsilons;
    LOG_(trace) << "Estimated noise borders: " << noise_borders;
    int epses_num = epsilons.get_size();
    Vec<float> matching_epses(samples_num);

    int noise_idx = 0, noise_sz = noise_borders.get_size();
    for (int i = 0; i < samples_num; i++) {  // acquring points to density levels
        if (kdist_data[i] >= epsilons[-1]) {  // condition of noise point
            target_class[i] = 0;
            was_processed[i] = 3;
        }
        noise_idx = 0;
        for (int j = 0; j < epses_num; j++) {
            if (kdist_data[i] < epsilons[j]) {
                matching_epses[i] = epsilons[j];
                break;
            }
            while (noise_borders[noise_idx] > epsilons[j] && noise_idx < noise_sz)
                noise_idx++;
            if (kdist_data[i] < noise_borders[noise_idx]) {
                target_class[i] = 0;
                was_processed[i] = 3;
                break;
            }
        }
    }

    int class_label = 1;
    for (int i = 0; i < samples_num; i++) {  // need to process core points only
        if (was_processed[i]) {
            LOG_(trace) << "Object " << i << " was processed already";
            continue;
        }
        LOG_(trace) << "Processing object " << i;

        int my_neighs = 0;
        for (int j = 0; j < samples_num; j++) {  // detecting, if current point is a core one
            if (dist_mat[i][j] <= matching_epses[i])
                my_neighs++;
        }

        if (my_neighs <= min_neigh) {
            LOG_(trace) << "Object " << i << " is a border point...";
            target_class[i] = -1;
            was_processed[i] = 2;
            continue;
        }

        target_class[i] = class_label;  // it is a core point, continue
        for (int j = 0; j < samples_num; j++) {  // matching all points from eps
            if (dist_mat[i][j] <= matching_epses[i])
                target_class[j] = class_label;
        }

        bool was_changed = 1;  // indicator of processing other points
        while (was_changed) {
            LOG_(trace) << "TARGET" << target_class;
            was_changed = 0;
            for (int j = 0; j < samples_num; j++) {
                if (target_class[j] == class_label && !was_processed[j]) {  // just border or core?
                    LOG_(trace) << "Inner processing object " << j << "...";
                    my_neighs = 0;
                    for (int k = 0; k < samples_num; k++)
                        if (dist_mat[j][k] <= matching_epses[j]) {
                            LOG_(trace) << "Object " << k << " within eps radius!";
                            my_neighs++;
                        }

                    if (my_neighs >= min_neigh) {  // is a core point
                        was_processed[j] = 1;
                        for (int k = 0; k < samples_num; k++)
                            if (dist_mat[j][k] <= matching_epses[j]) {
                                target_class[k] = class_label;
                                if(!was_processed[k])
                                    was_changed = 1;
                            }
                    } else {  // is a border point
                        target_class[i] = -1;
                        was_processed[j] = 2; //matching
                    }
                }
            }
        }
        LOG_(trace) << "Exited processing loop";
        LOG_(trace) << "Was processed:" << was_processed;
        LOG_(trace) << "Labels:" << target_class;
        class_label++;
    }

    for (int i = 0; i < samples_num; i++) {  // processing border points
        if (was_processed[i] != 2)
            continue;

        float min_dist = -1.0;
        target_class[i] = -1;
        for (int j = 0; j < samples_num; j++) {  // matching all points from eps
            if (dist_mat[i][j] <= matching_epses[i] && was_processed[j] == 1) {
                if (target_class[i] == -1 || dist_mat[i][j] < min_dist) {
                    min_dist = dist_mat[i][j];
                    target_class[i] = target_class[j];
                }
            }
        }
    }

    return target_class;
}

}
#endif  // INCLUDE_GENETIC_TYPES_H_
