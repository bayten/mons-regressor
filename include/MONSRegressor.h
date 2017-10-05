/* Copyright 2017 Baytekov Nikita */

#include "default_types.h"
#include "MONSClassifier.h"
#include "SampleHandler.h"
#include "clustering_algorithms.h"
#ifndef INCLUDE_MONSREGRESSOR_H_
#define INCLUDE_MONSREGRESSOR_H_

enum RegressionFeatureType {
    kPhenoType  = 0,
    kGenoType   = 1,
    kGenoPhenoType = 2
};
enum RegressionCastType {
    kAverage = 0,
    kMedian = 1,
    kNearestNeighboors = 2
};
enum ClusterAlgoType {
    kDMDBSCAN = 0
    // TODO: Implement more algorithms and find out the best one
};

template<typename S, typename T, typename U>  // S - features, T - answers, U - coverings
class MONSRegressor {
    MONSClassifier<S, U> mons_instance;
    RegressionFeatureType reg_features;
    RegressionCastType reg_cast;
    ClusterAlgoType cluster_algorithm;
    Vec<float> ca_params;

    Mat<S> saved_data;
    Vec<T> saved_target;
    Vec<int> cluster_target;
 public:
    MONSRegressor(GeneticDualizer<S, U> init_gen_dual,
                  ClusterAlgoType init_ca, Vec<float> init_ca_params,
                  LBBuilder<S, int>* init_lb_builder = new RandomLBBuilder<S>(),
                  int init_miter = 1000, float init_eps = 0.0001,
                  RegressionFeatureType init_rf = kGenoType, RegressionCastType init_rc = kAverage);
     ~MONSRegressor() {};

    void fit(const Mat<S>& X, const Vec<T>& y);
    Vec<T> predict(const Mat<S>& X);

 private:
    Mat<S> get_feature_data(const Mat<S>& X, const Vec<T>& y);
    Vec<T> knn(const Mat<S>& X, const Vec<int>& cl_y, cluster_algos::MetricType metric, int k_param);
};


template<typename S, typename T, typename U>
MONSRegressor<S, T, U>::MONSRegressor(GeneticDualizer<S, U> init_gen_dual,
                                      ClusterAlgoType init_ca, Vec<float> init_ca_params,
                                      LBBuilder<S, int>* init_lb_builder,
                                      int init_miter, float init_eps, RegressionFeatureType init_rf,
                                      RegressionCastType init_rc) :
                    mons_instance(init_gen_dual, init_lb_builder, init_miter, init_eps),
                    reg_features(init_rf), reg_cast(init_rc),
                    cluster_algorithm(init_ca), ca_params(init_ca_params) {
}

template<typename S, typename T, typename U>
void MONSRegressor<S, T, U>::fit(const Mat<S>& X, const Vec<T>& y) {
    saved_data = get_feature_data(X, y);
    saved_target = y;
    switch (cluster_algorithm) {
        case kDMDBSCAN:
            cluster_target = cluster_algos::dmdbscan<S>(saved_data,
                                                        cluster_algos::MetricType(ca_params[0]),
                                                        int(ca_params[1]) );
            break;

        default:
            LOG_(error) << "Unknown Cluster Algorithm Type(code:" << cluster_algorithm << ")";
    }
    LOG_(trace) << "DBSCAN worked correctly: " << cluster_target;
    LOG_(trace) << "y:" << y;
    mons_instance.fit(X, cluster_target);
}

template<typename S, typename T, typename U>
Vec<T> MONSRegressor<S, T, U>::predict(const Mat<S>& X) {
    int obj_num = X.get_sx();
    int target_size = cluster_target.get_size();
    Vec<int> cluster_pred = mons_instance.predict(X);
    Vec<T> pred_vec(X.get_sx());
    LOG_(trace) << "pred_vec size:" << pred_vec.get_size();
    LOG_(trace) << "Predictions from MONSClassifier: " << cluster_pred;
    switch (reg_cast) {
        case kAverage: {
            LOG_(trace) << "kAverage activated";
            LOG_(trace) << "cluster_target:" << cluster_target;
            LOG_(trace) << "saved_target:" << saved_target;
            for (int i = 0; i < obj_num; i++) {
                T group_sum = static_cast<T>(0.0);
                float denominator = 0.0;
                for (int j = 0; j < target_size; j++) {
                    if (cluster_target[j] == cluster_pred[i]) {
                        group_sum += saved_target[j];
                        denominator += 1.0;
                    }
                }
                LOG_(trace) << "Denominator:" << denominator;
                pred_vec[i] = group_sum / denominator;
            }
            break;
        }
        case kMedian: {
            // TODO: Realise median
            break;
        }
        case kNearestNeighboors: {
            LOG_(trace) << "kNearestNeighboors activated.";
            pred_vec = knn(X, cluster_pred, cluster_algos::MetricType(ca_params[0]), int(ca_params[2]));
            break;
        }
        default:
            LOG_(error) << "Unknown Regression Cast Type(code:" << reg_cast << ")";
            break;
    }

    return pred_vec;
}

template<typename S, typename T, typename U>
Vec<T> MONSRegressor<S, T, U>::knn(const Mat<S>& X, const Vec<int>& cl_y,
                                   cluster_algos::MetricType metric, int k_param) {
    LOG_(trace) << "KNN was launched";
    int obj_num = cl_y.get_size();
    int saved_num = saved_target.get_size();
    if (k_param > saved_num)
        k_param = saved_num;

    LOG_(trace) << "cluster_y:" << cl_y;
    LOG_(trace) << "saved_target:" << saved_target;
    LOG_(trace) << "k_param:" << k_param;

    Vec<T> out_vec(obj_num);
    Vec<int> sort_idx;
    for (int i = 0; i < obj_num; i++) {
        Vec<float> dist_vec(saved_num);
        for (int j = 0; j < saved_num; j++)
            dist_vec[j] = cluster_algos::find_dist(X[i], saved_data[j], metric);
        LOG_(trace) << "Dist vec " << i << ":" << dist_vec;
        sort_idx = dist_vec.sort_indices();
        LOG_(trace) << "Indices for sort:" << sort_idx;
        out_vec[i] = static_cast<T>(0.0);
        for (int j = 0; j < k_param; j++) {
            LOG_(trace) << "Summation launched!";
            out_vec[i] += saved_target[sort_idx[j]];
        }
        LOG_(trace) << "Out_vec[" << i << "]:" << out_vec[i];
        out_vec[i] /= float(k_param);
    }
    return out_vec;
}
template<typename S, typename T, typename U>
Mat<S> MONSRegressor<S, T, U>::get_feature_data(const Mat<S>& X, const Vec<T>& y) {
    int sx = X.get_sx();

    switch (reg_features) {
        case kPhenoType: {
            Mat<S> out_mat(sx, 1);
            for (int i = 0; i < sx; i++)
                out_mat[i][0] = static_cast<S>(y[i]);
            return out_mat;
        }

        case kGenoType: {
            return X;
        }

        case kGenoPhenoType: {
            Mat<S> out_mat = X;
            Vec<S> casted_vec(sx);
            for (int i = 0; i < sx; i++)
                casted_vec[i] = static_cast<S>(y[i]);
            return out_mat.vadd(casted_vec);
        }

        default: {
            LOG_(error) << "Unknown Regression Cast Type(code:" << reg_features << ")";
            break;
        }
    }

    return Mat<S>();
}
#endif  // INCLUDE_MONSREGRESSOR_H_
