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
    // kCustomType = 3  - TODO: make custom selection of features by bool vec
};
enum RegressionCastType {
    kAverage = 0,  // also median, etc...
    kNearestNeighboors = 1
    // kRecursive = 2
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

    SampleSet<int, T> cluster_matches;
 public:
    MONSRegressor(GeneticDualizer<S, U> init_gen_dual,
                  LBBuilder<S, int>* init_lb_builder,
                  ClusterAlgoType init_ca, Vec<float> init_ca_params,
                  int init_miter = 1000, float init_eps = 0.0001,
                  RegressionFeatureType init_rf = kGenoType, RegressionCastType init_rc = kAverage);
     ~MONSRegressor() {};

    void fit(const Mat<S>& X, const Vec<T>& y);
    Vec<T> predict(const Mat<S>& X);

 private:
    Mat<S> get_feature_data(const Mat<S>& X, const Vec<T>& y);
};


template<typename S, typename T, typename U>
MONSRegressor<S, T, U>::MONSRegressor(GeneticDualizer<S, U> init_gen_dual,
                                      LBBuilder<S, int>* init_lb_builder,
                                      ClusterAlgoType init_ca, Vec<float> init_ca_params,
                                      int init_miter, float init_eps, RegressionFeatureType init_rf,
                                      RegressionCastType init_rc) :
                    mons_instance(init_gen_dual, init_lb_builder, init_miter, init_eps),
                    reg_features(init_rf), reg_cast(init_rc),
                    cluster_algorithm(init_ca), ca_params(init_ca_params) {
}

template<typename S, typename T, typename U>
void MONSRegressor<S, T, U>::fit(const Mat<S>& X, const Vec<T>& y) {
    Mat<S> f_data = get_feature_data(X, y);
    Vec<int> cluster_y;
    switch (cluster_algorithm) {
        case kDMDBSCAN:
            cluster_y = cluster_algos::dmdbscan<S>(f_data, cluster_algos::MetricType(ca_params[0]), int(ca_params[1]) );
            break;

        default:
            LOG_(error) << "Unknown Cluster Algorithm Type(code:" << cluster_algorithm << ")";
    }
    cluster_matches = SampleSet<T, int>(y, cluster_y);
    mons_instance.fit(X, cluster_y);
}

template<typename S, typename T, typename U>
Vec<T> MONSRegressor<S, T, U>::predict(const Mat<S>& X) {
    int obj_num = X.get_sx();
    Vec<int> cluster_pred = mons_instance.predict(X);
    Vec<T> pred_vec(X.get_sx());

    switch (reg_cast) {
        case kAverage: {
            for (int i = 0; i < obj_num; i++) {
                GroupSamples<int, T> my_group = cluster_matches.get_group(cluster_pred[i]);
                int my_group_size = my_group.get_size();
                T group_sum = my_group[0][0];
                for (int j = 1; j < my_group_size; j++)
                    group_sum = group_sum + my_group[i][0];
                pred_vec[i] = group_sum / my_group_size;
            }
            break;
        }
        case kNearestNeighboors: {
            // TODO: Realise kNN
            break;
        }
        default:
            LOG_(error) << "Unknown Regression Cast Type(code:" << reg_cast << ")";
            break;
    }

    return pred_vec;
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
