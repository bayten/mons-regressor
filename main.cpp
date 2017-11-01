/* Copyright 2017 Baytekov Nikita */

#include <iostream>
#include "log_wrapper.h"
#include "ml_metrics.h"

#include "MONSRegressor.h"
#include "DatasetManager.h"

#include "SampleSet.h"
#include "SampleHandler.h"

using namespace ml_metrics;
using namespace cluster_algos;

int main(int argc, char* argv[]) {
    init_logging();

    DatasetManager<float> data_reader;
    char path[1000] = "datasets/09_cancer_wpbc.csv";
    if (argc > 1) {
        strcpy(path, argv[1]);
    }
    Mat<float> in_data = data_reader.read_csv(path);
    Vec<float> in_target = in_data.get_col(-1);
    in_data = in_data.get_rect(0, 0, -1, in_data.get_sy()-1);

    CrossValSH<float, float> cv_sh(10);
    Mat< SampleSet<float, float> > data_mat = cv_sh.make_samples(in_data, in_target, true, false);
    // int data_num = data_mat.get_sx();
    // LOG_(trace) << "Made sample sets successfully!";
    // LOG_(trace) << "Here:" << data_mat;
    Vec<float> ca_params(3, {cluster_algos::kEuclidean, 3, 3});
    GeneticDualizer<float, bool> gen_dual(3, 0.01, true, kDataScore, kDelayConverged, 0.0005);
    MONSRegressor<float, float, bool> mons(gen_dual, kDMDBSCAN, ca_params, \
                                           new RandomLBBuilder<float>(),
                                           1000, 0.0001, kGenoType, kNearestNeighboors);

    Mat<float> dummy_data;
    Vec<float> dummy_target;
    Vec<float> dummy_answer;

    Vec<float> quality(10);

    for (int i = 0; i < 10; i++) {
        data_mat[i][0].get_data(&dummy_data, &dummy_target);
        mons.fit(dummy_data, dummy_target);
        data_mat[i][1].get_data(&dummy_data, &dummy_target);
        LOG_(trace) << "Going to predict everything";
        dummy_answer = mons.predict(dummy_data);
        LOG_(trace) << "Answer:" << dummy_answer;
        LOG_(trace) << "Real:" << dummy_target;
        quality[i] = ml_mse<float>(dummy_answer, dummy_target);
        LOG_(info) << "MSE " << i << ": " << quality[i];
    }
    float avg_quality = 0.0;
    for (int i = 0; i < 10; i++)
        avg_quality += quality[i];
    LOG_(info) << "Quality:" << quality;
    LOG_(info) << "Avg quality:" << avg_quality/10;
    std::cout << "Quality for folds:" << quality;
    std::cout << "Average quality:" << avg_quality/10;
    
    return 0;
}
