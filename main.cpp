/* Copyright 2017 Baytekov Nikita */

#include <iostream>
#include "log_wrapper.h"
#include "ml_metrics.h"
#include "MONSRegressor.h"
#include "DatasetManager.h"

using namespace ml_metrics;
using namespace cluster_algos;

int main(int argc, char* argv[]) {
    init_logging();

    // LOG_(trace) << "A trace severity message";
    // LOG_(debug) << "A debug severity message";
    // LOG_(info) << "An informational severity message";
    // LOG_(warning) << "A warning severity message";
    // LOG_(error) << "An error severity message";
    // LOG_(fatal) << "A fatal severity message";

    DatasetManager<float> data_reader;
    Mat<float> data = data_reader.read_csv("two-moons.csv");
    Vec<float> float_target = data.get_col(-1);
    int float_target_sz = float_target.get_size();
    Vec<int> target(float_target_sz);

    for (int i = 0; i < float_target_sz; i++)
        target[i] = static_cast<int>(float_target[i]);

    data = data.get_rect(0, 0, -1, data.get_sy()-1);
    LOG_(info) << "Data matrix: " << data;
    LOG_(info) << "Target vector: " << target;

    float train_frac = 0.8;
    Mat<float> train_data = data.get_rect(0, 0, data.get_sx()*train_frac, -1);
    Vec<int> train_target = target.slice(0, target.get_size()*train_frac);
    Mat<float> test_data = data.get_rect(data.get_sx()*train_frac+1, 0, -1, -1);
    Vec<int> test_target = target.slice(target.get_size()*train_frac+1, -1);


    LOG_(info) << "TRAIN Data matrix: " << train_data;
    LOG_(info) << "TRAIN Target vector: " << train_target;


    LOG_(info) << "TEST Data matrix: " << test_data;
    LOG_(info) << "TEST Target vector: " << test_target;
    Vec<float> ca_params(2, {cluster_algos::kEuclidean, 7});
    GeneticDualizer<float, bool> gen_dual(5, 0.1, true, kDataScore, kDelayConverged, 0.001);
    MONSRegressor<float, int, bool> mons(gen_dual, new RandomLBBuilder<float>(), kDMDBSCAN, ca_params);
    mons.fit(train_data, train_target);

    LOG_(info) << "MONS was fitted successfully.";
    Vec<int> answ_target = mons.predict(test_data);
    LOG_(info) << "Answer was predicted successfully...";
    LOG_(info) << "Metric: " << ml_accuracy<int>(answ_target, test_target);
    LOG_(debug) << "Program ended successfully.";
    return 0;
}
