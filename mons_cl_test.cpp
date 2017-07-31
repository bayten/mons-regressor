/* Copyright 2017 Baytekov Nikita */

#include <iostream>
#include "log_wrapper.h"
#include "MONSClassifier.h"
#include "DatasetManager.h"
#include "ml_metrics.h"

using namespace ml_metrics;

int main(int argc, char* argv[]) {
    init_logging();

    DatasetManager<int> data_reader;
    Mat<int> data = data_reader.read_csv("datasets/audiology.csv");
    Vec<int> target = data.get_col(-1);

    data = data.get_rect(0, 0, -1, data.get_sy()-1);
    LOG_(info) << "Data matrix: " << data;
    LOG_(info) << "Target vector: " << target;

    const float train_frac = 0.9;

    Mat<int> train_data = data.get_rect(0, 0, data.get_sx()*train_frac, -1);
    Vec<int> train_target = target.slice(0, target.get_size()*train_frac);
    Mat<int> test_data = data.get_rect(data.get_sx()*train_frac, 0, -1, -1);
    Vec<int> test_target = target.slice(target.get_size()*train_frac, -1);

    LOG_(info) << "TRAIN Data matrix: " << train_data;
    LOG_(info) << "TRAIN Target vector: " << train_target;


    LOG_(info) << "TEST Data matrix: " << test_data;
    LOG_(info) << "TEST Target vector: " << test_target;
    GeneticDualizer<int, bool> gen_dual(5, 0.1, true, kDataScore, kDelayConverged, 0.001);
    MONSClassifier<int, bool> mons(gen_dual);
    mons.fit(train_data, train_target);

    LOG_(info) << "MONS was fitted successfully.";
    Vec<int> answ_target = mons.predict(test_data);
    LOG_(info) << "Answer was predicted successfully...";
    LOG_(debug) << "My answer: " << answ_target;
    LOG_(debug) << "True answer: " << test_target;
    LOG_(info) << "Metric: " << ml_accuracy<int>(answ_target, test_target);
    LOG_(debug) << "Program ended successfully.";
    return 0;
}
