/* Copyright 2017 Baytekov Nikita */

//=================================
// MONS Classifier test
//---------------------------------
// Testing MONSClassifier's quality
//=================================

#include <iostream>
#include "log_wrapper.h"
#include "MONSClassifier.h"
#include "DatasetManager.h"

#include "SampleSet.h"
#include "SampleHandler.h"

#include "ml_metrics.h"


using namespace ml_metrics;

int main(int argc, char* argv[]) {
    init_logging();

    DatasetManager<int> data_reader;
    Mat<int> in_data = data_reader.read_csv("datasets/audiology.csv");
    Vec<int> in_target = in_data.get_col(-1);
    in_data = in_data.get_rect(0, 0, -1, in_data.get_sy()-1);

    CrossValSH<int, int> cv_sh(10);
    Mat< SampleSet<int, int> > data_mat = cv_sh.make_samples(in_data, in_target, true, true);
    // int data_num = data_mat.get_sx();
    // LOG_(trace) << "Made sample sets successfully!";
    // LOG_(trace) << "Here:" << data_mat;

    GeneticDualizer<int, bool> gen_dual(3, 0.01, true, kDataScore, kDelayConverged, 0.005);
    MONSClassifier<int, bool> mons(gen_dual);

    Mat<int> dummy_data;
    Vec<int> dummy_target;
    Vec<int> dummy_answer;

    Vec<float> quality(10);

    for (int i = 0; i < 10; i++) {
        data_mat[i][0].get_data(&dummy_data, &dummy_target);
        mons.fit(dummy_data, dummy_target);

        mons.print_colls();

        data_mat[i][1].get_data(&dummy_data, &dummy_target);
        dummy_answer = mons.predict(dummy_data);
        LOG_(trace) << "Answer:" << dummy_answer;
        LOG_(trace) << "Real:" << dummy_target;
        quality[i] = ml_accuracy<int>(dummy_answer, dummy_target);
        LOG_(info) << "Accuracy " << i << ": " << quality[i];
    }
    float avg_quality = 0.0;
    for (int i = 0; i < 10; i++)
        avg_quality += quality[i];
    LOG_(info) << "Quality:" << quality;
    LOG_(info) << "Avg quality:" << avg_quality/10;
    // LOG_(info) << "MONS was fitted successfully.";
    // mons.save_mons_data("st2.mcl");
    // mons.load_mons_data("st2.mcl");
    //
    // Vec<int> answ_target = mons.predict(test_data);
    // LOG_(info) << "Answer was predicted successfully...";
    // LOG_(debug) << "My answer: " << answ_target;
    // LOG_(debug) << "True answer: " << test_target;
    // LOG_(info) << "Metric: " << ml_accuracy<int>(answ_target, test_target);
    // LOG_(debug) << "Program ended successfully.";
    return 0;
}
