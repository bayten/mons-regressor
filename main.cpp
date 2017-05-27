/* Copyright 2017 Baytekov Nikita */

#include <iostream>
#include "log_wrapper.h"
#include "MONSClassifier.h"
#include "DatasetReader.h"

int main(int argc, char* argv[]) {
    init_logging();

    // LOG_(trace) << "A trace severity message";
    // LOG_(debug) << "A debug severity message";
    // LOG_(info) << "An informational severity message";
    // LOG_(warning) << "A warning severity message";
    // LOG_(error) << "An error severity message";
    // LOG_(fatal) << "A fatal severity message";

    DatasetReader<int> data_reader;
    Mat<int> data = data_reader.read_csv("datasets/classification/ttt.csv");
    Vec<int> target = data.get_col(-1);
    data = data.get_rect(0, 0, -1, data.get_sy()-1);

    LOG_(info) << "Data matrix: " << data;
    LOG_(info) << "Target vector: " << target;

    MONSClassifier<int, bool> mons_classifier(GeneticDualizer<int, bool>(3, 0.1));
    mons_classifier.fit(data, target);

    LOG_(debug) << "Program ended successfully.";
    return 0;
}
