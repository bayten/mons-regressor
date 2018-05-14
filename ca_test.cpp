/* Copyright 2017 Baytekov Nikita */

//========================================
// Clustering algorithms test
//----------------------------------------
// Testing different clustering algorithms
//========================================

#include <iostream>
#include "log_wrapper.h"
#include "clustering_algorithms.h"
#include "DatasetManager.h"

using namespace cluster_algos;

int main(int argc, char* argv[]) {
    init_logging();

    DatasetManager<int> data_reader;
    Mat<int> data = data_reader.read_csv("datasets/02_yachts.csv");
    if(data.get_sx() <= 0 || data.get_sy() <= 0) {
        LOG_(error) << "Incorrect data format.";
        return 1;
    }

    Vec<int> target = data.get_col(-1);
    // data = data.get_rect(0, 0, -1, data.get_sy()-1);

    LOG_(info) << "Data matrix: " << data;
    LOG_(info) << "Target vector: " << target;

    Vec<int> clusters = vparzen(Mat<int>(target), 13);

    LOG_(info) << "Clustered:" << clusters;
    LOG_(debug) << "Program ended successfully.";
    return 0;
}
