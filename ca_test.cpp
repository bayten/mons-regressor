/* Copyright 2017 Baytekov Nikita */

#include <iostream>
#include "log_wrapper.h"
#include "clustering_algorithms.h"
#include "DatasetReader.h"

using namespace cluster_algos;

int main(int argc, char* argv[]) {
    init_logging();

    DatasetReader<int> data_reader;
    Mat<int> data = data_reader.read_csv("datasets/classification/clustering.csv");
    Vec<int> target = data.get_col(-1);
    data = data.get_rect(0, 0, -1, data.get_sy()-1);

    LOG_(info) << "Data matrix: " << data;
    LOG_(info) << "Target vector: " << target;

    SampleSet<int, int> sset;
    sset.append(data, target);

    SampleSet<int, int> clusters = dbscan(sset, kEuclidean, 2, 3);

    LOG_(info) << "Clustered:" << clusters;
    LOG_(debug) << "Program ended successfully.";
    return 0;
}
