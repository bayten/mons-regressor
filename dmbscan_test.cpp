/* Copyright 2017 Baytekov Nikita */

#include "default_types.h"
#include "log_wrapper.h"
#include "clustering_algorithms.h"
#include "DatasetManager.h"

using namespace cluster_algos;

int main(void) {
    init_logging();
    DatasetManager<float> data_manager;
    Mat<float> data = data_manager.read_csv("two-moons.csv");
    data = data.get_rect(0, 0, -1, data.get_sy()-1);
    Vec<int> clusters = dmdbscan(data, kEuclidean, 5);
    int size = clusters.get_size();

    Vec<float> float_clusters(size);
    for (int i = 0; i < size; i++)
        float_clusters[i] = static_cast<float>(clusters[i]);

    LOG_(trace) << "data:" << data;
    data_manager.write_csv("two-moons-out.csv", data.vadd(float_clusters));
    return 0;
}
