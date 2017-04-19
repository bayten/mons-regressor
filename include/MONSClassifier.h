/* Copyright 2017 Baytekov Nikita */

#include "../include/default_types.h"
#include "../include/CollFamily.h"
#include "../include/SampleHandler.h"
#include "../include/LBBuilder.h"
#include "../include/GeneticDualizer.h"

#ifndef INCLUDE_MONSCLASSIFIER_H_
#define INCLUDE_MONSCLASSIFIER_H_

template<typename S, typename T>
class MONSClassifier {
    GeneticDualizer<S, T> genetic_dualizer;
    LBBuilder<S> lb_builder;
    SampleHandler<S> sample_handler;

    Vec< CollFamily<S> > coll_sets;
    SampleSet<S> train_set;

    int max_iter;
    float eps;

 public:
    explicit MONSClassifier(GeneticDualizer<S, T> init_gen_dual,
                            LBBuilder<S> init_lb_builder = RandomLBBuilder<S>(),
                            int init_miter = 1000, float init_eps = 0.0001);
    ~MONSClassifier() {}

    void fit(const Mat<S>& X, const Vec<int>& y);
    Vec<int> predict(const Mat<S>& X);

 private:
    bool check_margin(const SampleSet<S>& valid, const SampleSet<S>& train,
                      Vec<float>* avg_margins, int curr_iter);
    float get_class_estim(const SampleSet<S>& train, const Vec<S>& x, int class_tag);
    Vec<float> get_class_estim(const SampleSet<S>& train, const Mat<S>& X, int class_tag);
};

#endif  // INCLUDE_MONSCLASSIFIER_H_
