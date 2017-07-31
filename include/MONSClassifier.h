/* Copyright 2017 Baytekov Nikita */

#include "default_types.h"
#include "CollFamily.h"
#include "SampleHandler.h"
#include "LBBuilder.h"
#include "GeneticDualizer.h"

#ifndef INCLUDE_MONSCLASSIFIER_H_
#define INCLUDE_MONSCLASSIFIER_H_

template<typename S, typename T>  // S - features, T - coverings
class MONSClassifier {
    GeneticDualizer<S, T> genetic_dualizer;
    LBBuilder<S, int>* lb_builder;
    SampleHandler<S, int> sample_handler;

    Vec< CollFamily<S> > coll_sets;
    Vec<int> class_tags;
    SampleSet<S, int> train_set;

    int max_iter;
    float eps;

 public:
    explicit MONSClassifier(GeneticDualizer<S, T> init_gen_dual,
                            LBBuilder<S, int>* init_lb_builder = new RandomLBBuilder<S>(),
                            int init_miter = 1000, float init_eps = 0.001);
    ~MONSClassifier();

    void fit(const Mat<S>& X, const Vec<int>& y);
    Vec<int> predict(const Mat<S>& X);

    bool save_coll_sets(const char path[]) const;
    bool load_coll_sets(const char path[]);

 private:
    bool check_margin(const SampleSet<S, int>& valid, const SampleSet<S, int>& train,
                      Vec<float>* avg_margins, int curr_iter);
    float get_class_estim(const SampleSet<S, int>& train, const Vec<S>& x, int class_tag);
    Vec<float> get_class_estim(const SampleSet<S, int>& train, const Mat<S>& X, int class_tag);
};



template<typename S, typename T>
MONSClassifier<S, T>::MONSClassifier(GeneticDualizer<S, T> init_gen_dual,
                                     LBBuilder<S, int>* init_lb_builder,
                                     int init_miter, float init_eps):
        genetic_dualizer(init_gen_dual),
        lb_builder(init_lb_builder),
        sample_handler(),
        coll_sets(),
        class_tags(),
        max_iter(init_miter),
        eps(init_eps) {
    srand(time(0));
    LOG_(debug) << "MONSClassifier instance was created (by usual constructor).";
}

template<typename S, typename T>
MONSClassifier<S, T>::~MONSClassifier() {
    delete lb_builder;
    LOG_(debug) << "MONSClassifier was destroyed.";
}

template<typename S, typename T>
void MONSClassifier<S, T>::fit(const Mat<S>& X, const Vec<int>& y) {
    LOG_(trace) << "Process of model fitting has begun...";
    train_set = sample_handler.make_samples(X, y, true);
    LOG_(trace) << "Initial data set:" << train_set;

    SampleSet<S, int> train, valid;
    sample_handler.make_train_and_valid(train_set, &train, &valid);
    int class_num = train_set.get_group_num();
    coll_sets = Vec<CollFamily<S>>(class_num);
    class_tags = train_set.get_tags();

    int train_num = train.get_group_num();
    Vec<int> train_tags = train.get_tags();

    LOG_(trace) << "Train set: " << train;
    LOG_(trace) << "Valid set: " << valid;

    Vec<float> avg_margins(valid.get_total_size());
    for (int i = 0; i < valid.get_total_size(); i++)
        avg_margins[i] = 0.0;

    LOG_(trace) << "Starting main fitting loop(class num: " << class_num << ")...";

    for (int i = 1; i <= max_iter; i++) {
        LOG_(trace) << "Loop iteration: " << i;
        for (int k = 0; k < train_num; k++) {
            LOG_(trace) << "For class: " << train_tags[k];
            ElColl<S> local_basis = lb_builder->build_lb(train, train_tags[k]);
            LOG_(trace) << "Local basis: " << local_basis;

            LOG_(trace) << "Processing Genetic Algorithm...";
            genetic_dualizer.set_init_data(train, local_basis, train_tags[k]);
            Vec< Vec<T> > encoded_colls = genetic_dualizer.execute_ga();
            Vec< ElColl<S> > new_colls = genetic_dualizer.decode_collections(encoded_colls);

            LOG_(trace) << "Adding new collections...:" << new_colls;
            coll_sets[k].add(new_colls);
            LOG_(trace) << "New collections were successfully added.";
        }
        if (check_margin(valid, train, &avg_margins, i)) {
            LOG_(trace) << "Main loop was stopped because of margins' condition";
            break;
        }
    }
}

template<typename S, typename T>
Vec<int> MONSClassifier<S, T>::predict(const Mat<S> & X) {
    LOG_(trace) << "Predicting stuff! Get prepared!";
    LOG_(trace) << "Collection sets: " << coll_sets;
    int class_num = coll_sets.get_size();
    int obj_num = X.get_sx();
    Vec<int> class_preds(obj_num);
    Vec<float> estim_vec = get_class_estim(train_set, X, class_tags[0]);
    for (int i = 0; i < obj_num; i++)
        class_preds[i] = 0;

    for (int i = 1; i < class_num; i++) {
        Vec<float> new_estim_vec = get_class_estim(train_set, X, class_tags[i]);
        for (int j = 0; j < obj_num; j++) {
            if (new_estim_vec[j] > estim_vec[j]) {
                estim_vec[j] = new_estim_vec[j];
                class_preds[j] = i;
            }
        }
    }

    return class_preds;
}

template<typename S, typename T>
bool MONSClassifier<S, T>::check_margin(const SampleSet<S, int>& valid,
                                        const SampleSet<S, int>& train,
                                        Vec<float>* avg_margins,
                                        int curr_iter) {
    LOG_(trace) << "Checking margin...";
    LOG_(trace) << "Prev margins:" << *avg_margins;

    int valid_num = valid.get_group_num();
    Vec<int> valid_tags = valid.get_tags();

    int train_num = train.get_group_num();
    Vec<int> train_tags = train.get_tags();

    int total_size = valid.get_total_size();
    bool was_satisfied = 1;
    int curr_idx = 0;
    float curr_margin = 0.0;
    Vec<float> new_margins(total_size);
    Vec<float> delta;

    for (int i = 0; i < valid_num; i++) {
        LOG_(trace) << "Checking class " << valid_tags[i] << "...";
        const GroupSamples<S, int>& valid_samples = valid.get_group(valid_tags[i]);
        int valid_samples_size = valid_samples.get_size();
        for (int j = 0; j < valid_samples_size; j++) {
            float my_estim = 0.0, max_estim = 0.0;
            for (int k = 0; k < train_num; k++) {
                float k_estim = get_class_estim(train, valid_samples[j], train_tags[k]);
                if (class_tags[k] == class_tags[i])
                    my_estim = k_estim;
                else if (max_estim < k_estim)
                    max_estim = k_estim;
            }
            curr_margin = my_estim - max_estim;
            new_margins[curr_idx] = ((*avg_margins)[curr_idx]*(curr_iter-1)+curr_margin)/curr_iter;
            delta.append(new_margins[curr_idx] - (*avg_margins)[curr_idx]);
            if (new_margins[curr_idx] - (*avg_margins)[curr_idx] >= eps)
                was_satisfied = 0;
            curr_idx++;
        }
    }
    LOG_(trace) << "Delta vec:  " << delta;
    LOG_(trace) << "New margins:" << new_margins;
    *avg_margins = new_margins;

    return was_satisfied;
}

template<typename S, typename T>
float MONSClassifier<S, T>::get_class_estim(const SampleSet<S, int>& train,
                                            const Vec<S>& x, int class_tag) {
    const GroupSamples<S, int>& needed_class = train.get_group(class_tag);
    CollFamily<S>& needed_family = coll_sets[class_tags.where(class_tag)];
    int class_size = needed_class.get_size();
    int coll_num = needed_family.get_size();

    float estim_sum = 0.0;
    for (int i = 0; i < class_size; i++)
        for (int j = 0; j < coll_num; j++)
            estim_sum += needed_family[j].vote_func(x, needed_class[i]);
    return estim_sum/class_size/coll_num;
}

template<typename S, typename T>
Vec<float> MONSClassifier<S, T>::get_class_estim(const SampleSet<S, int>& train,
                                                 const Mat<S>& X, int class_tag) {
    LOG_(trace) << "Launching class estim for tag:" << class_tag;
    LOG_(trace) << "Class tags:" << class_tags;

    const GroupSamples<S, int>& needed_class = train.get_group(class_tag);
    if (class_tags.where(class_tag) == -1) {
        LOG_(error) << "Unknown class tag!";
        return Vec<float>();
    }

    CollFamily<S>& needed_family = coll_sets[class_tags.where(class_tag)];
    int class_size = needed_class.get_size();
    int coll_num = needed_family.get_size();
    int mat_sx = X.get_sx();

    Vec<float> estim_vec_sum(mat_sx);
    for (int i = 0; i < class_size; i++)
        for (int j = 0; j < coll_num; j++)
            for (int k = 0; k < mat_sx; k++)
                estim_vec_sum[k] += needed_family[j].vote_func(X[k], needed_class[i]);

    float divisor = class_size * coll_num;
    for (int i = 0; i < mat_sx; i++)
        estim_vec_sum[i] /= divisor;
    return estim_vec_sum;
}


#endif  // INCLUDE_MONSCLASSIFIER_H_
