/* Copyright 2017 Baytekov Nikita */
#include "./include/MONSClassifier.h"

template<typename S, typename T>
MONSClassifier<S, T>::MONSClassifier(GeneticDualizer<S, T> init_gen_dual,
                                     LBBuilder<S> init_lb_builder,
                                     int init_miter, float init_eps):
        genetic_dualizer(init_gen_dual),
        lb_builder(init_lb_builder),
        sample_handler(),
        coll_sets(),
        max_iter(init_miter),
        eps(init_eps) {
}

template<typename S, typename T>
void MONSClassifier<S, T>::fit(const Mat<S>& X, const Vec<int>& y) {
    train_set = sample_handler.make_samples(X, y, true);
    SampleSet<S> train, valid;
    sample_handler.make_train_and_valid(train_set, &train, &valid);
    int class_num = train_set.get_class_num();

    Vec<float> avg_margins(valid.get_total_size());
    for (int i = 1; i <= max_iter; i++) {
        for (int k = 0; k < class_num; k++) {
            ElColl<S> local_basis = lb_builder.build_lb(train, k);

            genetic_dualizer.set_init_data(train, local_basis, k);
            Vec< Vec<T> > encoded_colls = genetic_dualizer.execute_ga();
            Vec< ElColl<S> > new_colls = genetic_dualizer.decode_collections(encoded_colls);

            coll_sets[k].add(new_colls);
        }

        if (check_margin(valid, &avg_margins, i))
            break;
    }
}

template<typename S, typename T>
Vec<int> MONSClassifier<S, T>::predict(const Mat<S> & X) {
    int class_num = coll_sets.get_size();
    int obj_num = X.get_sx();
    Vec<int> class_preds(obj_num);
    Vec<float> estim_vec = get_class_estim(train_set, X, 0);
    for (int i = 0; i < obj_num; i++)
        class_preds[i] = 0;

    for (int i = 1; i < class_num; i++) {
        Vec<float> new_estim_vec = get_class_estim(train_set, X, i);
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
bool MONSClassifier<S, T>::check_margin(const SampleSet<S>& valid,
                                        const SampleSet<S>& train,
                                        Vec<float>* avg_margins,
                                        int curr_iter) {
    int class_num = valid.get_class_num();
    int total_size = valid.get_total_size();
    bool was_satisfied = 1;
    int curr_idx = 0;
    float curr_margin = 0.0;
    Vec<float> new_margins(total_size);

    for (int i = 0; i < class_num; i++) {
        ClassSamples<S>& class_samples = valid.get_class(i);
        int class_samples_size = class_samples.get_size();
        for (int j = 0; j < class_samples_size; j++) {
            float my_estim = 0.0, max_estim = 0.0;
            for (int k = 0; k < class_num; k++) {
                if (k == i)
                    my_estim = get_class_estim(class_samples[j], k);
                else if (max_estim < get_class_estim(class_samples[j], k))
                    max_estim = get_class_estim(class_samples[j], k);
            }
            curr_margin = my_estim - max_estim;
            new_margins[curr_idx] = ((*avg_margins)[curr_idx]*(curr_iter-1)+curr_margin)/curr_iter;
            if (fabs(new_margins[curr_idx] - (*avg_margins)[curr_idx]) >= eps)
                was_satisfied = 0;
            curr_idx++;
        }
    }
    *avg_margins = new_margins;

    return was_satisfied;
}

template<typename S, typename T>
float MONSClassifier<S, T>::get_class_estim(const SampleSet<S>& train,
                                            const Vec<S>& x, int class_tag) {
    ClassSamples<S>& needed_class = train.get_class(class_tag);
    CollFamily<S>& needed_family = coll_sets[class_tag];
    int class_size = needed_class.get_size();
    int coll_num = needed_family.get_size();

    float estim_sum = 0.0;
    for (int i = 0; i < class_size; i++)
        for (int j = 0; j < coll_num; j++)
            estim_sum += needed_family.vote_func(x, needed_class[i]);
    return estim_sum/class_size/coll_num;
}

template<typename S, typename T>
Vec<float> MONSClassifier<S, T>::get_class_estim(const SampleSet<S>& train,
                                                 const Mat<S>& X, int class_tag) {
    ClassSamples<S>& needed_class = train.get_class(class_tag);
    CollFamily<S>& needed_family = coll_sets[class_tag];
    int class_size = needed_class.get_size();
    int coll_num = needed_family.get_size();
    int mat_sx = X.get_sx();

    Vec<float> estim_vec_sum;
    for (int i = 0; i < class_size; i++)
        for (int j = 0; j < coll_num; j++)
            for (int k = 0; k < mat_sx; k++)
                estim_vec_sum[k] += needed_family.vote_func(X[k], needed_class[i]);

    float divisor = class_size * coll_num;
    for (int i = 0; i < mat_sx; i++)
        estim_vec_sum[i] /= divisor;
    return estim_vec_sum;
}
