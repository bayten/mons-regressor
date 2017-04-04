/* Copyright 2017 Baytekov Nikita */
#include "./include/Population.h"

template<class T>
Individual<ElColl<T> > CoveringIndividual<T>::operator*(const Individual<ElColl<T> >& cross_obj) {
}

template<class T>
float CoveringIndividual<T>::get_quality(const SampleSet<T>& basic, \
                                         const SampleSet<T>& valid, \
                                         int target_tag) {
    ClassSamples<T>& basic_class = basic.get_class(target_tag);
    ClassSamples<T>& valid_class = valid.get_class(target_tag);

    int basic_class_num = basic_class.get_size();
    int valid_class_num = valid_class.get_size();

    float quality_sum = 0.0;
    for (int i = 0; i < basic_class_num; i++)
        for (int j = 0; j < valid_class_num; j++)
            quality_sum += this->data.corr_func(basic_class[i], valid_class[j]);

    return quality_sum / static_cast<float>(valid_class_num);
}

template<class T>
void Population<T>::update_costs(const SampleSet<T>& basic, \
                                 const SampleSet<T>& valid, \
                                 int target_tag) {
    int ind_num = ind_vec.get_size();
    float qualities[ind_num] = {};
    float min_val = -1.0;

    for (int i = 0; i < ind_num; i++) {
        qualities[i] = ind_vec[i].get_quality(basic, valid, target_tag);
        if (qualities[i] < min_val || min_val < 0.0)
            min_val = qualities[i];
    }

    for (int i = 0; i < ind_num; i++)
        ind_vec[i].set_cost(qualities[i]-min_val+1);
}
