/* Copyright 2017 Baytekov Nikita */

#include "../include/SampleSet.h"

template<class T>
void ClassSamples<T>::tearup(int tear_num, Mat<T>* torn_x, Vec<int>* torn_y) {
    std::random_shuffle(&(objs[0]), &(objs[-1]));
    *torn_x = objs.cutslice(0, tear_num);
    *torn_y = Vec<int>(tear_num);
    for (int j = 0; j < tear_num; j++)
        (*torn_y)[j] = class_tag;
}

template<class T>
SampleSet<T>::SampleSet(Vec< ClassSamples<T> > init_samples) :
        samples(init_samples), total_size(0) {
    int class_num = samples.get_size();
    for (int i = 0; i < class_num; i++)
        total_size += samples[i].get_size();
}

template<class T>
void SampleSet<T>::append(Mat<T> X, Vec<int> y) {
    int obj_num = y.get_size();
    int samples_size = samples.get_size();
    int i = 0;

    if (!samples_size) {  // if appending for the first time...
        samples.append(ClassSamples<T>(X[0], y[0]));
        i++;
    }

    for (; i < obj_num; i++) {
        bool was_found = 0;

        for (int j = 0; j < samples_size; j++) {  // searching within existing classes
            if (samples[j].class_tag == y[i]) {
                samples[j].objs.append(X[i]);
                was_found = 1;
                break;
            }
        }
        if (!was_found) {  // otherwise adding new class container =)
            samples.append(ClassSamples<T>(X[i], y[i]));
        }
    }
}

template<class T>
void SampleSet<T>::shuffle() {
    int samples_size = samples.get_size();
    for (int i = 0; i < samples_size; i++) {
        std::random_shuffle(&(samples[i].obj[0]), &(samples[i].obj[-1]));
    }
}

template<class T>
ClassSamples<T>& SampleSet<T>::get_class(int index_tag) {
    int samples_size = samples.get_size();
    for (int i = 0; i < samples_size; i++)
        if (samples[i].class_tag == index_tag)
            return samples[i];

    return samples[0];
}

template<class T>
Vec<T>& SampleSet<T>::operator[](int abs_index) {
    int samples_num = samples.get_size();
    for (int i = 0; i < samples_num; i++) {
        if (abs_index < samples[i].get_size())
            return samples[i][abs_index];

        abs_index -= samples[i].get_size();
    }

    return samples[-1][-1];
}
