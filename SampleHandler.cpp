/* Copyright 2017 Baytekov Nikita */

#include "./include/SampleHandler.h"

template<class T>
SampleHandler<T>::SampleHandler() {
}

template<class T>
SampleHandler<T>::~SampleHandler() {
}

template<class T>
SampleSet<T> SampleHandler<T>::make_samples(const Mat<T>& X, const Vec<int>& y) {
    SampleSet<T> sample_set;

    sample_set.append(X, y);
    sample_set.shuffle();

    return sample_set;
}

template<class T>
int SampleHandler<T>::make_train_and_valid(const SampleSet<T>& sample_set, \
                                                SampleSet<T>* train, \
                                                SampleSet<T>* valid) {
    Mat<T> torn_X;
    Vec<int> torn_y;
    *train = sample_set;  // copying sample set to train
    int class_num = train->get_size();  // doing this to have samples of every
    for (int i = 0; i < class_num; i++) {  // class in validation set
        int tear_num = 2*static_cast<int>(log(train->get_class(i).get_size(), 10))+1;
        train->get_class(i).tearup(tear_num, &torn_X, &torn_y);
        valid->append(torn_X, torn_y);
    }

    return class_num;
}
