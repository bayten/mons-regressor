/* Copyright 2017 Baytekov Nikita */
#include "./include/CollFamily.h"

template<class T>
ElClass<T>::ElClass(int rank):
        cols(rank), vals(rank) {
    for (int i = 0; i < rank; i++)
        cols[i] = -1;
}

template<class T>
ElClass<T>::ElClass(int rank, int init_cols[], T init_vals[]):
        cols(rank), vals(rank) {
    for (int i = 0; i < rank; i++) {
        cols[i] = init_cols[i];
        vals[i] = init_vals[i];
    }
}

template<class T>
bool ElClass<T>::operator==(const ElClass<T>& comp_obj) {
    int my_size = cols.get_size();
    if (my_size != comp_obj.cols.get_size())
        return 0;
    cols = cols.sort();
    vals = vals.sort();
    Vec<int> cols2 = comp_obj.cols.sort();
    Vec<T> vals2 = comp_obj.vals.sort();

    for (int i = 0; i < my_size; i++) {
        if (cols[i] != cols2[i] ||
            vals[i] != vals2[i])
            return 0;
    }
}


template<class T>
ElColl<T>::ElColl(int size, ElClass<T> init_ecs[]) : ecs(size) {
    for (int i = 0; i < size; i++)
        ecs[i] = init_ecs[i];
}

template<class T>
bool ElColl<T>::add(ElClass<T> ec_obj) {
    int ec_num = ecs.get_size();
    for (int i = 0; i < ec_num; i++)  // checking if there any similar elementary classifiers
        if (ecs[i] == ec_obj)
            return 0;
    ecs.push_back(ec_obj);
}
