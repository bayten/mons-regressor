/* Copyright 2017 Baytekov Nikita */
#ifndef INCLUDE_CONTAINERS_H_
#define INCLUDE_CONTAINERS_H_

#include <algorithm>
#include "./include/default_types.h"

template<class T>  // type of features
struct ElClass {
 public:
    Vec<int> cols;
    Vec<T> vals;

    explicit ElClass(int rank):
        cols(rank), vals(rank) {
        for (int i = 0; i < rank; i++)
            cols[i] = -1;
    }

    ElClass(int rank, int init_cols[], T init_vals[]) : cols(rank), vals(rank) {
        for (int i = 0; i < rank; i++) {
            cols[i] = init_cols[i];
            vals[i] = init_vals[i];
        }
    }

    ~ElClass() {
    }
};


template<class T>
struct ElColl {
    Vec< ElClass<T> > ecs;

    explicit ElColl(int size) : ecs(size) {
    }

    ElColl(int size, ElClass<T> init_ecs[]) : ecs(size) {
        for (int i = 0; i < size; i++)
            ecs[i] = init_ecs[i];
    }

    ~ElColl() {
    }
};


template<class T>
struct CollFamily {
    Vec< ElColl<T> > colls;

    CollFamily(int size, ElColl<T> init_colls[]) : colls(size) {
        for (int i = 0; i < size; i++)
            colls[i] = init_colls[i];
    }

    ~CollFamily() {
    }

    CollFamily<T>& operator+(const CollFamily<T>& W1);  // operator of uniting two families
};

#endif  // INCLUDE_CONTAINERS_H_
