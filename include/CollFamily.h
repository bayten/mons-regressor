/* Copyright 2017 Baytekov Nikita */
#ifndef INCLUDE_COLLFAMILY_H_
#define INCLUDE_COLLFAMILY_H_

#include <algorithm>
#include "./include/default_types.h"

template<class T>  // type of features
struct ElClass {
 public:
    Vec<int> cols;
    Vec<T> vals;

    explicit ElClass(int rank);

    ElClass(int rank, int init_cols[], T init_vals[]);

    ~ElClass() {}

    bool apply_to_object(const Vec<T>& object);
    bool operator==(const ElClass<T>& comp_obj);
};


template<class T>
struct ElColl {
    Vec< ElClass<T> > ecs;

    explicit ElColl(int size = 0) : ecs(size) {}
    ElColl(int size, ElClass<T> init_ecs[]);
    ~ElColl() {}

    bool add(ElClass<T> ec_obj);  // checks for similar elementary classifiers

    Vec<bool> apply_to_object(const Vec<T>& object);
};


template<class T>
struct CollFamily {
    Vec< ElColl<T> > colls;

    CollFamily(int size, ElColl<T> init_colls[]) : colls(size) {
        for (int i = 0; i < size; i++)
            colls[i] = init_colls[i];
    }
    ~CollFamily() {}

    CollFamily<T>& operator+(const CollFamily<T>& W1);  // operator of uniting two families
};

#endif  // INCLUDE_COLLFAMILY_H_
