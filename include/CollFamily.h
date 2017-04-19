/* Copyright 2017 Baytekov Nikita */
#ifndef INCLUDE_COLLFAMILY_H_
#define INCLUDE_COLLFAMILY_H_

#include <algorithm>
#include "./include/default_types.h"

template<typename T>  // type of features
class ElClass {
    Vec<int> cols;
    Vec<T> vals;

 public:
    explicit ElClass(int rank);

    ElClass(int rank, int init_cols[], T init_vals[]);

    ~ElClass() {}

    bool apply_to_object(const Vec<T>& object);
    int get_size() const { return cols.get_size(); }
    const Vec<int> get_cols() const { return cols; }
    const Vec<T> get_vals() const { return vals; }

    bool operator==(const ElClass<T>& comp_obj);
};


template<typename T>
class ElColl {
    Vec< ElClass<T> > ecs;

 public:
    explicit ElColl(int size = 0) : ecs(size) {}
    ElColl(int size, ElClass<T> init_ecs[]);
    ~ElColl() {}

    bool add(ElClass<T> ec_obj);  // checks for similar elementary classifiers

    Vec<bool> apply_to_object(const Vec<T>& object);
    bool vote_func(const Vec<T>& fst_obj, const Vec<T>& sec_obj);

    int get_size() const { return ecs.get_size(); }

    bool operator==(const ElColl<T>& comp_obj);
    ElClass<T>& operator[](int idx) { return ecs[idx]; }
    const ElClass<T>& operator[](int idx) const{ return ecs[idx]; }
};


template<typename T>
class CollFamily {
    Vec< ElColl<T> > colls;

 public:
    CollFamily(int size, ElColl<T> init_colls[]) : colls(size) {
        for (int i = 0; i < size; i++)
            colls[i] = init_colls[i];
    }
    ~CollFamily() {}
    Vec<bool> add(const Vec<ElColl<T> >& new_colls);  // ...

    int get_size() const { return colls.get_size(); }
};

#endif  // INCLUDE_COLLFAMILY_H_
