/* Copyright 2017 Baytekov Nikita */
#include "./include/CollFamily.h"

template<typename T>
ElClass<T>::ElClass(int rank):
        cols(rank), vals(rank) {
    for (int i = 0; i < rank; i++)
        cols[i] = -1;
}

template<typename T>
ElClass<T>::ElClass(int rank, int init_cols[], T init_vals[]):
        cols(rank), vals(rank) {
    for (int i = 0; i < rank; i++) {
        cols[i] = init_cols[i];
        vals[i] = init_vals[i];
    }
}

template<typename T>
bool ElClass<T>::apply_to_object(const Vec<T>& object) {
    int ec_rank = cols.get_size();
    for (int i = 0; i < ec_rank; i++)
        if (object[cols[i]] == vals[i])
            return 0;
    return 1;
}

template<typename T>
bool ElClass<T>::operator==(const ElClass<T>& comp_obj) {
    int my_size = cols.get_size();
    if (my_size != comp_obj.get_size())
        return 0;
    cols = cols.sort();
    vals = vals.sort();
    Vec<int> cols2 = comp_obj.get_cols().sort();
    Vec<T> vals2 = comp_obj.get_vals().sort();

    for (int i = 0; i < my_size; i++) {
        if (cols[i] != cols2[i] ||
            vals[i] != vals2[i])
            return 0;
    }
}


template<typename T>
ElColl<T>::ElColl(int size, ElClass<T> init_ecs[]) : ecs(size) {
    for (int i = 0; i < size; i++)
        ecs[i] = init_ecs[i];
}

template<typename T>
bool ElColl<T>::add(ElClass<T> ec_obj) {
    int ec_num = ecs.get_size();
    for (int i = 0; i < ec_num; i++)  // checking if there any similar elementary classifiers
        if (ecs[i] == ec_obj)
            return 0;
    ecs.append(ec_obj);
}

template<typename T>
Vec<bool> ElColl<T>::apply_to_object(const Vec<T>& object) {
    int ec_num = ecs.get_size();
    Vec<bool> bin_coll(ec_num);

    for (int i = 0; i < ec_num; i++)
        bin_coll[i] = ecs[i].apply_to_object(object);

    return bin_coll;
}

template<typename T>
bool ElColl<T>::operator==(const ElColl<T>& comp_obj) {
    int ec_num = ecs.get_size();
    if (ec_num != comp_obj.get_size())
        return 0;

    bool was_found = 0;
    for (int i = 0; i < ec_num; i++) {
        for (int j = 0; j < ec_num; j++) {
            if (comp_obj[i] == ecs[j]) {
                was_found = 1;
                break;
            }
        }
        if (!was_found)
            return 0;
        was_found = 1;
    }
    return 1;
}


template<typename T>
Vec<bool> CollFamily<T>::add(const Vec<ElColl<T> >& new_colls) {
    int new_coll_num = new_colls.get_size();
    Vec<bool> add_state(new_coll_num);

    for (int i = 0; i < new_coll_num; i++) {
        add_state[i] = 1;
        int coll_num = colls.get_size();  // size of coll vector differs from iter to iter

        for (int j = 0; j < coll_num; j++) {
            if (colls[j] == new_colls[i]) {
                add_state[i] = 0;
                break;
            }
        }

        if (add_state[i])
            colls.append(new_colls[i]);
    }
    return add_state;
}
