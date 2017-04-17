/* Copyright 2017 Baytekov Nikita */
#include "../include/default_types.h"

#ifndef INCLUDE_COVERINGHANDLER_H_
#define INCLUDE_COVERINGHANDLER_H_

template<typename T>
class CoveringHandler {
 public:
    bool is_covering(const Mat<bool>& in_mat, const Vec<T>& in_vec);
    Vec<T> build_covering(const Mat<bool>& in_mat);
    Vec<T> complete_to_covering(const Mat<bool>& in_mat, const Vec<T>& in_vec);
    Vec<T> make_covering_deadend(const Mat<bool>& in_mat, const Vec<T>& in_vec);
    Vec<T> standardize_int_form(const Mat<bool>& in_mat, const Vec<T>& in_vec);
    Vec<T> recover_admissibility(const Mat<bool>& in_mat, const Vec<T>& in_vec);
    Vec<int> get_ones_places(const Vec<bool>& in_vec);
};

#endif  // INCLUDE_COVERINGHANDLER_H_
