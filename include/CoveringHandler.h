/* Copyright 2017 Baytekov Nikita */
#include <type_traits>
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



template<typename T>
bool CoveringHandler<T>::is_covering(const Mat<bool>& in_mat, const Vec<T>& in_vec) {
    int mat_sx = in_mat.get_sx();
    int mat_sy = in_mat.get_sy();

    if (std::is_same<T, bool>::value) {
        bool was_found = 0;
        for (int i = 0; i < mat_sx; i++) {
            for (int j = 0; j < mat_sy; j++) {
                if (in_vec[j] && in_mat[i][j]) {
                    was_found = 1;
                    break;
                }
            }
            if (!was_found)
                return 0;

            was_found = 0;
        }
    } else if (std::is_same<T, int>::value) {
        for (int i = 0; i < mat_sx; i++)
            if (!in_mat[i][in_vec[i]])
                return 0;
    }
    return 1;
}

template<typename T>
Vec<T> CoveringHandler<T>::build_covering(const Mat<bool>& in_mat) {
    unsigned int rand_seed = static_cast<unsigned int>(time(0));
    if (std::is_same<T, bool>::value) {
        int mat_sy = in_mat.get_sy();
        Vec<T> new_vec(mat_sy);
        for (int i = 0; i < mat_sy; i++)
            new_vec[i] = rand_r(&rand_seed) % 2;

        if (!is_covering(in_mat, new_vec))
            new_vec = complete_to_covering(in_mat, new_vec);

        return new_vec;
    } else if (std::is_same<T, int>::value) {
        int mat_sx = in_mat.get_sx();
        Vec<T> new_vec(mat_sx);
        for (int i = 0; i < mat_sx; i++) {
            Vec<int> ones_vec = get_ones_places(in_mat[i]);
            new_vec[i] = ones_vec[rand_r(&rand_seed) % ones_vec.get_size()];
        }
        return standardize_int_form(in_mat, new_vec);
    }
}

template<typename T>
Vec<T> CoveringHandler<T>::complete_to_covering(const Mat<bool>& in_mat,
                                                const Vec<T>& in_vec) {
    int mat_sx = in_mat.get_sx();
    int mat_sy = in_mat.get_sy();
    unsigned int rand_seed = static_cast<unsigned int>(time(0));
    Vec<T> cover_vec = in_vec;
    if (std::is_same<T, bool>::value) {
        bool was_found = 0;

        for (int i = 0; i < mat_sx; i++) {
            for (int j = 0; j < mat_sy; j++) {
                if (cover_vec[j] && in_mat[i][j]) {
                    was_found = 1;
                    break;
                }
            }
            if (!was_found) {
                Vec<int> ones_vec = get_ones_places(in_mat[i]);
                cover_vec[ones_vec[rand_r(&rand_seed) % ones_vec.get_size()]] = 1;
            }
        }
    } else if (std::is_same<T, int>::value) {
        for (int i = 0; i < mat_sx; i++) {
            if (!in_mat[i][cover_vec[i]]) {
                Vec<int> ones_vec = get_ones_places(in_mat[i]);
                cover_vec[i] = ones_vec[rand_r(&rand_seed) % ones_vec.get_size()];
            }
        }
    }
    return cover_vec;
}

template<typename T>
Vec<T> CoveringHandler<T>::make_covering_deadend(const Mat<bool>& in_mat,
                                                 const Vec<T>& in_vec) {
    if (std::is_same<T, bool>::value) {
        Vec<T> deadend_vec = in_vec;
        Vec<T> prev_vec = deadend_vec;
        int mat_sy = in_mat.get_sy();
        int curr_idx = 0;

        for (int i = 0; i < mat_sy; i++) {
            deadend_vec.erase(curr_idx);
            if (is_covering(in_mat, deadend_vec)) {
                prev_vec = deadend_vec;
            } else {
                deadend_vec = prev_vec;
                curr_idx++;
            }
        }
        return deadend_vec;
    } else if (std::is_same<T, int>::value) {
        return standardize_int_form(in_mat, in_vec);
    }
}

template<typename T>
Vec<T> CoveringHandler<T>::recover_admissibility(const Mat<bool>& in_mat,
                                                 const Vec<T>& in_vec) {
    Vec<T> new_vec = complete_to_covering(in_mat, in_vec);
    return make_covering_deadend(in_mat, new_vec);
}

template<typename T>
Vec<T> CoveringHandler<T>::standardize_int_form(const Mat<bool>& in_mat,
                                                const Vec<T>& in_vec) {
    if (std::is_same<T, int>::value != true)
        return in_vec;  // TODO(Baytekov): interpret as error!

    Vec<T> sorted_vec = in_vec.sort();
    int vec_len = sorted_vec.get_size();
    Vec<T> unique_vals(1, &(sorted_vec[0]));

    for (int i = 1; i < vec_len; i++)
        if (unique_vals[-1] != sorted_vec[i])
            unique_vals.append(sorted_vec[i]);

    int unique_size = unique_vals.get_size();
    int mat_sx = in_mat.get_sx();
    Vec<int> standard_vec(mat_sx);

    for (int i = 0; i < mat_sx; i++) {
        Vec<int> ones_vec = get_ones_places(in_mat[i]);
        int ones_size = ones_vec.get_size();
        bool was_set = 0;
        for (int j = 0; j < unique_size; j++) {
            for (int k = 0; k < ones_size; k++)
                if (unique_vals[j] == ones_vec[k]) {
                    standard_vec[i] = unique_vals[j];
                    was_set = 1;
                    break;
                }
            if (was_set)
                break;
        }
    }
}

template<typename T>
Vec<int> CoveringHandler<T>::get_ones_places(const Vec<bool>& in_vec) {
    Vec<int> out_vec;
    int vec_size = in_vec.get_size();

    for (int i = 0; i < vec_size; i++)
        if (in_vec[i])
            out_vec.append(i);

    return out_vec;
}

#endif  // INCLUDE_COVERINGHANDLER_H_
