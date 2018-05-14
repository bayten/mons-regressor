/* Copyright 2017 Baytekov Nikita */

//==========================================
// Genetic Dualizer test 1
//------------------------------------------
// Testing GD for working with WeightedScore
//==========================================

#include <iostream>
#include <ctime>
#include <vector>
#include "log_wrapper.h"
#include "default_types.h"
#include "LBBuilder.h"
#include "GeneticDualizer.h"
#include "CollFamily.h"

#include "SampleSet.h"
#include "SampleHandler.h"

void get_rand_matrix(Mat<bool>* mat, int ones_per_row);

Vec< Vec<bool> > get_all_covers(const Mat<bool>& mat);

bool check_results(const Vec< Vec<bool> >& true_covs, const Vec< Vec<bool> >& pred_covs);

bool merge_colls(Vec< Vec<bool> >* covs, const Vec< Vec<bool> >& iter_covs);

int main(void) {
    srand(time(0));
    init_logging(kDebug);

    GeneticDualizer<int, bool> test_object(120, 0.27, false, kWeightedScore, kDelayConverged, 0.001);
    LOG_(debug) << "Test of GeneticDualizer has been initiated.";

    for (int i = 0; i < 10; i++) {
        Mat<bool> rand_mat(30, 20);
        get_rand_matrix(&rand_mat, 7);
        LOG_(debug) << "Random L matrix:";
        LOG_(debug) << rand_mat;

        Vec< Vec<bool> > true_covs = get_all_covers(rand_mat);
        // LOG_(debug) << "True covs:" << true_covs;

        test_object.set_matrix(rand_mat);

        clock_t begin_tm = clock();
        int term_crit = 0;
        // std::vector<clock_t> iter_cuts(0);
        Vec< Vec<bool> > covs(0);

        while (term_crit < 50) {
            Vec< Vec<bool> > iter_covs = test_object.execute_ga();
            // iter_cuts.push_back(clock());
            if(merge_colls(&covs, iter_covs))
                term_crit = 0;
            else
                term_crit++;
        }
        clock_t total_end = clock();

        // LOG_(debug) << "Estimated coverages: " << covs;
        check_results(true_covs, covs);

        // for (unsigned int i = 0; i < iter_cuts.size(); i++) {
        //     double delta_tm = double(total_end - begin_tm) / CLOCKS_PER_SEC;
        //     LOG_(debug) << "Iter #" << i << " elapsed time: " << delta_tm;
        // }

        double total_delta_tm = double(total_end - begin_tm) / CLOCKS_PER_SEC;
        LOG_(debug) << "Total elapsed time: " << total_delta_tm << " s";
    }
    return 0;
}

void get_rand_matrix(Mat<bool>* mat, int ones_per_row) {
    int sx = mat->get_sx();
    int sy = mat->get_sy();

    for (int i = 0; i < sx; i++) {
        for (int j = 0; j < sy; j++)
            (*mat)[i][j] = 0;

        for (int j = 0; j < ones_per_row; j++) {
            int rand_place = rand() % (sy-j);
            int k = 0, ind = 0;
            while (k != rand_place) {
                if ((*mat)[i][ind] == 0)
                    k++;
                ind++;
            }
            (*mat)[i][ind-1] = 1;
        }
    }
    return;
}

Vec< Vec<bool> > get_all_covers(const Mat<bool>& mat) {
    int sy = mat.get_sy();

    Vec< Vec<bool> > covs(0);
    Vec<bool> cov(sy);
    CoveringHandler<bool> cov_handler;
    unsigned long long limit = std::pow(2, sy);
    LOG_(debug) << "Total test subjects:" << limit;
    for (unsigned long code_cov = 1; code_cov < limit; code_cov++) {
        unsigned long code_dummy = code_cov;

        for (int i = sy-1; i >= 0; i--) {
            cov[i] = code_dummy & 1;
            code_dummy >>= 1;
        }
        // LOG_(debug) << "Code:" << code_cov << " vec:" << cov;

        if (cov_handler.is_covering(mat, cov)) {
            Vec<bool> dead_cov = cov_handler.make_covering_deadend(mat, cov);
            int ready_covs = covs.get_size();
            bool trig = true;

            for (int i = 0; i < ready_covs; i++) {
                if (dead_cov == covs[i]) {
                    trig = false;
                    break;
                }
            }

            if (trig) {
                covs.append(dead_cov);
            }
        }
    }
    return covs;
}

bool check_results(const Vec< Vec<bool> >& true_covs, const Vec< Vec<bool> >& pred_covs) {
    int true_num = true_covs.get_size();
    int pred_num = pred_covs.get_size();
    Vec<bool> matched_covs(true_num);
    Vec<bool> incorrect_covs(pred_num);
    int matched = 0;
    int incorrect = 0;
    for (int i = 0; i < true_num; i++)
        matched_covs[i] = false;

    for (int i = 0; i < pred_num; i++) {
        bool found_flag = false;
        incorrect_covs[i] = false;

        for (int j = 0; j < true_num; j++) {
            if (pred_covs[i] == true_covs[j]) {
                matched_covs[j] = true;
                matched++;
                found_flag = true;
                break;
            }
        }
        if (!found_flag) {
            incorrect++;
            incorrect_covs[i] = true;
        }
    }

    if (matched == true_num && true_num == pred_num) {
        LOG_(debug) << "Test was passed successfully!";
        return 1;
    }

    LOG_(debug) << "Stats:";
    LOG_(debug) << "Accuracy:" << float(matched) / float(true_num);
    LOG_(debug) << "Incorrect:" << incorrect;
    LOG_(debug) << "Not found:" << true_num - matched;
    return 0;
}

bool merge_colls(Vec< Vec<bool> >* covs, const Vec< Vec<bool> >& iter_covs) {
    int icovs_sz = iter_covs.get_size();
    int covs_sz = covs->get_size();
    bool global_flag = false;

    for (int i = 0; i < icovs_sz; i++) {
        bool trigger = true;
        for (int j = 0; j < covs_sz; j++) {
            if (iter_covs[i] == (*covs)[j]) {
                trigger = false;
                break;
            }
        }
        if (trigger) {
            global_flag = true;
            covs->append(iter_covs[i]);
        }
    }

    return global_flag;
}
