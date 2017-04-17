/* Copyright 2017 Baytekov Nikita */
#include <type_traits>
#include "./include/SampleHandler.h"
#include "./include/GeneticDualizer.h"

template<typename S, typename T>
GenDualInitiator<S, T>::GenDualInitiator(GeneticDualizer<S, T>* init_parent_ptr, int init_num):
        GeneticInitiator<S, T>(init_num), my_parent(init_parent_ptr), covering_handler() {
}

template<typename S, typename T>
Population<T> GenDualInitiator<S, T>::get_init_population(const SampleSet<S>& sample_set) {
    Population<T> out_popul;
    Mat<bool> gen_matrix = my_parent->get_gen_matrix();
    for (int i = 0; i < this->popul_num; i++) {
        Vec<T> new_genes = covering_handler.build_covering(gen_matrix);
        new_genes = covering_handler.make_covering_deadend(gen_matrix, new_genes);
        if (!out_popul.add_chromo(Chromosome<T>(new_genes)))
            i--;
    }
    return out_popul;
}


template<typename S, typename T>
GenDualMutator<S, T>::GenDualMutator(GeneticDualizer<S, T>* init_parent_ptr,
                                  float init_mrate,
                                  float init_pfrac):
        GeneticMutator<T>(init_pfrac),
        my_parent(init_parent_ptr),
        mutation_rate(init_mrate),
        curr_iter(0) {
}

template<typename S, typename T>
Population<T> GenDualMutator<S, T>::mutate_population(const Population<T>& in_popul) {
    int mut_frac = this->popul_frac * (1.0 - 1.0/(mutation_rate * curr_iter + 1.0));
    int popul_size = in_popul.get_size();
    Population<T> out_popul;
    if (std::is_same<T, bool>::value) {
        for (int i = 0; i < popul_size; i++) {
            int chromo_len = in_popul[i].get_size();
            Vec<T> mut_chromo(chromo_len);
            for (int j = 0; j < chromo_len; j++) {
                if (static_cast<float>(rand_r(time(0))) / (RAND_MAX) > mut_frac)
                    mut_chromo[i] = in_popul[i][j];
                else
                    mut_chromo[i] = !in_popul[i][j];
            }
            out_popul.add_chromo(mut_chromo);
        }
    } else if (std::is_same<T, int>::value) {
        int sample_set_size = this->sample_set.get_size();

        for (int i = 0; i < popul_size; i++) {
            int chromo_len = in_popul[i].get_size();
            Vec<T> mut_chromo(chromo_len);
            for (int j = 0; j < chromo_len; j++) {
                if (static_cast<float>(rand_r(time(0))) / (RAND_MAX) > mut_frac) {
                    mut_chromo[i] = in_popul[i][j];
                } else {
                    int rand_val = rand_r(time(0)) % sample_set_size;
                    while (rand_val == in_popul[i][j])
                        rand_val = rand_r(time(0)) % sample_set_size;
                    mut_chromo[i] = rand_val;
                }
            }
            out_popul.add_chromo(mut_chromo);
        }
    }
    return recover_admissibility(out_popul);
}

template <typename S, typename T>
Population<T> GenDualMutator<S, T>::recover_admissibility(const Population<T>& in_popul) {
    Vec< Vec<T> > in_popul_data = in_popul.get_popul_data();
    Mat<bool> gen_matrix = my_parent->get_gen_matrix();
    Population<T> out_popul;
    Vec<T> new_vec;

    int data_size = in_popul_data.get_size();
    for (int i = 0; i < data_size; i++) {
        new_vec = covering_handler.recover_admissibility(gen_matrix, in_popul_data[i]);
        out_popul.add_chromo(Chromosome<T>(new_vec));
    }
}


template<typename S, typename T>
GeneticDualizer<S, T>::GeneticDualizer(int popul_size, float sel_frac, float max_mut,
                                       float mut_frac, TerminationCriterionType init_tcrit,
                                       float init_tcrit_val):
            GeneticAlgorithm<S, T>(GenDualInitiator<S, T>(this, popul_size),
                                   RouletteSelector<T>(sel_frac),
                                   PanmixiaBreeder<T>(popul_size, kUniform),
                                   GenDualMutator<S, T>(this, popul_size, mut_frac),
                                   init_tcrit,
                                   init_tcrit_val),
            sample_handler(),
            local_basis(),
            basic(),
            valid(),
            target_tag(-1) {
}

template<typename S, typename T>
Vec<ElColl<S> > GeneticDualizer<S, T>::decode_collections(Vec< Vec<T> > in_code) {
    int code_len = in_code.get_size();
    int vec_len = in_code[0].get_size();
    Vec<ElColl<S> > answer_colls(code_len);
    if (std::is_same<T, bool>::value) {
        for (int i = 0; i < code_len; i++)
            for (int j = 0; j < vec_len; j++)
                if (in_code[i][j])
                    answer_colls[i].add(local_basis[j]);
    } else if (std::is_same<T, int>::value) {
        for (int i = 0; i < code_len; i++)
            for (int j = 0; j < vec_len; j++)
                answer_colls[i].add(local_basis[in_code[i][j]]);
    }

    return answer_colls;
}

template<typename S, typename T>
void GeneticDualizer<S, T>::update_costs(Population<T>* in_popul) {
    int chromo_num = in_popul->get_size();
    float qualities[chromo_num] = {};
    float min_val = -1.0;

    for (int i = 0; i < chromo_num; i++) {
        qualities[i] = get_quality(in_popul[i]);
        if (qualities[i] < min_val || min_val < 0.0)
            min_val = qualities[i];
    }

    for (int i = 0; i < chromo_num; i++)
        in_popul[i].set_cost(qualities[i]-min_val+1);
    return 1;
}

template<typename S, typename T>
float GeneticDualizer<S, T>::get_quality(const Chromosome<T>& chromo) {
    ClassSamples<T>& basic_class = basic.get_class(target_tag);
    ClassSamples<T>& valid_class = valid.get_class(target_tag);
    Vec<T> chromo_genes = chromo.get_genes();

    int basic_class_num = basic_class.get_size();
    int valid_class_num = valid_class.get_size();

    float quality_sum = 0.0;
    for (int i = 0; i < basic_class_num; i++)
        for (int j = 0; j < valid_class_num; j++) {
            int genes_len = basic_class.get_size();
            for (int k = 0; k < genes_len; k++)
                if (basic_class[i][k] < valid_class[j][k]) {
                    quality_sum--;
                    break;
                }
            quality_sum++;
        }

    return quality_sum / static_cast<float>(valid_class_num);
}

template<typename S, typename T>
void GeneticDualizer<S, T>::set_init_data(const SampleSet<S>& init_set,
                                          const ElColl<S>& init_lb,
                                          int init_target_tag) {
    sample_handler.make_train_and_valid(init_set, &basic, &valid);
    local_basis = init_lb;
    target_tag = init_target_tag;

    ClassSamples<S> tag_class = init_set.get_class(target_tag);
    SampleSet<S> tag_anticlass = init_set.get_anticlass(target_tag);

    int tclass_size = tag_class.get_size();
    int tanticlass_size = tag_anticlass.get_total_size();

    int lb_size = local_basis.get_size();

    gen_matrix = Mat<bool>(tclass_size * tanticlass_size, local_basis.get_size());
    for (int i = 0; i < tclass_size; i++) {
        for (int j = 0; j < tanticlass_size; j++) {
            Vec<bool> class_vec = local_basis.apply_to_object(tag_class[i]);
            Vec<bool> aclass_vec = local_basis.apply_to_object(tag_anticlass[j]);

            for (int k = 0; k < lb_size; k++)
                gen_matrix[i*tanticlass_size + j] = class_vec[k] && !aclass_vec[k];
        }
    }
    get_gen_matrix();
}
