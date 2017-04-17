/* Copyright 2017 Baytekov Nikita */
#include "./include/GeneticAlgorithm.h"

template<typename T>
Chromosome<T>::Chromosome(Vec<T> init_genes, float init_score):
        genes(init_genes), score(init_score) {
}

template<typename T>
Chromosome<T>::Chromosome(T init_gene, float init_score):
        genes(1, &init_gene), score(init_score) {
}

template<typename T>
Chromosome<T>::Chromosome(const Chromosome<T>& chromo_obj):
        genes(chromo_obj.get_genes()), score(chromo_obj.get_score()) {
}

template<typename T>
Chromosome<T>& Chromosome<T>::operator=(const Chromosome<T>& chromo_obj) {
    genes = chromo_obj.get_genes();
    score = chromo_obj.get_score();
    return (*this);
}

template<typename T>
bool Chromosome<T>::operator==(const Chromosome<T>& chromo_obj) {
    int chromo_size = genes.get_size();
    if (chromo_size != chromo_obj.get_size() || score != chromo_obj.get_score())
        return 0;

    Vec<T> obj_genes = chromo_obj.get_genes();
    for (int i = 0; i < chromo_size; i++)
        if (obj_genes[i] != genes[i])
            return 0;
    return 1;
}


template<typename T>
bool Population<T>::add_chromo(Chromosome<T> add_obj) {
    int chromo_num = chromo_vec.get_size();
    T add_obj_data = add_obj.get_data();

    for (int i = 0; i < chromo_num; i++)
        if (add_obj_data == chromo_vec[i].get_genes())
            return 0;
    chromo_vec.append(add_obj);

    return 1;
}

template<typename T>
bool Population<T>::del_chromo(int idx) {
    return chromo_vec.erase(idx);
}

template<typename T>
float Population<T>::get_avg_score() const {
    int chromo_num = chromo_vec.get_size();
    float total_score = 0.0;
    for (int i = 0; i < chromo_num; i++)
        total_score += chromo_vec[i].get_score();

    return total_score;
}

template<typename T>
const Chromosome<T>& Population<T>::get_best_ind() const {
    int chromo_num = chromo_vec.get_size();
    int idx = 0;
    float max_score = chromo_vec[idx];
    for (int i = 1; i < chromo_num; i++)
        if (chromo_vec[i].get_score() > max_score) {
            idx = i;
            max_score = chromo_vec[idx];
        }

    return chromo_vec[idx];
}

template<typename T>
Vec< Vec<T> > Population<T>::get_popul_data() const {
    int chromo_len = chromo_vec.get_size();
    Vec< Vec<T> > out_vec(chromo_len);

    for (int i = 0; i < chromo_len; i++)
        out_vec[i] = chromo_vec[i].get_genes();
    return out_vec;
}


template<typename T>
Population<T> TournamentSelector<T>::select_population(const Population<T>& in_popul) {
    Population<T> out_popul;
    int popul_size = in_popul.get_size();
    int needed_size = in_popul * this->popul_frac;

    for (int i = 0; i < needed_size; i++) {
        Chromosome<T> fst_chromo = in_popul[rand_r(time(0)) % popul_size];
        Chromosome<T> sec_chromo = in_popul[rand_r(time(0)) % popul_size];

        if (fst_chromo == sec_chromo) {
            i--;
            continue;
        } else if (fst_chromo.get_score() >= sec_chromo.get_score()) {
            if (!out_popul.add_chromo(fst_chromo)) {
                i--;
                continue;
            }
        } else if (!out_popul.add_chromo(sec_chromo)) {
            i--;
            continue;
        }
    }
    return out_popul;
}


template<typename T>
Population<T> RouletteSelector<T>::select_population(const Population<T>& in_popul) {
    Population<T> out_popul;
    int popul_size = in_popul.get_size();
    int needed_size = in_popul * this->population_fraction;
    Vec<float> probs(popul_size);  // probs will be stored here
    float total_score = 0;
    for (int i = 0; i < popul_size; i++) {
        probs[i] = out_popul[i].get_score();
        total_score += probs[i];
    }
    for (int i = 0; i < popul_size; i++)
        probs[i] /= total_score;

    for (int i = 0; i < needed_size; i++) {
        float rand_frac = static_cast<float>(rand_r(time(0))) / (RAND_MAX);
        int chromo_idx = 0;
        for (; chromo_idx < popul_size; chromo_idx++) {
            if (rand_frac < probs[chromo_idx])
                break;
            rand_frac -= probs[chromo_idx];
        }
        if (!out_popul.add_chromo(in_popul[chromo_idx]))
            i--;
    }
    return out_popul;
}

template<typename T>
RankingSelector<T>::RankingSelector(int init_pfrac, int init_uniform) :
        GeneticSelector<T>(init_pfrac), uniform_thresh(init_uniform) {
}

template<typename T>
Population<T> RankingSelector<T>::select_population(const Population<T>& in_popul) {
    Population<T> out_popul;
    int popul_size = in_popul.get_size();
    int needed_size = in_popul * this->popul_frac;
    Vec<float> probs(popul_size);

    for (int i = 0; i < popul_size; i++)
        probs[i] = in_popul.ind_vec[i].get_score();
    Vec<int> indices = probs.sort_indices();

    if (uniform_thresh > 1) {
        float val = 1.0/uniform_thresh;
        for (int i = 0; i < popul_size; i++)
            probs[i] = (popul_size-indices[i] <= uniform_thresh) ? val : 0;
    } else {
        float divisor = popul_size*(popul_size+1)/2.0;
        for (int i = 0; i < popul_size; i++)
            probs[i] = (popul_size-indices[i])/divisor;
    }

    for (int i = 0; i < needed_size; i++)
        if (!out_popul.add_chromo(in_popul[rand_r(time(0)) % popul_size]))
            i--;
    return out_popul;
}


template<typename T>
Population<T> SigmaTruncSelector<T>::select_population(const Population<T>& in_popul) {
    Population<T> out_popul;
    int popul_size = in_popul.get_size();
    int needed_size = in_popul * this->popul_frac;
    Vec<float> probs(popul_size);

    for (int i = 0; i < popul_size; i++)
        probs[i] = in_popul.ind_vec[i].get_score();
    Vec<int> indices = probs.sort_indices();

    for (int i = 0; i < popul_size; i++)
        if (popul_size-indices[i] <= needed_size)
            out_popul.add_chromo(in_popul[i]);
    return out_popul;
}


template<typename T>
GeneticBreeder<T>::GeneticBreeder(int init_pnum, CrossoverType init_ctype, float init_cparam):
        popul_num(init_pnum), cross_type(init_ctype), cross_param(init_cparam) {
}

template<typename T>
Vec< Chromosome<T> > GeneticBreeder<T>::apply_crossover(const Chromosome<T>& fst,
                                                        const Chromosome<T>& sec) {
    int chromo_len = fst.get_size();
    int rand_val = rand_r(time(0)) % chromo_len;
    Vec<T> fst_child_genes(chromo_len);
    Vec<T> sec_child_genes(chromo_len);
    Vec< Chromosome<T> > output_vec;

    int fst_score = fst.get_score();
    int total_score = fst_score + sec.get_score();

    switch (cross_type) {
        case kOnePoint:
            for (int i = 0; i <= rand_val; i++) {
                fst_child_genes[i] = fst[i];
                sec_child_genes[i] = sec[i];
            }
            for (int i = rand_val+1; i < chromo_len; i++) {
                fst_child_genes[i] = sec[i];
                sec_child_genes[i] = fst[i];
            }
            output_vec.append(Chromosome<T>(fst_child_genes));
            output_vec.append(Chromosome<T>(sec_child_genes));
            break;

        case kMultiPoint:
            for (int i = 0; i < cross_param; i++) {
                for (int j = 0; j <= rand_val; j++) {
                    fst_child_genes[j] = fst[j];
                    sec_child_genes[j] = sec[j];
                }
                for (int j = rand_val+1; j < chromo_len; j++) {
                    fst_child_genes[j] = sec[j];
                    sec_child_genes[j] = fst[j];
                }
                output_vec.append(Chromosome<T>(fst_child_genes));
                output_vec.append(Chromosome<T>(sec_child_genes));

                rand_val = rand_r(time(0)) % chromo_len;
            }
            break;

        case kUniform:
            for (int i = 0; i < chromo_len; i++)
                fst_child_genes[i] = (rand_r(time(0)) % 2) ? fst[i] : sec[i];
            output_vec.append(Chromosome<T>(fst_child_genes));
            break;

        case kScoredUniform:
            for (int i = 0; i < chromo_len; i++) {
                if (rand_r(time(0)) % total_score > fst_score)
                    fst_child_genes[i] = fst[i];
                else
                    fst_child_genes[i] = sec[i];
            }
            output_vec.append(Chromosome<T>(fst_child_genes));
            break;
    }

    return output_vec;
}


template<typename T>
PanmixiaBreeder<T>::PanmixiaBreeder(int init_pnum, CrossoverType init_ctype, float init_cparam):
        GeneticBreeder<T>(init_pnum, init_ctype, init_cparam) {
}

template<typename T>
Population<T> PanmixiaBreeder<T>::breed_new_population(const Population<T>& in_popul) {
    Population<T> out_popul;
    int in_popul_size = in_popul.get_size();
    if (in_popul_size <= 1)
        return in_popul;

    for (int i = 0; i < this->popul_num;) {
        Chromosome<T> fst_parent = in_popul[rand_r(time(0)) % in_popul_size];
        Chromosome<T> sec_parent = in_popul[rand_r(time(0)) % in_popul_size];
        while (sec_parent == fst_parent)
            sec_parent = in_popul[rand_r(time(0)) % in_popul_size];
        Vec< Chromosome<T> > children = apply_crossover(fst_parent, sec_parent);
        int children_num = children.get_size();

        for (int j = 0; j < children_num; j++)
            if (out_popul.add_chromo(children[i])) {
                i++;
                if (i >= this->popul_num)
                    break;
            }
    }
    return out_popul;
}

template<typename T>
InOutBreeder<T>::InOutBreeder(int init_pnum, InOutBreederType init_btype,
                              CrossoverType init_ctype, float init_cparam):
        GeneticBreeder<T>(init_pnum, init_ctype, init_cparam), breeder_type(init_btype) {
}

template<typename T>
Population<T> InOutBreeder<T>::breed_new_population(const Population<T>& in_popul) {
    Population<T> out_popul;
    int in_popul_size = in_popul.get_size();
    if (in_popul_size <= 1)
        return in_popul;

    for (int i = 0; i < this->popul_num; i++) {
        Chromosome<T> fst_parent = in_popul[rand_r(time(0)) % in_popul_size];
        Chromosome<T> sec_parent = find_match(in_popul, fst_parent);
        Vec< Chromosome<T> > children = apply_crossover(fst_parent, sec_parent);
        int children_num = children.get_size();

        for (int j = 0; j < children_num; j++)
            if (out_popul.add_chromo(children[i])) {
                i++;
                if (i >= this->popul_num)
                    break;
            }
    }
    return out_popul;
}

template<typename T>
const Chromosome<T>& InOutBreeder<T>::find_match(const Population<T>& in_popul,
                                                 const Chromosome<T>& suitor) {
    int in_popul_size = in_popul.get_size();
    int partner_idx = (in_popul[0] == suitor) ? 1 : 0;
    float param_diff = 0.0, curr_param_diff = 0.0;
    float suitor_score = suitor.get_score();

    switch (breeder_type) {
        case kInPhenoType:
            param_diff = fabs(in_popul[partner_idx].get_score() - suitor_score);
            for (int i = (partner_idx+1); i < in_popul_size; i++) {
                curr_param_diff = fabs(in_popul[i].get_score() - suitor_score);
                if (curr_param_diff < param_diff && !(in_popul[i] == suitor)) {
                    partner_idx = i;
                    param_diff = curr_param_diff;
                }
            }
            break;

        case kInGenoType:
            param_diff = find_geno_diff(in_popul[partner_idx], suitor);
            for (int i = (partner_idx+1); i < in_popul_size; i++) {
                curr_param_diff = in_popul[i] % suitor;
                if (curr_param_diff < param_diff && !(in_popul[i] == suitor)) {
                    partner_idx = i;
                    param_diff = curr_param_diff;
                }
            }
            break;

        case kOutPhenoType:
            param_diff = fabs(in_popul[partner_idx].get_score() - suitor_score);
            for (int i = (partner_idx+1); i < in_popul_size; i++) {
                curr_param_diff = fabs(in_popul[i].get_score() - suitor_score);
                if (curr_param_diff > param_diff && !(in_popul[i] == suitor)) {
                    partner_idx = i;
                    param_diff = curr_param_diff;
                }
            }
            break;

        case kOutGenoType:
            param_diff = find_geno_diff(in_popul[partner_idx], suitor);
            for (int i = (partner_idx+1); i < in_popul_size; i++) {
                curr_param_diff = in_popul[i] % suitor;
                if (curr_param_diff > param_diff && !(in_popul[i] == suitor)) {
                    partner_idx = i;
                    param_diff = curr_param_diff;
                }
            }
            break;
    }
    return in_popul[partner_idx];
}

template<typename T>
float InOutBreeder<T>::find_geno_diff(const Chromosome<T>& fst, const Chromosome<T>& sec) {
    int chromo_len = fst.get_size();
    float sum = 0.0;
    for (int i = 0; i < chromo_len; i++)
        sum += (fst[i]-sec[i])*(fst[i]-sec[i]);
    return sqrt(sum);
}

template<typename S, typename T>
GeneticAlgorithm<S, T>::GeneticAlgorithm(GeneticInitiator<S, T> init_initiator, \
                                         GeneticSelector<T>  init_selector, \
                                         GeneticBreeder<T>   init_breeder, \
                                         GeneticMutator<T>   init_mutator, \
                                         TerminationCriterionType init_term_crit,
                                         float init_term_crit_val) :
                initiator(init_initiator),
                selector(init_selector),
                breeder(init_breeder),
                mutator(init_mutator),
                term_crit(init_term_crit),
                term_crit_val(init_term_crit_val) {
}

template<typename S, typename T>
Vec< Vec<T> > GeneticAlgorithm<S, T>::execute_ga() {
    Population<T> curr_population = initiator->get_init_population(sample_set);
    int64_t iter_count = 0L;

    while (true) {
        curr_population.update_costs();

        Population<T> breed_population = selector.select_population(curr_population);
        Population<T> new_population = breeder.breed_new_population(breed_population);

        new_population = mutator.mutate_population(new_population);
        new_population.update_costs();

        if (check_term_crit(curr_population, new_population, iter_count))
            return new_population.get_popul_data();
    }
}

template<typename S, typename T>
bool GeneticAlgorithm<S, T>::check_term_crit(Population<T> curr_popul, \
                                             Population<T> new_popul, \
                                             int64_t iter_cnt) {
    switch (term_crit) {
        case kPopulConverged:
            if (new_popul.get_avg_score() - curr_popul.get_avg_score() < term_crit_val)
                return 1;
            return 0;

        case kBestConverged:
            float new_best_score = new_popul.get_best_ind().get_score();
            float curr_best_score = curr_popul.get_best_ind().get_score();
            if (new_best_score - curr_best_score < term_crit_val)
                return 1;
            return 0;

        case kMaxPopulNum:
            if (iter_cnt < static_cast<int64_t>(term_crit_val))
                return 1;
            return 0;
    }
}
