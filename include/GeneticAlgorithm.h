/* Copyright 2017 Baytekov Nikita */
#ifndef INCLUDE_GENETICALGORITHM_H_
#define INCLUDE_GENETICALGORITHM_H_

#include "../include/default_types.h"
#include "../include/SampleHandler.h"

// TODO(Baytekov): Add weights!

template<typename T>
class Chromosome {
    Vec<T> genes;
    float score;

 public:
    explicit Chromosome(Vec<T> init_genes, float init_score = 0.0);
    explicit Chromosome(T init_gene, float init_score = 0.0);
    explicit Chromosome(const Chromosome<T>& chromo_obj);
    ~Chromosome();

    int get_size() const { return genes.get_size(); }
    const Vec<T>& get_genes() const { return genes; }
    float get_score() const { return score; }
    void set_score(float score_val) { score = score_val; }

    Chromosome<T>& operator=(const Chromosome<T>& chromo_obj);
    bool operator==(const Chromosome<T>& chromo_obj);
    T& operator[] (int idx) { return genes[idx]; }
};

template<typename T>
class Population {
    Vec< Chromosome<T> > chromo_vec;

 public:
    explicit Population(Vec< Chromosome<T> > init_vec) : chromo_vec(init_vec) {}
    Population() {}
    ~Population() {}

    bool add_chromo(Chromosome<T> add_obj);
    bool del_chromo(int idx);

    const Vec< Chromosome<T> >& get_chromos() const { return chromo_vec; }
    Vec< Vec<T> > get_popul_data() const;
    int get_size() const { return chromo_vec.get_size(); }

    float get_avg_score() const;
    const Chromosome<T>& get_best_ind() const;

    Chromosome<T>& operator[] (int idx) { return chromo_vec[idx]; }
};


template<typename S, typename T>
class GeneticInitiator {
    int popul_num;

 public:
    explicit GeneticInitiator(int init_num = 0) : popul_num(init_num) {}
    virtual ~GeneticInitiator() = 0;
    virtual Population<T> get_init_population(const SampleSet<S>& sample_set) = 0;
};


template<typename T>
class GeneticSelector {
    float popul_frac;
 public:
    explicit GeneticSelector(int init_pfrac = 0.5) : popul_frac(init_pfrac) {}
    virtual ~GeneticSelector() = 0;
    virtual Population<T> select_population(const Population<T>& in_popul) = 0;
};

template<typename T>
class TournamentSelector : public GeneticSelector<T> {
 public:
    explicit TournamentSelector(int init_pfrac = 0.5) : GeneticSelector<T>(init_pfrac) {}
    virtual ~TournamentSelector();
    virtual Population<T> select_population(const Population<T>& in_popul);
};

template<typename T>
class RouletteSelector : public GeneticSelector<T> {
 public:
    explicit RouletteSelector(int init_pfrac = 0.5) : GeneticSelector<T>(init_pfrac) {}
    virtual ~RouletteSelector();
    virtual Population<T> select_population(const Population<T>& in_popul);
};

template<typename T>
class RankingSelector : public GeneticSelector<T> {
    int uniform_thresh;
 public:
    explicit RankingSelector(int init_pfrac = 0.5, int init_uniform = 0);
    virtual ~RankingSelector();
    virtual Population<T> select_population(const Population<T>& in_popul);
};

template<typename T>
class SigmaTruncSelector : public GeneticSelector<T> {
 public:
    explicit SigmaTruncSelector(int init_pfrac = 0.5) : GeneticSelector<T>(init_pfrac) {}
    virtual ~SigmaTruncSelector();
    virtual Population<T> select_population(const Population<T>& in_popul);
};

enum CrossoverType {
    kOnePoint = 0,
    kMultiPoint = 1,
    kUniform = 2,
    kScoredUniform = 3
};
template<typename T>
class GeneticBreeder {
    int popul_num;
    CrossoverType cross_type;
    float cross_param;

 public:
    explicit GeneticBreeder(int init_pnum = 0, CrossoverType init_ctype = kOnePoint,
                            float init_cparam = 1.0);
    virtual ~GeneticBreeder() = 0;
    virtual Population<T> breed_new_population(const Population<T>& in_popul) = 0;

 protected:
    virtual Vec< Chromosome<T> >apply_crossover(const Chromosome<T>& fst,
                                                const Chromosome<T>& sec);
};

template<typename T>
class PanmixiaBreeder : public GeneticBreeder<T> {
 public:
    explicit PanmixiaBreeder(int init_pnum = 0, CrossoverType init_ctype = kOnePoint,
                             float init_cparam = 1.0);
    ~PanmixiaBreeder();
    virtual Population<T> breed_new_population(const Population<T>& in_popul);
};

enum InOutBreederType {
    kInPhenoType  = 0,
    kInGenoType   = 1,
    kOutPhenoType = 2,
    kOutGenoType  = 3
};
template<typename T>
class InOutBreeder : public GeneticBreeder<T> {
    InOutBreederType breeder_type;

 public:
    explicit InOutBreeder(int init_pnum = 0, InOutBreederType init_btype = kInPhenoType,
                          CrossoverType init_ctype = kOnePoint, float init_cparam = 1.0);
    ~InOutBreeder();
    virtual Population<T> breed_new_population(const Population<T>& in_popul);

 private:
    virtual const Chromosome<T>& find_match(const Population<T>& in_popul,
                                            const Chromosome<T>& suitor);
    float find_geno_diff(const Chromosome<T>& fst, const Chromosome<T>& sec);
};


template<typename T>
class GeneticMutator {
    float popul_frac;
 public:
    explicit GeneticMutator(int init_pfrac = 0.0) : popul_frac(init_pfrac) {}
    virtual ~GeneticMutator() = 0;
    virtual Population<T> mutate_population(const Population<T>& in_popul) = 0;
};


enum TerminationCriterionType {
    kPopulConverged = 0,  // term_crit_val = min delta between two populations
    kBestConverged  = 1,  // term_crit_val = min delta between two fitnesses of best individual
    kMaxPopulNum    = 2   // term_crit_val = amount of populations
};
template<typename S, typename T>
class GeneticAlgorithm {
 protected:
    SampleSet<S> sample_set;

    GeneticInitiator<S, T> initiator;
    GeneticSelector<T> selector;
    GeneticBreeder<T> breeder;
    GeneticMutator<T> mutator;

    TerminationCriterionType term_crit;
    float term_crit_val;

 public:
    GeneticAlgorithm(GeneticInitiator<S, T> init_initiator, GeneticSelector<T> init_selector,
                     GeneticBreeder<T> init_breeder, GeneticMutator<T> init_mutator,
                     TerminationCriterionType init_tcrit = kPopulConverged,
                     float init_tcrit_val = 1.0);
    virtual ~GeneticAlgorithm();

    virtual Vec< Vec<T> > execute_ga();
    virtual bool check_term_crit(Population<T> curr_popul, Population<T> new_popul,
                                 int64_t iter_cnt);

     virtual void set_sample_set(SampleSet<T> init_set) { sample_set = init_set; }
     virtual void update_costs(Population<T>* in_popul) = 0;
};
#endif  // INCLUDE_GENETICALGORITHM_H_
