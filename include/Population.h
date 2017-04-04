/* Copyright 2017 Baytekov Nikita */
#include "../include/CollFamily.h"
#include "../include/SampleHandler.h"

#ifndef INCLUDE_POPULATION_H_
#define INCLUDE_POPULATION_H_

template<class T>
class Individual {
 protected:
    T data;
    float cost;

 public:
    explicit Individual(T init_data) : data(init_data), cost(0.0) {}
    ~Individual() {}

    float get_cost() const { return cost; }
    T get_data() const { return data; }
    void set_cost(float cost_val) { cost = cost_val; }

    virtual Individual<T> operator*(const Individual<T>& cross_obj);  // crossover operator
};

template<class T>
class CoveringIndividual : public Individual < ElColl<T> > {
 public:
    explicit CoveringIndividual(ElColl<T> init_coll) : Individual< ElColl<T> >(init_coll) {}
    ~CoveringIndividual() {}

    float get_quality(const SampleSet<T>& basic, \
                      const SampleSet<T>& valid, \
                      int target_tag);

    Individual<ElColl<T> > operator*(const Individual<ElColl<T> >& cross_obj);
};

// TODO(Baytekov): Solve Population class problems with CoveringIndividual - templates!
template<class T>
class Population {
    Vec<Individual<T> > ind_vec;

 public:
    explicit Population(Vec<Individual<T> > init_vec);
    ~Population();

    void update_costs(const SampleSet<T>& basic, const SampleSet<T>& valid, int target_tag);
    // bool add_individual()
};

#endif  // INCLUDE_POPULATION_H_
