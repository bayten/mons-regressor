#include "default_types.h"
#include "containers.cpp"

class Individual
{

};

class Population
{
};

template<class T>
class GeneticInitiator
{
  public:
    virtual GeneticInitiator(T population);
    virtual ~GeneticInitiator();
    virtual Population get_init_population();
};

class GeneticCrossover
{

};

class GeneticMutator
{

};

class GeneticDualizer
{
  private:
    Population population;
    GeneticInitiator initiator;
    GeneticCrossover cross_op;
    GeneticMutator mutator;
};