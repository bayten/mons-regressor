#include <iostream>
#include <set>
#include <vector>

#include
#include "default_types.h"


template<class T1, class T2>  // T1 - type of features, T2 - type of class tags
class MONSClassifier
{
    vector<Family> Wsets;  // This is where all families of correctors are stored =)
    GeneticDuali
    int max_iter;
    float eps;

  public:
    MONSClassifier(int max_i=1000, float epsilon=0.0001) : max_iter(max_i), eps(epsilon) {};
    ~MONSClassifier(){};

    void fit(const Mat_t<T1> & X, const Vec_t<T2> & y);
    Vec_t<T2> predict(const Mat_t<T1> & X);

  private:
    int GetTrainAndValidSets(const Mat_t<T1> & X, const Vec_t<T2> & y,\
                             SampleSet_t<T1, T2> & Train,\
                             SampleSet_t<T1, T2> & Valid);
    ElColl_t<T1> & BuildLB(const SampleSet_t<T1, T2> & Train, T2 class_tag);  // TODO: BOOSTLO
    Vec_t< ElColl_t<T1> > & GeneticAlgorithm(const Mat_t<T1> & Train, T2 class_tag,\
                                             const ElColl_t<T1> & LB);
    Vec_t< ElColl_t<T1> > & UniteFamilies(const Vec_t< ElColl_t<T1> > & W1,\
                                          const Vec_t< ElColl_t<T1> > & W2);
    bool CheckMargin(const Mat_t<T1> & Valid);
};

template<class T1, class T2>
void MONSClassifier<T1, T2>::fit(const Mat_t<T1> & X, const Vec_t<T2> & y)
{
    Mat_t<T1> Train, Valid;
    int class_num = GetTrainAndValidSets(X, Train, Valid);

    for(int i = 0; i < max_iter; i++){
        for(int k = 0; k < class_num; k++){
            Vec_t< ElClass_t<T1> > IterLB = BuildLB(Train, k);
            Vec_t< ElColl_t<T1> > IterWk = GeneticAlgorithm(Train, k, IterLB);
            Wsets[k] = UniteFamilies(Wsets[k], IterWk);
        }

        if(CheckMargin(Valid))
            break;
    }
}

template<class T1, class T2>
int MONSClassifier<T1, T2>::GetTrainAndValidSets(const Mat_t<T1> & X, const Vec_t<T2> & y,\
                                                 SampleSet_t<T1, T2> & Train,\
                                                 SampleSet_t<T1, T2> & Valid)
{
    Mat_t<T1> tornX;
    Vec_t<T2> torny;

    Train.append(X, y);
    Train.shuffle();

    int class_num = Train.samples.sz;  // doing this to have samples of every
    for(int i = 0; i < class_num; i++){  // class in validation set
        Train.samples[i].tearup(2*int(log(Train[i].sz, 10))+1, tornX, torny);
        Valid.append(tornX, torny);
    }

    return class_num;
}

template<class T1, class T2>
ElColl_t<T1> & MONSClassifier<T1, T2>::BuildLB(const SampleSet_t<T1, T2> & Train, T2 class_tag)
{
    // Generate all elementary classifiers for this class
    // Get random subset 2*log2(m)+3


}
    

Vec_t< ElColl_t<T1> > & GeneticAlgorithm(const Mat_t<T1> & Train, T2 class_tag,\
                                         const ElColl_t<T1> & LB)
{
    
}


int main(void)
{
    return 0;
}