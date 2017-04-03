/* Copyright 2017 Baytekov Nikita */
#include "GeneticDualizer.cpp"
#include "SampleHandler.cpp"
#include "LBBuilder.cpp"

template<class T>
class MONSClassifier {
    Vec< CollFamily<T> > Wsets;  // This is where all families of correctors are stored =)
    SampleHandler<T> sample_handler;
    LBBuilder<T> lb_builder;
    GeneticDualizer<T> genetic_dualizer;
    int max_iter;
    float eps;

 public:
    explicit MONSClassifier(int max_i = 1000, float epsilon = 0.0001):
            max_iter(max_i), eps(epsilon) {
    }

    ~MONSClassifier() {
    }

    void fit(const Mat<T> & X, const Vec<int> & y);
    Vec<int> predict(const Mat<T> & X);

 private:
    Vec< CollFamily<T> > & GeneticAlgorithm(const Mat<T> & Train, int class_tag, \
                                            const ElColl<T> & LB);
    bool CheckMargin(const Mat<T> & Valid);
};

template<class T>
void MONSClassifier<T>::fit(const Mat<T>& X, const Vec<int>& y) {
    // Mat<T> Train, Valid;
    // int class_num = GetTrainAndValidSets(X, Train, Valid);
    //
    // for (int i = 0; i < max_iter; i++) {
    //     for (int k = 0; k < class_num; k++) {
    //         Vec< ElClass_t<T1> > IterLB = BuildLB(Train, k);
    //         Vec_t< ElColl_t<T1> > IterWk = GeneticAlgorithm(Train, k, IterLB);
    //         Wsets[k] = UniteFamilies(Wsets[k], IterWk);
    //     }
    //
    //     if (CheckMargin(Valid))
    //         break;
    // }
}

int main(void) {
    return 0;
}
