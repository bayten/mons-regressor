/* Copyright 2017 Baytekov Nikita */
#include "GeneticDualizer.cpp"
#include "SampleHandler.cpp"
#include "LBBuilder.cpp"

template<typename S, typename T>
class MONSClassifier {
    Vec< CollFamily<T> > coll_sets;
    SampleHandler<T> sample_handler;
    LBBuilder<T> lb_builder;
    GeneticDualizer<S, T> genetic_dualizer;

    int max_iter;
    float eps;

 public:
    explicit MONSClassifier(int max_i = 1000, float epsilon = 0.0001);
    ~MONSClassifier() {}

    void fit(const Mat<T>& X, const Vec<int>& y);
    Vec<int> predict(const Mat<T> & X);

 private:
    bool check_margin();
};

template<typename S, typename T>
void MONSClassifier<S, T>::fit(const Mat<T>& X, const Vec<int>& y) {
    // Mat<T> train, valid;
    // int class_num = sampleGetTrainAndValidSets(X, Train, Valid);
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
