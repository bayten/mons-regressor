/* Copyright 2017 Baytekov Nikita */

#include <iostream>
#include "./include/MONSClassifier.h"

int main(int argc, char* argv[]) {
    MONSClassifier<int, bool> mons_classifier(GeneticDualizer<int, bool>(5, 0.33, 3, 0.1));
    return 0;
}
