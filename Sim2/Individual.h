#ifndef INDIVIDUAL_H // include guard
#define INDIVIDUAL_H
#include<vector>
#include "random.h"

class Individual{
    public:
    Individual(const int &genotype);
    Individual(Individual *parent1, Individual *parent2);

    bool selection();
    void drive();
    int return_type1();
    int return_type2();

    private:
    bool drive1[3];
    bool drive2[3];
};

#endif