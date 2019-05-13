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
    int return_type();

    private:
    int type;
};

#endif