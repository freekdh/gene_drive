#ifndef INDIVIDUAL_H // include guard
#define INDIVIDUAL_H
#include<vector>
#include "random.h"

class Individual{
    public:
    Individual(const int &genotype);
    Individual(Individual &parent1, Individual &parent2, std::vector<std::vector<int> > &Recombination);

    void mutate(std::vector<std::vector<int > > &mutationlist, const double &m);

    int return_type();
    private:
    int type;
};

void set_RecombDistribution(const double &r);


#endif