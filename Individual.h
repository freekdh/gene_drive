#ifndef INDIVIDUAL_H // include guard
#define INDIVIDUAL_H

#include<vector>

class Individual{
    public:
    Individual(const int&genotype);
    Individual::Individual(Individual &parent1, Individual &parent2, std::vector<std::vector<int> > &Recombination);

    void mutate(std::vector<std::vector<int > > &mutationlist, const double &m);

    private:
    int type;
};

#endif