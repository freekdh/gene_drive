#ifndef POPULATION_H // include guard
#define POPULATION_H
#include  <vector>
#include "random.h"
#include "Individual.h"

class Population{

    public:
    Population(std::vector<int> &initIndividuals);
    
    private:
    std::vector<Individual*> parents;
    std::vector<Individual*> offspring;
};

#endif