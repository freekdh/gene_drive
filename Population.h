#ifndef POPULATION_H // include guard
#define POPULATION_H
#include  <vector>
#include "random.h"
#include "Individual.h"

class Population{

    public:
    Population(std::vector<int> &initIndividuals);
    
    void mutate(std::vector<std::vector<int > > &mutationlist, const double &m);
    void drive(std::vector<std::vector<int > > &drivelist, const double &d);

    void selection(std::vector<double> &fitnesslist);
    void reproduce(std::vector<std::vector<int> > &recombinationlist);

    private:
    std::vector<Individual*>::iterator it;
    std::vector<Individual*> parents;
    std::vector<Individual*> offspring;

    double calculate_avgfitness(std::vector<double> fitnesslist);
};

#endif