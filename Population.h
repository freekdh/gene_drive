#ifndef POPULATION_H 
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

    std::vector<int> get_gametes(std::vector<std::vector<int> > &TypeToGametes, std::vector<int> &GameteCount);
    private:
    std::vector<Individual*>::iterator it;
    std::vector<Individual*> parents;
    std::vector<Individual*> offspring;

    double calculate_avgfitness(std::vector<double> fitnesslist);
};

#endif