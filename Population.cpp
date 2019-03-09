#include "Population.h"

Population::Population(std::vector<int> &initIndividuals){
    for(int i = 0; i < initIndividuals.size(); ++i){
        for(int j = 0; j < initIndividuals[i]; ++j){
            parents.push_back(new Individual(i));
        }
    }
}

void Population::mutate(std::vector<std::vector<int > > &mutationlist, const double &m){
    for(it = parents.begin(); it < parents.end(); ++it){
        (*it)->mutate(mutationlist,m);
    }
}

void Population::drive(std::vector<std::vector<int > > &drivelist, const double &d){
    for(it = parents.begin(); it < parents.end(); ++it){
        (*it)->drive(drivelist,d);
    }  
}

void Population::selection(std::vector<double> &fitnesslist){
    //double avg_fitness = calculate_avgfitness(fitnesslist); Assuming we have already normalized fitness.
    offspring.clear();
    for(it = parents.begin(); it < parents.end(); ++it){
        rnd::uniform() <= fitnesslist[(*it)->return_type()] ? offspring.push_back(*it) : delete *it;
    }
    parents.clear();
}

void Population::reproduce(std::vector<std::vector<int> > &recombinationlist){
    parents.clear();
    const int totaloffspring = offspring.size();
    const int totalparents = 100;
    for(int i = 0; i < totalparents; ++i){
        parents.push_back(new Individual(offspring[rnd::integer(totaloffspring)],offspring[rnd::integer(totaloffspring)], recombinationlist));
    }
    for(it = offspring.begin(); it < offspring.end(); ++it){
        delete *it;
    }
    offspring.clear();
}

double Population::calculate_avgfitness(std::vector<double> fitnesslist){
    // calculate average fitness
    double avg_fitness = 0.0;
    for(it = parents.begin(); it < parents.end(); ++it){
        avg_fitness += fitnesslist[(*it)->return_type()];
    }
    return avg_fitness/(double)parents.size();
}