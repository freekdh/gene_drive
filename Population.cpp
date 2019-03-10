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
    parents = offspring;
    offspring.clear();
}

void Population::reproduce(std::vector<std::vector<int> > &recombinationlist){
    offspring.clear();
    const int totalparents = parents.size();
    const int totaloffspring = 100;
    for(int i = 0; i < totaloffspring; ++i){
        offspring.push_back(new Individual(parents[rnd::integer(totalparents)],parents[rnd::integer(totalparents)], recombinationlist));
    }

    for(it = parents.begin(); it < parents.end(); ++it){
        delete *it;
    }
    parents = offspring;
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

std::vector<int> Population::get_gametes(std::vector<std::vector<int> > &TypeToGametes, std::vector<int> &GameteCount){
    GameteCount.resize(16,0);
    for(it = parents.begin(); it < parents.end(); ++it){
        ++GameteCount[TypeToGametes[(*it)->return_type()][0]];
        ++GameteCount[TypeToGametes[(*it)->return_type()][1]];
    }
    return GameteCount;
}
