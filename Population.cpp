#include "Population.h"

Population::Population(std::vector<int> &initIndividuals){
        for(int i = 0; i < initIndividuals.size(); ++i){
            for(int j = 0; j < initIndividuals[i]; ++j){
                parents.push_back(new Individual(i));
            }
        }
}

