/* 3 locus 3 allele IBS */
/* Freek de Haas, UBC */

/* headers */
#include <bitset>
#include <iostream>
#include <vector>
#include "random.h"
#include "utils.h"
#include "getpars.h"
#include <assert.h>

/*global variables*/
std::vector<std::vector<double> > pars_fitness;
std::vector<std::vector<double> > pars_mutation;
std::vector<std::vector<double> > pars_drive;
std::vector<double> recombinationrate;

enum allele_type {variant0, variant1, variant2};

class Individual{
    public:
    Individual(Individual &parent1, Individual &parent2){
        haplotype[0] = parent1.creategamete();
        haplotype[1] = parent2.creategamete();
    }
    Individual operator+(Individual &mate){
        return Individual(*this, mate);
    }

    std::vector<int> creategamete(){
        std::vector<int> gamete;
        bool focal = rnd::bernoulli(0.5);
        for(int i = 0; i < 3; ++i){
            gamete[i] = haplotype[focal][i];
            focal = rnd::uniform() <= recombinationrate[i] ?  !focal : focal;
        }
        return gamete;
    }

    void mutate(){
        
    }

    double fitness(){
        return 1.0;
    }

    private:
    std::vector<int> haplotype[2];
};

auto it = std::vector<Individual*>::iterator();

class Population{
    public:
    Population(const int &n){
        for(int i = 0; i < n; ++i){
            population.push_back(new Individual(0));
        }
    }
    void in_migration(Individual* migrant){
        population.push_back(migrant);
    }

    void next(){
        for(auto it = population.begin(); it <= population.end(); ++it){
            
        }
    }

    private:
    std::vector<Individual*> population;
};

/* main */
int main(int argc, char *argv[]){
    //hyperparameters
     pars_fitness = get_fitnesspars("fitness.csv");
     pars_mutation = get_mutationpars("mutation.csv");
     pars_drive = get_drivepars("drive.csv");
    
    //initial conditions
    
    Individual a(1),b(2);
    Individual c = a+b;

    Population aa(10);
    Population bb(10);
    return 0;
}