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
#include <fstream>
#include "Individual.h" 

// auto it = std::vector<Individual*>::iterator();

// class Population{
//     public:
//     Population(const int &n){
//         for(int i = 0; i < n; ++i){
//             population.push_back(new Individual(0));
//         }
//     }
//     void in_migration(Individual* migrant){
//         population.push_back(migrant);
//     }

//     void next(){
//         for(auto it = population.begin(); it <= population.end(); ++it){
            
//         }
//     }

//     private:
//     std::vector<Individual*> population;
// };

/* main */
int main(int argc, char *argv[]){
    //hyperparameters and initial conditions
    const double r = 0.5;
    std::vector<double> Fitness = get_fitnesspars("Fitness.csv");
    std::vector<std::vector<int > > Mutation = get_mutationpars("Mutation.csv");
    std::vector<std::vector<int > > Drive = get_drivepars("MeioticDrive.csv");
    std::vector<std::vector<std::string> > Genotypes = get_genotypes("Genotypes.csv");
    std::vector<int> initIndividuals = get_initalindividuals("initIndividuals.csv");
    std::vector<std::vector<int> > Recombination = get_recombinationpars("Recombination.csv");
    set_RecombDistribution(0.5);

    

    return 0;
}