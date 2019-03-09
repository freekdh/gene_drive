/* 3 locus 3 allele IBS */
/* Freek de Haas, UBC */

#include <bitset>
#include <iostream>
#include <vector>
#include "random.h"
#include "utils.h"
#include "getpars.h"
#include <assert.h>
#include <fstream>
#include "Individual.h"
#include "Population.h"

/* main */
int main(int argc, char *argv[]){
    rnd::set_seed();
    //hyperparameters and initial conditions
    const double r = 0.5;
    std::vector<double> Fitness = get_fitnesspars("Fitness.csv");
    std::vector<std::vector<int > > Mutation = get_mutationpars("Mutation.csv");
    std::vector<std::vector<int > > Drive = get_drivepars("MeioticDrive.csv");
    std::vector<std::vector<std::string> > Genotypes = get_genotypes("Genotypes.csv");
    std::vector<int> initIndividuals = get_initalindividuals("initIndividuals.csv");
    std::vector<std::vector<int> > Recombination = get_recombinationpars("Recombination.csv");
    set_RecombDistribution(0.5);

    Population initpop(initIndividuals);
    initpop.mutate(Mutation,0.01);
    initpop.drive(Drive,0.99);
    initpop.selection(Fitness);
    initpop.reproduce(Recombination);
    return 0;
}