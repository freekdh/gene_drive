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
    std::vector<std::vector<int> > TypeToGametes = get_typetogametes("typetogametes.csv");
    set_RecombDistribution(0.5);

    for(int i = 0; i < initIndividuals.size(); ++i){
        std::cout << initIndividuals[i] << std::endl;
    }

    std::vector<int> count1, count2, count3, count4;

    Population initpop(initIndividuals);
    count1 = initpop.get_gametes(TypeToGametes);
    initpop.mutate(Mutation,0.01);
    initpop.drive(Drive,0.99);
    initpop.selection(Fitness);
    initpop.reproduce(Recombination);
    count2 = initpop.get_gametes(TypeToGametes);
    for(int i = 0; i < count1.size(); ++i){
        std::cout << count1[i] << '\t' << count2[i] << '\n'; 
    }
    std::cout << "test";
    initpop.mutate(Mutation,0.01);
    std::cout << "test";
    initpop.drive(Drive,0.99);
    std::cout << "test";
    initpop.selection(Fitness);
    std::cout << "test";
    initpop.reproduce(Recombination);
    std::cout << "test";
    count3 = initpop.get_gametes(TypeToGametes);
    std::cout << "test";
    initpop.mutate(Mutation,0.01);
    initpop.drive(Drive,0.99);
    initpop.selection(Fitness);
    initpop.reproduce(Recombination);
    count4 = initpop.get_gametes(TypeToGametes);


    return 0;
}