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

void WriteGameteToFile(Population &pop, std::vector<std::vector<int> > TypeToGametes, std::ofstream &out){
    std::vector<int> Gametes;
    pop.get_gametes(TypeToGametes,Gametes);
    const int nGametes = Gametes.size();
    for(int i = 0; i < nGametes; ++i){
        out << Gametes[i] << out.fill();
    }
    out << std::endl;
}

/* main */
int main(int argc, char *argv[]){
    rnd::set_seed();
    //hyperparameters and initial conditions
    const double r = 0.0;
    const int ngen = 10000;

    std::vector<double> Fitness = get_fitnesspars("Fitness.csv");
    std::vector<std::vector<int > > Mutation = get_mutationpars("Mutation.csv");
    std::vector<std::vector<int > > Drive = get_drivepars("MeioticDrive.csv");
    std::vector<std::vector<std::string> > Genotypes = get_genotypes("Genotypes.csv");
    std::vector<int> initIndividuals = get_initalindividuals("initIndividuals.csv");
    std::vector<std::vector<int> > Recombination = get_recombinationpars("Recombination.csv");
    std::vector<std::vector<int> > TypeToGametes = get_typetogametes("typetogametes.csv");
    set_RecombDistribution(0.5);
    
    std::ofstream out("test.txt");
    out.fill(',');
    Population initpop(initIndividuals);
    WriteGameteToFile(initpop,TypeToGametes,out);
    for(int i = 0; i < ngen; ++i){
        initpop.mutate(Mutation,0.1);
        initpop.drive(Drive,0.0);
        initpop.selection(Fitness);
        initpop.reproduce(Recombination);
        WriteGameteToFile(initpop,TypeToGametes,out);
    }

    return 0;
}