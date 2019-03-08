#include "Individual.h"
#include "random.h"

Individual::Individual(const int&genotype):type(genotype){}

Individual::Individual(Individual &parent1, Individual &parent2, std::vector<std::vector<int> > &Recombination){
     const int type1 = parent1.type;
     const int type2 = parent2.type;
     const int index = (type1 - 1)*Length[Genotypes] + type2;
}

void Individual::mutate(std::vector<std::vector<int > > &mutationlist, const double &m){
    for(int i = 0; i < mutationlist[type].size(); ++i){
        type = rnd::uniform() <= m ? mutationlist[type][i] : type;
    }
}

