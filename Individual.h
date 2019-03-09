#ifndef INDIVIDUAL_H // include guard
#define INDIVIDUAL_H
#include<vector>

const int ngenotypes = 265;
rnd::discrete_distribution chooserecombinant(16);

class Individual{
    public:
    Individual(const int&genotype);
    Individual::Individual(Individual &parent1, Individual &parent2, const double &r, std::vector<std::vector<int> > &Recombination);

    void mutate(std::vector<std::vector<int > > &mutationlist, const double &m);

    private:
    int type;
};

void set_RecombDistribution(const double &r){
    //Set chooserecombinant distribution based on r. This is due to the structure of Recombination.csv. See Mathematica file.
    for(int i = 0; i<4; ++i){
        chooserecombinant[i] = (1-r)*(1-r);
    }
    for(int i = 4; i<12; ++i){
        chooserecombinant[i] = (1-r)*r;
    }
    for(int i = 12; i<16; ++i){
        chooserecombinant[i] = r*r;
    }
}

#endif