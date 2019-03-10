#include "Individual.h"
#include "random.h"

const int ngenotypes = 265;
rnd::discrete_distribution chooserecombinant(16);

Individual::Individual(const int&genotype){type = genotype; }

Individual::Individual(Individual* parent1, Individual* parent2, std::vector<std::vector<int> > &Recombination){
     const int index = (parent1->type)*ngenotypes + parent2->type;
     type = Recombination[index][chooserecombinant.sample()]-1; // assuming recombinationrate = 0.5 -> Make a discrete distribtuion based on r.

}

void Individual::mutate(std::vector<std::vector<int > > &mutationlist, const double &m){
    for(int i = 0; i < mutationlist[type].size(); ++i){
        type = rnd::uniform() <= m ? mutationlist[type][i]-1 : type;
    }
}

void Individual::drive(std::vector<std::vector<int > > &drivelist, const double &d){
    for(int i = 0; i < drivelist[type].size(); ++i){
        type = rnd::uniform() <= d ? drivelist[type][i]-1 : type;
    }
}

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

int Individual::return_type() {return type;}
