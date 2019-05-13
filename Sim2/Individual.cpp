#include "Individual.h"
#include "random.h"

Individual::Individual(const int&genotype){type = genotype;}

Individual::Individual(Individual* parent1, Individual* parent2){
     type=parent1->type;
}

void Individual::drive(){
    int i =0;
}

bool Individual::selection(){
    return false;
}

int Individual::return_type(){
    return 1;
}

