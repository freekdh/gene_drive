#include "Individual.h"
#include "random.h"

Individual::Individual(const int&genotype){
    if(genotype==0){
        drive1[0] = 0;
        drive1[1] = 0;
        drive1[2] = 0;
        drive2[0] = 0;
        drive2[1] = 0;
        drive2[2] = 0;
    }
    else{
        drive1[0] = 1;
        drive1[1] = 1;
        drive1[2] = 1;
        drive2[0] = 1;
        drive2[1] = 1;
        drive2[2] = 1;
    }
}

Individual::Individual(Individual* parent1, Individual* parent2){
     drive1[0] = rnd::uniform() <= 0.5 ? parent1->drive1[0] : parent1->drive2[0];
     drive1[1] = rnd::uniform() <= 0.5 ? parent1->drive1[1] : parent1->drive2[1];
     drive1[2] = rnd::uniform() <= 0.5 ? parent1->drive1[2] : parent1->drive2[2];

     drive2[0] = rnd::uniform() <= 0.5 ? parent2->drive1[0] : parent2->drive2[0];
     drive2[1] = rnd::uniform() <= 0.5 ? parent2->drive1[1] : parent2->drive2[1];
     drive2[2] = rnd::uniform() <= 0.5 ? parent2->drive1[2] : parent2->drive2[2];
}

void Individual::drive(){
    if(drive1[0] || drive2[0]){
        if(drive1[1] || drive2[1]){
            drive1[1]=true;
            drive2[1]=true;
        }
    }

    if(drive1[1] || drive2[1]){
        if(drive1[2] || drive2[2]){
            drive1[2]=true;
            drive2[2]=true;
        }
    }
}

bool Individual::selection(){
    const int n_eud = (int)drive1[0]+(int)drive1[1]+(int)drive1[2]+(int)drive2[0]+(int)drive2[1]+(int)drive2[2];
    const double s = pow((1-0.1),double(n_eud));
    return rnd::uniform() < s ? false : true;
}

int Individual::return_type1(){
    return drive1[0]*1+drive1[1]*2+drive1[2]*4;
}

int Individual::return_type2(){
    return drive2[0]*1+drive2[1]*2+drive2[2]*4;
}

