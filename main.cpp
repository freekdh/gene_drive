#include <bitset>
#include <iostream>
#include <vector>
#include "random.h"
#include "utils.h"

enum individual_constructor {randomi, wildtype, engineered};
enum allele_type {wild, engin1, engin2};

class Individual{
    public:
    Individual(int type){
        switch(type) {  
            case wildtype : 
                haplotype1.resize(3,wild);
                haplotype2.resize(3,wild);
                    break;
            case engineered :
                haplotype1.resize(3,engin1);
                haplotype2.resize(3,engin2);
                    break;
        }
    };
    Individual(Individual &parent1, Individual &parent2){
        haplotype1 = parent1.creategamete(0.1);
        haplotype2 = parent2.creategamete(0.1);
    }

    Individual operator+(Individual &mate){
        return Individual(*this, mate);
    }

    std::vector<int> creategamete(const double &drive){
        return haplotype1;
        //rnd::uniform() <= 0.5 ? haplotype1 : haplotype2; 
    }

    void mutate(){
        for(int i = 0; i < haplotype1.size(); ++i){
            haplotype1[i] = rnd::uniform() <= 0.5 ? haplotype1[i] : haplotype1[i];
            haplotype2[i] = rnd::uniform() <= 0.5 ? haplotype2[i] : haplotype2[i];
        }
    }

    double fitness(){
        return 1.0;
    }

    private:
    std::vector<int> haplotype1;
    std::vector<int> haplotype2;
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


int main(int argc, char *argv[]){
    //pars
    std::vector<std::vector<double>> fitness = get_fitnesspars("fitness.csv");
    std::vector<double> mutation = get_mutationpars("mutation.csv");
    std::vector<double> drive = get_drivepars("drive.csv");
    

    Individual a(1),b(2);
    Individual c = a+b;

    Population aa(10);
    Population bb(10);
    return 0;
}