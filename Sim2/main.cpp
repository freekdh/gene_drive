/* 3 locus 3 allele IBS */
/* Freek de Haas, UBC */

#include <iostream>
#include <vector>
#include "random.h"
#include "utils.h"
#include <assert.h>
#include <fstream>
#include "Individual.h"

const int n_wild=1000;
const int n_eud=50;
const int n_gen=100;
const int n_rep = 1;
const double avg_offspring = 1.01;

std::vector<Individual*> initialize_pop(const int &n_wild, const int &n_eud);
void selection(std::vector<Individual*> &);
void drive(std::vector<Individual*> &);
void recombination(std::vector<Individual*> &);
void print(std::vector<Individual*> & , std::ofstream &);
void print_header(std::ofstream &);

/* main */
int main(int argc, char *argv[]){
    rnd::set_seed();
    std::ofstream output_file("./Data/sim_results.csv");
    output_file.fill(',');
    print_header(output_file);

    for(int rep = 0; rep < n_rep; ++rep){
        std::vector<Individual*> pop1 = initialize_pop(n_wild,n_eud);
        print(pop1,output_file);
        for(int gen = 0; gen < n_gen; ++gen){
            selection(pop1);
            drive(pop1);
            recombination(pop1);
            print(pop1, output_file);
        }
    }
    
    return 0;
}

std::vector<Individual*> initialize_pop(const int &n_wild, const int &n_eud){
    std::vector<Individual*> pop(n_wild+n_eud);
    for(int i = 0; i < n_wild; ++i){
        pop[i] = new Individual(0);
    }
    for(int i = n_wild; i < n_wild+n_eud; ++i){
        pop[i] = new Individual(1);
    }

    return pop;
}

void selection(std::vector<Individual*> &pop){
    std::vector<Individual*> offspring;
    const int pop_size = pop.size();
    for(int i = 0; i < pop_size; ++i){
       if(pop[i]->selection() == false) offspring.push_back(pop[i]);
    }
    pop = offspring;
}

void drive(std::vector<Individual*> &pop){
    const int pop_size = pop.size();
    for(int i = 0; i < pop_size; ++i){
        pop[i]->drive();
    }
}

void recombination(std::vector<Individual*> &pop){
    const int pop_size = pop.size();
    const int n_offspring = rnd::poisson(avg_offspring*(double)pop_size);
    std::vector<Individual*> offspring(n_offspring);
    for(int i = 0; i < n_offspring; ++i){
        offspring[i] = new Individual(pop[rnd::integer(pop_size)],pop[rnd::integer(pop_size)]);
    }

    pop=offspring;
}

void print(std::vector<Individual*> &pop, std::ofstream &output_file){
    const int pop_size = pop.size();
    std::vector<int> count_types(8,0);
    for(int i = 0; i < pop_size; ++i){
        ++count_types[pop[i]->return_type1()];
        ++count_types[pop[i]->return_type2()];
    }

    for(int i = 0; i < count_types.size(); ++i){
        output_file << count_types[i] << output_file.fill();
    }
    output_file << '\n';
}

void print_header(std::ofstream &output_file){
    for(int i = 0; i < 8; ++i){
        output_file << "type" << i << output_file.fill();
    }
    output_file << '\n';
}