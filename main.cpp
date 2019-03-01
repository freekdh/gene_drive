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
        haplotype1 = parent1.creategamete();
        haplotype2 = parent2.creategamete();
    }

    Individual operator+(Individual &mate){
        return Individual(*this, mate);
    }

    std::vector<int> creategamete(){
        return rnd::uniform() <= 0.5 ? haplotype1 : haplotype2; 
        }
    
    int return_first_locus(){
        return haplotype1[0];
    }

    private:
    std::vector<int> haplotype1;
    std::vector<int> haplotype2;
};

int main(int argc, char *argv[]){
    Individual a(1),b(2);
    Individual c = a+b;
    std::cout << c.return_first_locus() << std::endl;
    return 0;
}