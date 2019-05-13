#ifndef GETPARS_H // include guard
#define GETPARS_H
#include <vector>
#include <string>
struct parameters{
    public:
        int n0;
        int n1;
        int n2;
    private:
};

std::vector<std::vector<std::string> > parseCSV(const std::string &filename);
std::vector<std::vector<std::string> > get_genotypes(const std::string &name);
std::vector<double> get_fitnesspars(const std::string &name);
std::vector<std::vector<int> > get_mutationpars(const std::string &name);
std::vector<std::vector<int > > get_drivepars(const std::string &name);
std::vector<int> get_initalindividuals(const std::string &name);
std::vector<std::vector<int> > get_recombinationpars(const std::string &name);
std::vector<std::vector<int> > get_typetogametes(const std::string &name);
#endif
