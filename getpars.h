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
std::vector<std::vector<double> > get_fitnesspars(const std::string &name);
std::vector<std::vector<double> > get_mutationpars(const std::string &name);
std::vector<std::vector<double> > get_drivepars(const std::string &name);
parameters get_initial(const std::string &name);

#endif
