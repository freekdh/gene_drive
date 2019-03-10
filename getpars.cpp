#include <fstream>
#include <vector>
#include <sstream>
#include <iostream>
#include "getpars.h"

std::vector<std::string> variablenames;

std::vector<std::vector<std::string> > parseCSV(const std::string &filename){
    std::ifstream  data(filename);
    std::string line;
    std::vector<std::vector<std::string> > parsedCsv;
    while(std::getline(data,line))
    {
        std::stringstream lineStream(line);
        std::string cell;
        std::vector<std::string> parsedRow;
        while(std::getline(lineStream,cell,','))
        {
            parsedRow.push_back(cell);
        }
        parsedCsv.push_back(parsedRow);
    }
    return parsedCsv;
};

std::vector<double> get_fitnesspars(const std::string &name){
    std::vector<std::vector<std::string> > data = parseCSV(name);
    const int nrows = data.size();
    std::vector<double> out(nrows, 0);
    for(int i = 0; i < nrows; ++i){
            out[i] = std::stod(data[i][0]);
    }
    return out;
}

std::vector<std::vector<std::string> > get_genotypes(const std::string &name){
    return parseCSV(name);
}

std::vector<std::vector<int> > get_mutationpars(const std::string &name){
    std::vector<std::vector<std::string> > data = parseCSV(name);
    const int nrows = data.size();
    std::vector<std::vector<int> > out(nrows);
    for(int i = 0; i < nrows; ++i){
        const int ncols = data[i].size();
        std::vector<int> temp(ncols);
        for(int j = 0; j < ncols; ++j){
            temp[j] = std::stoi(data[i][j]);
        }
        out[i]=temp;
    }
    return out;
}

std::vector<std::vector<int > > get_drivepars(const std::string &name){
    std::vector<std::vector<std::string> > data = parseCSV(name);
    const int nrows = data.size();
    std::vector<std::vector<int> > out(nrows);
    for(int i = 0; i < nrows; ++i){
        const int ncols = data[i].size();
        std::vector<int> temp(ncols);
        for(int j = 0; j < ncols; ++j){
            temp[j] = std::stoi(data[i][j]);
        }
        out[i]=temp;
    }
    return out;
}

std::vector<std::vector<int> > get_recombinationpars(const std::string &name){
    std::vector<std::vector<std::string> > data = parseCSV(name);
    const int nrows = data.size();
    std::vector<std::vector<int> > out(nrows);
    for(int i = 0; i < nrows; ++i){
        const int ncols = data[i].size();
        std::vector<int> temp(ncols);
        for(int j = 0; j < ncols; ++j){
            temp[j] = std::stoi(data[i][j]);
        }
        out[i]=temp;
    }
    return out;
}

std::vector<int> get_initalindividuals(const std::string &name){
    std::vector<std::vector<std::string> > data = parseCSV(name);
    const int nrows = data.size();
    std::vector<int> out(nrows, 0);
    for(int i = 0; i < nrows; ++i){
            out[i] = std::stoi(data[i][0]);
    }
    return out;
}

std::vector<std::vector<int> > get_typetogametes(const std::string &name){
    std::vector<std::vector<std::string> > data = parseCSV(name);
    const int nrows = data.size();
    std::vector<std::vector<int> > out(nrows);
    for(int i = 0; i < nrows; ++i){
        const int ncols = data[i].size();
        std::vector<int> temp(ncols);
        for(int j = 0; j < ncols; ++j){
            temp[j] = std::stoi(data[i][j]);
        }
        out[i]=temp;
    }
    return out;
}
