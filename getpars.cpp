#include <ofstream>
#include <vector>

std::vector<std::vector<std::string> > parseCSV(std::string &filename)
{
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

std::vector<std::vector<double> > get_fitnesspars(const std::string &name){
    std::vector<std::vector<std::string> > data = parseCSV(name);
    const int nrows = data.size();
    const int ncols = data[0].size();
    std::vector<std::vector<double> > out(nrows, std::vector<double>(ncols,0));
    for(int i = 0; i < nrows; ++i){
        for(int j = 0; j < ncols; ++j){
            out[i][j] = std::stod(data[i][j]);
        }
    }
    return out;
}

std::vector<std::vector<double> > get_mutationpars(const std::string &name){
    std::vector<std::vector<std::string> > data = parseCSV(name);
    const int nrows = data.size();
    const int ncols = data[0].size();
    std::vector<std::vector<double> > out(nrows, std::vector<double>(ncols,0));
    for(int i = 0; i < nrows; ++i){
        for(int j = 0; j < ncols; ++j){
            out[i][j] = std::stod(data[i][j]);
        }
    }
    return out;  
}