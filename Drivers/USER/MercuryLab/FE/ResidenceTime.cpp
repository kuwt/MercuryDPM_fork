#include <iostream>
#include <fstream>
#include "Logger.h"
#include "vector"
#include "string"

/**
 * extracts the file name from the command line input
 */
std::string getDataName(int argc, char *argv[]) {
    logger.assert_always(argc>1,"Please provide data file name, output file name, inflow height and outflow height as input argument, e.g.\n  ./ResidenceTime GCG.data GCG.out14 0.035 -0.035");
    return std::string(argv[1]);
}

std::string getOutName(int argc, char *argv[]) {
    logger.assert_always(argc>2,"Please provide data file name, output file name, inflow height and outflow height as input argument, e.g.\n  ./ResidenceTime GCG.data GCG.out14 0.035 -0.035");
    return std::string(argv[2]);
}

double getInflowHeight(int argc, char *argv[]) {
    double height = (argc>3)?atof(argv[3]):0.04;
    logger(INFO,"Inflow height %", height);
    return height;
}

double getOutflowHeight(int argc, char *argv[]) {
    double height = (argc>4)?atof(argv[4]):-0.04;
    logger(INFO,"Inflow height %", height);
    return height;
}

struct Out {
    Mdouble time, id, radius;
};

std::vector<Out> readOut(std::string name) {
    std::vector<Out> OutData;
    std::ifstream file(name);
    logger.assert_always(!file.fail(),"% could not be opened",name);
    logger(INFO,"Opened % for output",name);
    //get header
    std::string line;
    std::getline(file,line);
    Out out;
    while (!file.fail()) {
        file >> out.time >> out.id >> out.radius;
        OutData.emplace_back(out);
    }
    logger(INFO,"Read in % output particles",OutData.size());
}

class ResidenceTime {
    std::ifstream dataFile;
    std::ofstream resFile;

public:

    /**
     * opens the input and output data files
     */
    ResidenceTime(int argc, char *argv[]) {
        std::string dataName = getDataName(argc,argv);
        std::string outName = getOutName(argc,argv);
        std::string resName = dataName.substr(0,dataName.length()-5)+".res";
        logger(INFO,"Output file: %",resName);

        double inflowHeight = getInflowHeight(argc,argv);
        double outflowHeight = getOutflowHeight(argc,argv);

        std::vector<Out> outflowTimes = readOut(outName);

        exit(-1);
    }
};

int main (int argc, char *argv[])
{
    ResidenceTime residenceTime(argc, argv);
    return 0;
}
