#include <iostream>
#include <fstream>
#include <Math/Vector.h>
#include "Logger.h"
#include "vector"
#include "string"

/**
 * extracts the file name from the command line input
 */
std::string getName(int argc, char *argv[]) {
    logger.assert_always(argc>1,"Please provide root name (e.g. GCG1 as input argument");
    return std::string(argv[1]);
}

struct IOData {
    Mdouble aboveInflowLastTime = -1;
    Vec3D aboveInflowPosition;
};

void extractInOutData(std::string name, Mdouble inflowHeight, Mdouble outflowHeight) {
    std::ifstream dataFile(name+".data");
    std::ofstream inOutFile(name+".inOut");
    logger.assert(!dataFile.fail(),"data file could not be opened");
    logger.assert(!dataFile.fail(),"inOut file could not be opened");
    std::vector<IOData> data;

    unsigned n, species, index;
    Mdouble time, scalar;
    Vec3D position, vector;
    do {
        // read number of particles and time
        dataFile >> n;
        dataFile >> time;
        logger(INFO,"time % n % maxIndex %",time, n, data.size());
        dataFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
        for (unsigned i=0; i<n; ++i) {
            dataFile >> position >> vector >> scalar >> vector >> vector >> species >> index;
            //dataFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
            if (data.capacity()<index+1) data.reserve(2*index+1);
            if (data.size()<index+1) data.resize(index+1);
            IOData& val = data[index];
            if (position.Z>inflowHeight) {
                val.aboveInflowLastTime = time;
                val.aboveInflowPosition = position;
                //logger(INFO,"particle % flows in",index);
            } else if (position.Z<outflowHeight && val.aboveInflowLastTime!=-1) {
                inOutFile << index
                        << ' ' << species
                        << ' ' << val.aboveInflowLastTime
                        << ' ' << val.aboveInflowPosition
                        << ' ' << time
                        << ' ' << position << '\n';
                logger(INFO,"particle % flows out % > %",index,val.aboveInflowLastTime,time);
                val.aboveInflowLastTime = -1;
            }
            //std::cout << position << '\n';
        }
    } while (!dataFile.fail());
}

int main (int argc, char *argv[])
{
    std::string name = getName(argc, argv);
    Mdouble inflowHeight = 0.045;
    Mdouble outflowHeight = -0.045;
    extractInOutData(name,inflowHeight,outflowHeight);
    return 0;
}
