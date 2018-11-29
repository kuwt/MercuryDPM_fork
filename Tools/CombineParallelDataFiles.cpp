#include <iostream>
#include <fstream>
#include "Logger.h"
#include "vector"
#include "string"

/**
 * extracts the file name from the command line input
 */
std::string getName(int argc, char *argv[]) {
    logger.assert_always(argc>1,"Please provide file root name as input argument");
    return std::string(argv[1]);
}


struct IFile {
    IFile (std::ifstream* file_, unsigned n_) : file(file_), n(n_) {}
    std::ifstream* file;
    unsigned n;
};

class DataFiles {
    std::vector<IFile> iFiles;
    std::ofstream oFile;
    unsigned oCount = 0;

public:

    /**
     * opens the input and output data files
     */
    DataFiles(std::string name) {
        // opens the input data files
        for(unsigned iCount = 0; true; ++iCount) {
            std::string fileName = name + ".data" + std::to_string(iCount);
            auto* file = new std::ifstream(fileName.c_str());
            if (file->fail()) {
                break;
            } else {
                iFiles.emplace_back(file,0);
            }
            logger(INFO,"Opened % for input",fileName);
        }

        // opens the output data file
        std::string oFileName = name + ".data";
        //logger.assert_always(!exist(oFileName),"File % already exists",oFileName);
        oFile.open(oFileName);
        if (oFile.fail()) logger(ERROR,"% could not be opened",oFileName);
        logger(INFO,"Opened % for output",oFileName);
    }

    bool exist(std::string fileName) {
        std::ifstream file(fileName);
        return !file.fail();
    }

    /**
     * closes the input/output data files
     */
    ~DataFiles() {
        // closes the input/output data files
        //logger(INFO, "Closing % input files", iFiles.size());
        for (auto iFile : iFiles) {
            delete iFile.file;
        }
    }

    /**
     * write all time steps
     */
    void write() {
        while (writeTimeStep()) {};
        logger(INFO,"Written % time steps",oCount);
    }

    /**
     * writes single timestep
     */
    bool writeTimeStep() {
        // read number of particles
        unsigned n = 0;
        for(auto& iFile : iFiles) {
            *iFile.file >> iFile.n;
            n += iFile.n;
            if (iFile.file->fail()) return false;
        }
        //read header
        std::string line;
        for(auto& iFile : iFiles) {
            std::getline(*iFile.file,line);
            //logger(INFO,"i % %",iFile.n,line);
        }
        //write header
        oFile << n << line << '\n';
        logger(INFO,"Writing %%",n,line);
        // write content
        for(auto& iFile : iFiles) {
            for (unsigned i=0; i<iFile.n; ++i) {
                std::getline(*iFile.file, line);
                oFile << line << '\n';
            }
        }
        ++oCount;
        return true;
    }

};

int main (int argc, char *argv[])
{
    //std::string name = "../../parallel/Drivers/USER/MercuryLab/FE/GCG14"; 
    std::string name = getName(argc, argv);
    DataFiles dataFiles(name);
    dataFiles.write();
    return 0;
}
