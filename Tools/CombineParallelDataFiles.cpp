#include <iostream>
#include <fstream>
#include <File.h>
#include <Math/Helpers.h>
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
    Mdouble time;
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
        oFileName = name + ".data";
        //logger.assert_always(!exist(oFileName),"File % already exists",oFileName);
        logger.assert(fileType!=FileType::NO_FILE,"File type cannot be NO_FILE");
        if (fileType==FileType::ONE_FILE) oFile.open(oFileName);
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

private:

    /**
     * writes single timestep
     */
    bool writeTimeStep() {
        // increase the file counter
        if (fileType!=FileType::ONE_FILE) {
            oFile.close();
            if (fileType==FileType::MULTIPLE_FILES_PADDED) {
                oFile.open(oFileName + '.' + to_string_padded(oCount));
            } else if (fileType==FileType::MULTIPLE_FILES) {
                oFile.open(oFileName + '.' + std::to_string(oCount));
            }
            if (oFile.fail()) logger(ERROR,"% could not be opened",oFileName);
        }

        // read header line
        static Mdouble time = 0;
        unsigned n = 0;
        std::string line;
        for(auto& iFile : iFiles) {
            // read number of particles
            *iFile.file >> iFile.n;
            n += iFile.n;
            // read time stamp
            *iFile.file >> iFile.time;
            //stop if this was last time step
            if (iFile.file->fail()) return false;
            //read rest of line
            std::getline(*iFile.file, line);
            //check this is indeed a header line
            static std::string headerLine = line;
            // this is the error-catching routine:
            // if this is not a true header line, or the timestamp is in the past, keep reading
            if (line != headerLine || iFile.time < time) {
                do {
                    *iFile.file >> iFile.n;
                    n += iFile.n;
                    *iFile.file >> iFile.time;
                    if (iFile.file->fail()) return false;
                    std::getline(*iFile.file, line);
                } while (line != headerLine || iFile.time < time);
                logger(WARN,"Mistake detected; moving forward to time %",iFile.time);
                oCount--; //remove last-written timestep (if MULTIPLE_FILES)
            }
            time = iFile.time;
        }
        time *= 1.000000001; //so the next time step is surely bigger than the last

        //debug output
        //for(auto& iFile : iFiles) std::cout << iFile.time << ' ';
        //std::cout << time << std::endl;

        //write header
        oFile << n << ' ' << time << line << '\n';
        logger(INFO,"Writing % %%",n,time,line);
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

    // base output file name (e.g. "problem.data")
    std::string oFileName;

public:
    // type of output file (one_file or multiple_padded)
    FileType fileType = FileType::ONE_FILE;
};

int main (int argc, char *argv[])
{
    //std::string name = "../../parallel/Drivers/USER/MercuryLab/FE/GCG14"; 
    std::string name = getName(argc, argv);
    DataFiles dataFiles(name);
    if (helpers::readFromCommandLine(argc,argv,"-multipleFiles"))
        dataFiles.fileType = FileType::MULTIPLE_FILES_PADDED;
    dataFiles.write();
    return 0;
}
