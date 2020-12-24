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
    logger.assert_always(argc>1,"Please provide name of first data file as input argument");
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
    std::string ending = "";

public:

    /**
     * opens the input and output data files
     */
    DataFiles(int argc, char *argv[])
    {
        std::string name = getName(argc, argv);

        //check if there is an ending (the user might just specify the file root)
        size_t dot = name.find_last_of('.');
        // if no ending is found, assume .data0
        if (dot==std::string::npos) {
            name += ".data0";
        }

        //take off ending .0000 if necessary
        dot = name.find_last_of('.');
        char afterDot = name[dot + 1];
        //if last ending is a number
        if (afterDot >= '0' && afterDot <= '9') {
            ending = name.substr(dot);
            name.resize(dot);
        }

        // take off ending .data0
        {
            size_t dot = name.find_last_of('.');
            name.resize(dot);
        }

        // opens the input data files, if they exist
        for(unsigned processorID = 0; true; ++processorID) {
            //get input file name
            std::string fileName = name + ".data" + std::to_string(processorID) + ending;
            auto* file = new std::ifstream(fileName.c_str());
            if (file->fail()) {
                break;
            } else {
                iFiles.emplace_back(file,0);
                logger(INFO,"Opened % for input",fileName);
            }
        }
        logger.assert_always(!iFiles.empty(),"No input file found with name % ",name + ".data0" + ending);

        // opens output data file
        std::string oFileName = name + ".data" + ending;
        oFile.open(oFileName);
        if (oFile.fail()) {
            logger(ERROR,"% could not be opened",oFileName);
        } else {
            logger(INFO, "Opened % for output", oFileName);
        }

        write();
    }

    /**
     * closes the input/output data files
     */
    ~DataFiles() {
        for (auto iFile : iFiles) {
            delete iFile.file;
        }
    }

    /**
     * write all time steps
     */
    void write() {
        while (writeTimeStep()) {};
        if (oCount!=1) logger(INFO,"Written % time steps",oCount);
    }

private:

    /**
     * writes single timestep
     */
    bool writeTimeStep() {
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
};

int main (int argc, char *argv[])
{
    DataFiles dataFiles(argc,argv);
    return 0;
}

/*
 * Note: to merge muilti-file input data, use
   for i in 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15; do
       echo $i && ls -1 -t -r name.data$i.* | xargs cat > name.data$i;
   done
 */