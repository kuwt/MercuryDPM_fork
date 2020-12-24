//Copyright (c) 2015, The MercuryDPM Developers Team. All rights reserved.
//For the list of developers, see <http://www.MercuryDPM.org/Team>.
//
//Redistribution and use in source and binary forms, with or without
//modification, are permitted provided that the following conditions are met:
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//  * Neither the name MercuryDPM nor the
//    names of its contributors may be used to endorse or promote products
//    derived from this software without specific prior written permission.
//
//THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
//ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
//WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
//DISCLAIMED. IN NO EVENT SHALL THE MERCURYDPM DEVELOPERS TEAM BE LIABLE FOR ANY
//DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
//(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
//LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
//ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
//(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
//SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

#include <iostream>
#include "Mercury3D.h"

int main(int argc, char *argv[])
{
    // write manual if number of arguments is not right
    const std::string manual = "Enter the prefix of the simulation you want to restart\n"
        " and, optionally, the prefix of the output files.\n"
        " e.g. restart2data input output\n"
        " reads from input.restart and writes to output.data\n";
    if (argc<2) logger(ERROR,manual);

    // first argument is used as prefix for the input
    std::string in = argv[1];
    //append .restart if necessary
    if (in.find(".restart")==-1) in += ".restart";
    logger(INFO,"Reading from %", in);

    // if second argument is given, use it as prefix for the output; otherwise, use input prefix
    bool prefixGiven = argc>=3 && argv[2][0]!='-';
    char* out = prefixGiven?argv[2]:argv[1];
    logger(INFO,"Writing to %", out);

    //read in from restart file and output data file
    Mercury3D problem;
    if (problem.readRestartFile(in)) {
        problem.setName(out);
        problem.cgHandler.computeContactPoints();
        problem.writeXBallsScript();
        problem.writeDataFile();
        problem.writeFStatFile();
        problem.writeEneFile();
        logger(INFO,"Written to %", problem.dataFile.getFullName());
    } else if (problem.readRestartFile(std::string(in)+".restart.0")
               || problem.readRestartFile(std::string(in)+".restart.0000")) {
        do {
            problem.setName(out);
            problem.cgHandler.computeContactPoints();
            problem.dataFile.setCounter(problem.restartFile.getCounter() - 1);
            problem.fStatFile.setCounter(problem.restartFile.getCounter() - 1);
            problem.eneFile.setCounter(problem.restartFile.getCounter() - 1);
            //problem.writeXBallsScript();
            problem.writeDataFile();
            problem.writeFStatFile();
            problem.writeEneFile();
            logger(INFO,"Written to %", problem.dataFile.getFullName());
            problem.setName(in);
        } while (problem.readRestartFile());
    } else {
        logger(ERROR,"File % not found",in);
    }
}
