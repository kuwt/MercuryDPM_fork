//Copyright (c) 2013-2023, The MercuryDPM Developers Team. All rights reserved.
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

#include <Mercury3D.h>
#include <CG/TimeSmoothedCG.h>
#include "CG/TimeAveragedCG.h"

int main(int argc, char** argv)
{
    //determine restart file
    std::string restartFile;
    Mdouble dt = 1e-3;
    if (argc<2) {
        restartFile = "N25300/SiloBehzadMuR0.1.restart";
        logger(WARN,"Loading from default restart file %.\n"
                "You can also specify the restart file, e.g. ./SiloBehzadCG SiloBehzad.restart", restartFile);
        logger(WARN,"Using default time interval %", dt);
    } else if (argc<3) {
        restartFile = argv[1];
        logger(WARN,"Using default time interval %", dt);
    } else {
        restartFile = argv[1];
        dt = atof(argv[2]);
    }

    //restart existing simulation
    Mercury3D silo;
    silo.readRestartFile(restartFile);
    silo.setName(silo.getName()+"CG");
    silo.setFileType(FileType::NO_FILE);

    //create xz statistics
    {
        TimeAveragedCG<CGCoordinates::XZ, CGFunctions::Heaviside> c;
        c.setN({25, 1, 100});
        c.setWidth(3e-3);
        c.statFile.setSaveCount(25);
        c.statFile.setName(silo.getName() + "XZ.TA.stat");
        silo.cgHandler.copyAndAddObject(c);
    }

    //create xz statistics
    {
        TimeSmoothedCG<CGCoordinates::XZ, CGFunctions::Heaviside> c;
        c.setN({25, 1, 100});
        c.setWidth(3e-3);
        c.statFile.setSaveCount(25);
        c.setWidthTime(0.01);
        c.setTimeStep(0.005);
        c.statFile.setName(silo.getName() + ".XZ.TS.stat");
        silo.cgHandler.copyAndAddObject(c);
    }

    //create O statistics
    {
        CG<CGCoordinates::O> c;
        c.statFile.setSaveCount(25);
        c.statFile.setName(silo.getName() + "O.T.stat");
        silo.cgHandler.copyAndAddObject(c);
    }

    //create O statistics
    {
        TimeSmoothedCG<CGCoordinates::O> c;
        c.statFile.setSaveCount(25);
        c.setWidthTime(0.01);
        c.setTimeStep(0.005);
        c.statFile.setName(silo.getName() + ".O.TS.stat");
        silo.cgHandler.copyAndAddObject(c);
    }

    //continue run for 1ms and take statistics
    silo.setTimeMax(silo.getTime()+dt);
    silo.solve(argc-1, argv+1);
    return 0;
}
