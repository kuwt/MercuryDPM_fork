//Copyright (c) 2013-2020, The MercuryDPM Developers Team. All rights reserved.
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

#include "Mercury3D.h"
#include <CG/Fields/LiquidMigrationFields.h>
#include <CG/TimeAveragedCG.h>
#include "CG/CG.h"

int main(int argc, char *argv[])
{
    TimeAveragedCG<CGCoordinates::XZ,CGFunctions::Gauss,CGFields::LiquidMigrationFields> cg;
    //CG<CGCoordinates::O,CGFunctions::Gauss,CGFields::LiquidMigrationFields> cg;
    cg.setNZ(20);
    cg.setNX(5.31*cg.getNZ()); //to get the right ratio
    cg.setWidth(0.5);

    Mercury3D dpm;
    dpm.cgHandler.copyAndAddObject(cg);
    if (argc<=1) {
        logger(ERROR,"Please enter a filename as argument, e.g. ./CGForLiquids CSCRun_final.restart.0275");
    } else {
        dpm.cgHandler.restartAndEvaluateRestartFiles(argv[1]);
    }

    //write a m-file to read stat file into matlab and plot
    helpers::writeToFile(dpm.getName()+".m","%% read in stat file, reshape into 2D\n"
            "data = importdata('"+dpm.getName()+".0.stat',' ',2)\n"
            "nz = length(unique(data.data(:,3)))\n"
            "t = reshape(data.data(:,1),nz,[]);\n"
            "x = reshape(data.data(:,2),nz,[]);\n"
            "z = reshape(data.data(:,3),nz,[]);\n"
            "liquidBridgeVolume = reshape(data.data(:,4),nz,[]);\n"
            "liquidFilmVolume = reshape(data.data(:,5),nz,[]);\n"
            "%% plot\n"
            "subplot(2,1,1)\n"
            "contourf(x,z,liquidFilmVolume,20,'EdgeColor','none')\n"
            "axis image; title('liquidFilmVolume')\n"
            "subplot(2,1,2)\n"
            "contourf(x,z,liquidBridgeVolume,20,'EdgeColor','none')\n"
            "axis image; title('liquidBridgeVolume')");
    logger(INFO,"Run %.m to view output",dpm.getName());

    return 0;
}
