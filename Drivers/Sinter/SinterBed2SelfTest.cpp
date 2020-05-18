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

#include "SphericalIndenter.h"
#include "Species/SinterReversibleAdhesiveSpecies.h"
#include "Species/SinterFrictionSpecies.h"
#include "Walls/InfiniteWall.h"
#include "Boundaries/PeriodicBoundary.h"
#include "Walls/AxisymmetricIntersectionOfWalls.h"
#include "Logger.h"
using std::cout;
using std::endl;
using helpers::readFromFile;
using helpers::to_string;
using helpers::writeToFile;

/// Single particle, indented slowly by spherical indenter.

class Indenter : public SphericalIndenter
{
public:

    Indenter(Mdouble indenterDiameter, Mdouble indentationVelocity, Mdouble indentationForce)
        : SphericalIndenter(indenterDiameter, indentationVelocity, indentationForce)
    {
        unsigned restartNumber = readFromFile<unsigned>("in","restartNumber",5);
        restartFile.setName("SinterBed1.restart."+to_string(restartNumber));
        readRestartFile();
        setRestarted(false);
        setName("SinterBed2");

        unsigned locationNumber = readFromFile<unsigned>("in","locationNumber",0);
        const unsigned nCells = 3; //number of grid cells in x and y
        Mdouble dx = (getXMax()-getXMin())/ static_cast<Mdouble>(nCells)*static_cast<Mdouble>(locationNumber%nCells);
        Mdouble dy = (getYMax()-getYMin())/ static_cast<Mdouble>(nCells)*static_cast<Mdouble>(locationNumber/nCells);
        Vec3D dPos(dx,dy,0);
        logger(INFO,"Shifting location by %",dPos);
        for (auto p : particleHandler) {
            p->setPosition(p->getPosition()+dPos);
        }

        species = dynamic_cast<SinterFrictionSpecies*>(speciesHandler.getObject(0));
        logger.assert(species,"Species pointer not set");
        species->setSinterRate(0.0);
        //species->setSinterAdhesion(0.0);
    }

    SinterFrictionSpecies *species;
};

int main(int argc UNUSED, char *argv[] UNUSED)
{
    Mdouble timeMax = 200;
    Mdouble indenterDiameter = 127e-6;
    Mdouble indentationForce = 4e-3; //4mN
    Mdouble indentationVelocity = readFromFile<Mdouble>("in","indentationVelocity",2e-6);

    //Ered 2e8
    //F=4e-3 = 4/3 E c del
    //F=4e-3 = 4/3 0.2e9 1e-6 1e-6  * 15
    //k=267
    Indenter i(indenterDiameter, indentationVelocity, indentationForce);
	if (i.getSystemDimensions()==2) {
        i.setIndentationForce(i.getIndentationForce()/3.0);
    }
    Mdouble frictionCoefficient = readFromFile<Mdouble>("in","frictionCoefficient",10);
    i.species->setSlidingFrictionCoefficient(frictionCoefficient);
    i.species->setRollingFrictionCoefficient(.1*frictionCoefficient);
    i.setTimeMax(timeMax);

    writeToFile("SinterBed2.gnu",
                         "set xlabel 'displacement [um]'\n"
                          "set ylabel 'force [mN]'\n"
                          "p 'SinterBed2.ene' u (-$2*1e6):($3*1e3) w lp\n"
    );

    i.setSaveCount(50);
    i.fStatFile.setFileType(FileType::NO_FILE);
    i.restartFile.setFileType(FileType::ONE_FILE);
    //i.setParticlesWriteVTK(true);
    i.solve();
    
    std::cout << "Execute 'gnuplot SinterBed2.gnu' to view output" << std::endl;

    return 0;
}
