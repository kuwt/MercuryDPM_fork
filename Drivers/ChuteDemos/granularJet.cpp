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

#include <sstream>
#include <iostream>
#include <iomanip>
#include <cmath>

#include "ChuteWithHopper.h"

#include <sys/types.h>
#include <sys/stat.h>
#include <Species/LinearViscoelasticFrictionSpecies.h>
#include <cstring>

using namespace std;

class GranularJet : public ChuteWithHopper
{
public:
    
    void setSilbert()
    {
        //time stepping
        setTimeStep(1e-4);
        setTimeMax(1e20);
        setSaveCount(1e4); //save every unit time (\hat{t}=sqrt(d/g)=~0.008s)
        
        //particle radii
        setInflowParticleRadius(.5);
        setFixedParticleRadius(.5);//getInflowParticleRadius());
        setRoughBottomType(MULTILAYER);
        
        //particle properties
        LinearViscoelasticFrictionSpecies species;
        species.setDensity(6 / constants::pi);
        species.setStiffness(2e5);
        species.setSlidingStiffness(2.0 / 7.0 * species.getStiffness());
        species.setDissipation(25.0);
        species.setSlidingDissipation(species.getDissipation());
        species.setSlidingFrictionCoefficient(0.5);
        speciesHandler.copyAndAddObject(species);
        
        //chute properties
        setChuteAngleAndMagnitudeOfGravity(0.0, 1.0);
    }
    
    bool readNextArgument(int& i, int argc, char* argv[]) override
    {
        if (!strcmp(argv[i], "-ExitLength"))
        {
            Mdouble exitHeight = atof(argv[i + 1]);
            Mdouble exitLength = 1.0 * exitHeight;
            Mdouble hopperAngle = 45.0;
            Mdouble hopperLength = 4.0 * exitLength;
            Mdouble hopperHeight = hopperLength;
            setHopper(exitLength, exitHeight, hopperAngle, hopperLength, hopperHeight);
        }
        else
        {
            return Chute::readNextArgument(i, argc, argv);
        }
        return true; //returns true if argv is found
    }
    
    void actionsBeforeTimeLoop() override
    {
        write(std::cout, false);
        writeRestartFile();
    }
    
    
};

int main(int argc, char* argv[])
{
    GranularJet problem;
    
    // Particle properties
    problem.setSilbert();
    // Problem parameters
    problem.setName("GranularJet");
    
    //Corrections
    problem.setRoughBottomType(MONOLAYER_DISORDERED);
    
    // Chute properties
    problem.setChuteAngle(20);
    problem.setChuteLength(700 * .25); //700=40cm
    problem.setChuteWidth(400 * .25); //400=24cm
    problem.setMaxFailed(6);
    problem.makeChutePeriodic();
    problem.setHopperDimension(2);
    problem.setHopperLift(200); //250=15cm
    Mdouble exitHeight = 25;
    Mdouble exitLength = 3.0 * exitHeight;
    Mdouble hopperAngle = 45.0;
    Mdouble hopperLength = 4.0 * exitLength;
    Mdouble hopperHeight = hopperLength;
    problem.setHopper(exitLength, exitHeight, hopperAngle, hopperLength, hopperHeight);
    problem.setZMax(200);
    
    problem.dataFile.setFileType(FileType::NO_FILE);
    problem.fStatFile.setFileType(FileType::NO_FILE);
    problem.restartFile.setFileType(FileType::MULTIPLE_FILES_PADDED);
    
    //solve
    problem.readArguments(argc, argv);
    problem.solve();
}
