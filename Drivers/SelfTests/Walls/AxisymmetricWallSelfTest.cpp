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

#include <Species/LinearViscoelasticSpecies.h>
#include <Mercury3D.h>
#include "Walls/AxisymmetricIntersectionOfWalls.h"

/**
 * This class defines all properties ot the THZ 1125G Mixer; see THZ-BROCHURE.PDF
 * It then adds particles.
 */
class AxisymmetricWallSelfTest : public Mercury3D
{
public:

    /*
     * Things that should not be done when restarting, leave in setupInitialConditions
     */
    AxisymmetricWallSelfTest()
    {
        //set domain
        setDomain({-1,-1,-1},{1,1,1});
        //turn on VTK output
        setParticlesWriteVTK(true);
        setWallsWriteVTK(FileType::MULTIPLE_FILES);
        // set file outout
        fStatFile.setFileType(FileType::NO_FILE);
        //define species and time step
        setSpeciesAndTimeStep();
        //define mixer geometry
        setGeometry();
        // set process time
        setTimeMax(0.0);
        // set name
        setName("AxisymmetricWallSelfTest");
    }

    void setupInitialConditions() override
    {
        wallHandler.getLastObject()->addParticlesAtWall(100);
    }

    //create a species for all walls and particles
    void setSpeciesAndTimeStep()
    {
        LinearViscoelasticSpecies s;
        s.setCollisionTimeAndRestitutionCoefficient(1.0, 1.0, 1.0);
        auto species = speciesHandler.copyAndAddObject(s);
        setTimeStep(0.1);
    }

    //define base wall and add to wallHandler
    void setGeometry()
    {
        ParticleSpecies *wallSpecies = speciesHandler.getLastObject();
        AxisymmetricIntersectionOfWalls w;
        w.setSpecies(wallSpecies);
        w.setPosition({0,0,0});
        w.createPrism({{.1,0,.1}, {.1,0,.2}, {.2,0,.1}}, {0,1,0});
        w.setAxis({0,0,1});
        wallHandler.copyAndAddObject(w);
        logger(INFO, "Added % walls", wallHandler.getSize());
    }
};


/*
 * Solves the default Orange mixer aplication (5-10 min's runtime)
 */
int main(int argc, char** argv)
{
    AxisymmetricWallSelfTest dpm;
    dpm.solve();
    return 0;
}
