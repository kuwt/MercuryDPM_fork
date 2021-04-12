//Copyright (c) 2013-2021, The MercuryDPM Developers Team. All rights reserved.
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
#include "MercuryOS.h"
#include <Walls/InfiniteWall.h>

/**
 * Tests tangential forces in particle starting to roll over wall.
 */
class MindlinSelfTest : public MercuryOS {

public:
    
    // used to set the initial conditions of the particles, walls, species, etc
    void setupInitialConditions () override
    {
        // name
        setName("OSMindlinSelfTest");

        // gravity
        setGravity({5, 0, -5});

        // time step and maximum simulation time
        setTimeStep(5e-7);
        setTimeMax(5e-3);

        // output
        setSaveCount(10);
        //fStatFile.writeFirstAndLastTimeStep();
        restartFile.writeFirstAndLastTimeStep();

        // domain for visualisation
        setMax(Vec3D(1e-3,1e-3,2.2e-3));
        setMin(Vec3D(-1e-3,-1e-3,0));
    
        // define the material properties of M1, M2, steel (see MercuryOS.h)
        setMaterialProperties();

        // add particles
        SphericalParticle particle;
        particle.setSpecies(steel);
        particle.setRadius(0.001);
        particle.setPosition(Vec3D(0,0,0.001));
        particleHandler.copyAndAddObject(particle);

        // add walls
        InfiniteWall wall;
        wall.setSpecies(m1);
        wall.set(Vec3D(0,0,-1),Vec3D(0,0,0));
        wallHandler.copyAndAddObject(wall);
    }
};

int main(int argc, char** argv)
{
    // create an instance of the class
    MindlinSelfTest dpm;
    // call the solve routine
    dpm.solve();
    // create analysis script
    helpers::writeToFile("OSMindlinSelfTest.gnu",
                         "set xlabel 'time [s]'\n"
                         "set ylabel 'force [N]'\n"
                         "p 'OSMindlinSelfTest.fstat' u 1:9 t 'normal', '' u 1:10 t 'tangential'\n");
    logger(INFO,"Run 'gnuplot OSMindlinSelfTest.gnu --persist' to show resulting forces");
    return 0;
}
