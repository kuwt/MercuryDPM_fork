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

#include <iostream>
#include "Mercury3D.h"
#include "Boundaries/CubeInsertionBoundary.h"
#include "Species/LinearViscoelasticSpecies.h"


/**
 * An example code on how to add particles of a given particle size distribution and flow rate.
 */
class PSDSelfTest : public Mercury3D
{
public:

    /**
     * Define an insertion boundary, which will insert particles in a given region
     */
    void setupInitialConditions() override
    {
        setName("UniformPSDSelfTest");
        // make volume 2 times bigger to get a good test
        setDomain(Vec3D(0, 0, 0), Vec3D(.5, .5, .5));
        setTimeStep(1);
        setTimeMax(0);
        
        LinearViscoelasticSpecies species;
        species.setDensity(2000);
        species.setStiffness(10000);
        auto s = speciesHandler.copyAndAddObject(species);
        
        // define a cube insertion boundary, which inserts particles in a cuboid region
        CubeInsertionBoundary insertionBoundary;
        SphericalParticle particle(s);
        insertionBoundary.set(&particle, 10, getMin(), getMax(), {0,0,0}, {0,0,0}, 1, 1);
        insertionBoundary.setInitialVolume(1);

        //uniform distribution
        PSD psd;
        psd.setParticleSizeDistribution({{0, 0},{0.015, 1}});
        //PSD::convertCumulativeVolumeToNumber(psd); //PSND
        insertionBoundary.setPSD(psd);
        
        //add the insertion boundary to the handler
        auto i = boundaryHandler.copyAndAddObject(insertionBoundary);
        i->checkBoundaryBeforeTimeStep(this);
    
        logger(INFO, "N=%", particleHandler.getSize());
    }
};

int main(int argc UNUSED, char* argv[] UNUSED)
{
    PSDSelfTest problem;
    problem.setupInitialConditions();
    problem.writeDataFile();
    //system("python UniformPSDSelfTest.py");
}
