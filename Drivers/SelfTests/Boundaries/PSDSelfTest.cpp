//Copyright (c) 2013-2018, The MercuryDPM Developers Team. All rights reserved.
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
#include "Boundaries/PeriodicBoundary.h"
#include "Particles/BaseParticle.h"
#include "Species/LinearViscoelasticSpecies.h"
#include "Walls/InfiniteWall.h"


class PSDSelfTest : public Mercury3D
{
public:
    
    void setupInitialConditions() override
    {
        setName("PSDSelfTest");
        setSystemDimensions(3);
        setGravity(Vec3D(0, 0, 0));
        setTimeStep(1e-4);
        dataFile.setSaveCount(10);
        setTimeMax(2e-2);
        setHGridMaxLevels(2);
        
        setMin(Vec3D(0, 0, 0));
        setMax(Vec3D(.1, .1, .1));
        
        LinearViscoelasticSpecies species;
        species.setDensity(2000);
        species.setStiffness(10000);
        speciesHandler.copyAndAddObject(species);
        
        
        BaseParticle insertionBoundaryParticle;
        insertionBoundaryParticle.setSpecies(speciesHandler.getObject(0));
        
        CubeInsertionBoundary insertionBoundary;
        insertionBoundary.set(&insertionBoundaryParticle, 1, getMin(), getMax(),
                Vec3D(0, 0, 0), Vec3D(0, 0, 0), 1, 2);
        const std::vector<PSD> APAPM = {{0.000387645, 0},
                                        {0.000405195, 0.00641963},
                                        {0.000423541, 0.0157881},
                                        {0.000442717, 0.0289129},
                                        {0.000462761, 0.0475877},
                                        {0.000483713, 0.0702285},
                                        {0.000505613, 0.096661},
                                        {0.000528505, 0.12752},
                                        {0.000552434, 0.162984},
                                        {0.000577445, 0.199952},
                                        {0.00060359,  0.236205},
                                        {0.000630917, 0.270216},
                                        {0.000659483, 0.300988},
                                        {0.000689341, 0.328802},
                                        {0.000720551, 0.354678},
                                        {0.000753175, 0.379335},
                                        {0.000787275, 0.404133},
                                        {0.000822919, 0.430189},
                                        {0.000860177, 0.458372},
                                        {0.000899122, 0.489708},
                                        {0.000939831, 0.524177},
                                        {0.000982382, 0.561866},
                                        {0.00102686,  0.602098},
                                        {0.00107335,  0.643887},
                                        {0.00112195,  0.686123},
                                        {0.00117275,  0.727694},
                                        {0.00122584,  0.767495},
                                        {0.00128134,  0.804712},
                                        {0.00133936,  0.838663},
                                        {0.0014,      0.868961},
                                        {0.00146338,  0.895445},
                                        {0.00152964,  0.918157},
                                        {0.00159889,  0.937242},
                                        {0.00167128,  0.952947},
                                        {0.00174695,  0.965604},
                                        {0.00182605,  0.975541},
                                        {0.00190872,  0.983136},
                                        {0.00199514,  0.988764},
                                        {0.00208547,  0.992798},
                                        {0.00217989,  0.995601},
                                        {0.00227859,  0.997454},
                                        {0.00238175,  0.998635},
                                        {0.00248959,  0.999345},
                                        {0.0026023,   0.999741},
                                        {0.00272012,  0.999939},
                                        {0.00284328,  0.999995},
                                        {0.00297201,  1},};
        insertionBoundary.setPSD(APAPM);
        boundaryHandler.copyAndAddObject(insertionBoundary);
        
    }
    
    void printTime() const override
    {
        logger(INFO, "t=%, tMax=%, N=%", getTime(), getTimeMax(), particleHandler.getSize());
    }
};

int main(int argc UNUSED, char* argv[] UNUSED)
{
    logger(INFO, "Simple box for creating particles");
    
    PSDSelfTest insertionBoundary_problem;
    insertionBoundary_problem.solve();
}
