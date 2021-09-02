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



#include "Mercury2D.h"
#include "Walls/InfiniteWall.h"
#include <Species/HertzianViscoelasticMindlinSpecies.h>

/// This case does a single elastic particle performing an oblique impact on an infinite plane.
/// The k is chosen so that the maximum overlap with the wall is around 2% of the particles diameter;
/// whereas, the time step must be taken to ensure 50 steps with a collision.
class ObliqueImpactSelfTest : public Mercury2D
{
public:
    
    void setupInitialConditions() override
    {
        Mdouble radius = 0.005;
        Mdouble dist = radius + 0.000005;
        
        setMax({0.01,0.01,0.0});
        setMin(0.0, 0.0, 0.0);
//        setGravity({0.0,-9.8,0.0});
        setTimeMax(0.0005);
        
        InfiniteWall w0;
        w0.setSpecies(speciesHandler.getObject(0));
        w0.set(Vec3D(0,-1,0), Vec3D(0, getYMin(), 0));
        wallHandler.copyAndAddObject(w0);
        
        SphericalParticle p0;
        p0.setSpecies(speciesHandler.getObject(0));
        p0.setPosition(Vec3D(getXMin() + radius, dist,0.0));
        p0.setVelocity(Vec3D(10.0,-1.0,0.0));
        p0.setRadius(0.005);
        particleHandler.copyAndAddObject(p0);
    }
    
};

int main(int argc, char* argv[])
{
    logger(INFO, "Single particle hitting the bottom plate in an oblique impact");
    // Make the problem and set the name
    ObliqueImpactSelfTest ObliqueImpactSelfTestProblem;
    ObliqueImpactSelfTestProblem.setName("ObliqueImpactSelfTest");
    
    //Set the species of the particle and wall, and its properties
    HertzianViscoelasticMindlinSpecies species;
    species.setDensity(6./constants::pi);
    ObliqueImpactSelfTestProblem.setParticleDimensions(3);
    species.setEffectiveElasticModulusAndPoissonRatio(1e5, 0.3);
//    species.setDissipation(0.2);
//    species.setSlidingDissipation(0.9);
    species.setSlidingFrictionCoefficient(0.1);
    ObliqueImpactSelfTestProblem.speciesHandler.copyAndAddObject(species);
    
    //set the parameters for the solver
    ObliqueImpactSelfTestProblem.setSaveCount(500);
    ObliqueImpactSelfTestProblem.setFileType(FileType::ONE_FILE);
    ObliqueImpactSelfTestProblem.fStatFile.setFileType(FileType::ONE_FILE);
//    ObliqueImpactSelfTestProblem.setWallsWriteVTK(FileType::NO_FILE);
    ObliqueImpactSelfTestProblem.setTimeStep(0.25e-8);
    
    //solve the system, the single particle will now bounce on the plate
    ObliqueImpactSelfTestProblem.solve(argc, argv);
    helpers::writeToFile(ObliqueImpactSelfTestProblem.getName() + ".gnu", "plot '" + ObliqueImpactSelfTestProblem.getName() + ".fstat' u 8:($10*$14) w lp");
    logger(INFO, "finished oblique impact test: run 'gnuplot %.gnu' to view output", ObliqueImpactSelfTestProblem.getName());
}

