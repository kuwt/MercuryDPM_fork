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
#include "Species/LinearViscoelasticSpecies.h"
#include "Boundaries/FluxBoundary.h"

/*! 
 * \brief Tests if FluxBoundary correctly measures the flux of particles crossing the boundary.
 * \details The problem consists of three particle moving at constant velocity and one flux boundary.
 */
class FluxBoundaryUnitTest : public Mercury3D
{
public:

    FluxBoundaryUnitTest()
    {
        setName("FluxBoundaryUnitTest");

        LinearViscoelasticSpecies species;
        species.setDensity(1.0);
        species.setCollisionTimeAndRestitutionCoefficient(.01, .1, 1.);

        speciesHandler.copyAndAddObject(species);
        speciesHandler.copyAndAddObject(species);
    }

    void setupInitialConditions() override {
        //define the domain
        setMax(Vec3D(1, 0.5, 0.5));
        setMin(Vec3D(-1, -0.5, -0.5));

        //define moving particles
        SphericalParticle p0;

        p0.setSpecies(speciesHandler.getObject(0));
        p0.setRadius(0.1);

        p0.setPosition(Vec3D(-0.1, 0.0, 0.0));
        p0.setVelocity(Vec3D(0.2, 0.0, 0.0));
        particleHandler.copyAndAddObject(p0);

        p0.setPosition(Vec3D(-0.1, 0.3, 0.0));
        p0.setVelocity(Vec3D(0.2, 0.0, 0.0));
        particleHandler.copyAndAddObject(p0);

        p0.setPosition(Vec3D(0.1, -0.3, 0.0));
        p0.setVelocity(Vec3D(-0.2, 0.0, 0.0));
        particleHandler.copyAndAddObject(p0);

        // Define a flux boundary
        FluxBoundary fluxBoundary;
        fluxBoundary.set(Vec3D(1, 0, 0), 0.0);

        boundaryHandler.copyAndAddObject(fluxBoundary);
    }
};

int main(int argc, char *argv[])
{
    //declare the DPM problem and set the name
    FluxBoundaryUnitTest problem;
    problem.setTimeMax(1);
    problem.setTimeStep(0.1);
    problem.setFileType(FileType::NO_FILE);
    
    //run the simulation
    problem.solve(argc, argv);

    // Check the value of the fluxBoundary
    FluxBoundary* fluxBoundary = static_cast<FluxBoundary*>(problem.boundaryHandler.getObject(0));

    logger.assert_always(fluxBoundary->getNumberOfParticlesCrossedNet() == 1, "Wrong number of particles crossed the flux boundary");
    logger.assert_always(fluxBoundary->getNumberOfParticlesCrossedForw() == 2, "Wrong number of particles crossed the flux boundary (Forward)");
    logger.assert_always(fluxBoundary->getNumberOfParticlesCrossedForw() == 2, "Wrong number of particles crossed the flux boundary (Backward)");

    // Get and ckeck volume of particles
    Mdouble particleVolume = problem.particleHandler.getObject(0)->getVolume();
    logger.assert_always(fluxBoundary->getNumberOfParticlesCrossedNet()*particleVolume == fluxBoundary->getVolumeOfParticlesCrossedNet(), "Wrong volume of particles crossed the flux boundary");
    logger.assert_always(fluxBoundary->getNumberOfParticlesCrossedForw()*particleVolume == fluxBoundary->getVolumeOfParticlesCrossedForw(), "Wrong volume of particles crossed the flux boundary (Forward)");
    logger.assert_always(fluxBoundary->getNumberOfParticlesCrossedBack()*particleVolume == fluxBoundary->getVolumeOfParticlesCrossedBack(), "Wrong volume of particles crossed the flux boundary (Backward)");

    // Get and check mass of particles
    Mdouble particleMass = problem.particleHandler.getObject(0)->getMass();
    logger.assert_always(fluxBoundary->getNumberOfParticlesCrossedNet()*particleMass == fluxBoundary->getMassOfParticlesCrossedNet(), "Wrong mass of particles crossed the flux boundary");
    logger.assert_always(fluxBoundary->getNumberOfParticlesCrossedForw()*particleMass == fluxBoundary->getMassOfParticlesCrossedForw(), "Wrong mass of particles crossed the flux boundary (Forward)");
    logger.assert_always(fluxBoundary->getNumberOfParticlesCrossedBack()*particleMass == fluxBoundary->getMassOfParticlesCrossedBack(), "Wrong mass of particles crossed the flux boundary (Backward)");
}

