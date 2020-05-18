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

#include "Chute.h"
#include <Species/LinearViscoelasticFrictionSpecies.h>
#include <Walls/InfiniteWall.h>

/*!
 * This code simulates monodisperse flow over smooth inclined channel.
 * \todo decide what to do with commented out code
 * \bug this code gives an HGrid warning in the beginning
 */

class SmoothChute : public Chute
{
public:

    void actionsBeforeTimeStep() override {
        Chute::cleanChute();
        if (getTime() < insertTime && getTime() + getTimeStep() > insertTime)
        {
            ///\todo this message was commented out in a cout, maybe make a more sensible message
            logger(VERBOSE, "hello % : %", getTime(), getTime() + getTimeStep());
            SphericalParticle p0;
            p0.setSpecies(speciesHandler.getObject(0));
            p0.setRadius(0.0015);
            Vec3D pos;
            unsigned int numberOfParticlesInserted = 0;
            while (numberOfParticlesInserted < numberOfParticlesToBeInserted)
            {
                pos.X = p0.getRadius();
                pos.Y = random.getRandomNumber(getYMin(), getYMax());
                pos.Z = random.getRandomNumber(getZMin(), getZMax());
                p0.setPosition(pos);
                p0.setVelocity(Vec3D(0.5, 0.0, 0.0));
                if (checkParticleForInteraction(p0))
                {
                    particleHandler.copyAndAddObject(p0);
                    ++numberOfParticlesInserted;
                }
            }
            logger(VERBOSE, "Number of particles inserted at time %: %", insertTime, particleHandler.getNumberOfObjects());
            for (unsigned int i = particleHandler.getNumberOfObjects() - numberOfParticlesToBeInserted;
                 i < particleHandler.getNumberOfObjects(); ++i)
            {
                particleHandler.getObject(i)->setForce(Vec3D(0, 0, 0));
                particleHandler.getObject(i)->setTorque(Vec3D(0, 0, 0));
            }
            insertTime = insertTime + insertTimeInterval;
        }
    }

    void setupInitialConditions() override {
        //Rectangular Chute geometry
        setXMax(300 * 0.006);
        setXMin(0.0);
        setYMax(0.08);
        setYMin(0.0);
        setZMax(10. * 0.006);
        setZMin(0.0);

        InfiniteWall w0;

        w0.set(Vec3D(0.0, 0.0, -1.0), Vec3D(0.0, 0.0, getZMin()));
        w0.setSpecies(speciesHandler.getObject(1));
        wallHandler.copyAndAddObject(w0);

        w0.set(Vec3D(0.0, 1.0, 0.0), Vec3D(0.0, getYMax(), 0.0));
        w0.setSpecies(speciesHandler.getObject(1));
        wallHandler.copyAndAddObject(w0);

        w0.set(Vec3D(0.0, -1.0, 0.0), Vec3D(0.0, getYMin(), 0.0));
        w0.setSpecies(speciesHandler.getObject(1));
        wallHandler.copyAndAddObject(w0);

        // Initial particle position and velocity
        SphericalParticle p0;
        Vec3D pos;
        p0.setSpecies(speciesHandler.getObject(0));
        p0.setRadius(0.0015);
        int numberOfParticlesInserted = 0;
        while (numberOfParticlesInserted < numberOfParticlesToBeInserted)
        {
            pos.X = p0.getRadius();
            pos.Y = random.getRandomNumber(getYMin(), getYMax());
            pos.Z = random.getRandomNumber(getZMin(), getZMax());
            p0.setPosition(pos);
            p0.setVelocity(Vec3D(0.5, 0.0, 0.0));
            //
            if (checkParticleForInteraction(p0))
            {
                particleHandler.copyAndAddObject(p0);
                ++numberOfParticlesInserted;
            }
        }
        logger(INFO, "Particles inserted while setting up initial conditions: %", numberOfParticlesInserted);
    }

    double insertTime;
    double insertTimeInterval;
    unsigned int numberOfParticlesToBeInserted;
};

int main(int argc, char *argv[])
{
    logger(INFO, "Chute flow over a smooth bottom");

    //set the problem parameters: name, simulation time and time step
    SmoothChute problem;
    problem.setName("MonodisperseSmoothInclinedChute");
    problem.setTimeMax(100);
    problem.setTimeStep(1e-4);

    //set the direction of the gravity as 30 degrees, and magnitude 9.81 m/s^2
    problem.setChuteAngleAndMagnitudeOfGravity(30.0, 9.81);
    
    problem.numberOfParticlesToBeInserted = 286;
    problem.insertTime = 0.006;
    problem.insertTimeInterval = problem.insertTime;
    //
    auto particleSpecies = problem.speciesHandler.copyAndAddObject(LinearViscoelasticFrictionSpecies());
    auto wallSpecies = problem.speciesHandler.copyAndAddObject(LinearViscoelasticFrictionSpecies());
    auto wallParticleSpecies = problem.speciesHandler.getMixedObject(wallSpecies, particleSpecies);
    
    // Mdouble effectiveMass = 0.5*particleSpecies->getMassFromRadius(0.5*0.003);
    // std::cout << "Mass" << 2.0*effectiveMass << std::endl;
    // helpers::KAndDisp kAndDispPP = helpers::computeKAndDispFromCollisionTimeAndRestitutionCoefficientAndEffectiveMass(0.005,0.97,effectiveMass);
    //helpers::KAndDisp kAndDispPW = helpers::computeKAndDispFromCollisionTimeAndRestitutionCoefficientAndEffectiveMass(0.005,0.95,effectiveMass);
    
    //std::cout << "Stifness" << kAndDispPP.k << std::endl;

    //set parameters for the particle species: density, stiffness, dissipation and sliding friction coefficient
    particleSpecies->setDensity(2550.0);
    particleSpecies->setStiffness(1000.0);
    particleSpecies->setSlidingStiffness(2.0 / 7.0 * particleSpecies->getStiffness());
    particleSpecies->setDissipation(0.0);
    particleSpecies->setSlidingFrictionCoefficient(0.12);

    //set the parameters for the wall species
    wallParticleSpecies->setStiffness(979.38);
    wallParticleSpecies->setDissipation(particleSpecies->getDissipation());
    wallParticleSpecies->setSlidingStiffness(2. / 7.0 * wallParticleSpecies->getStiffness());
    wallParticleSpecies->setSlidingFrictionCoefficient(0.22);

    //set the output parameters: how many time steps should be skipped before a next time step is written to the files,
    //and which files should be written.
    problem.setSaveCount(100);
    problem.dataFile.setFileType(FileType::ONE_FILE);
    problem.restartFile.setFileType(FileType::ONE_FILE);
    problem.fStatFile.setFileType(FileType::NO_FILE);
    problem.eneFile.setFileType(FileType::NO_FILE);
    //
    //problem.autoNumber();
    problem.setXBallsAdditionalArguments("-solidf -v0");

    //run the simulation with all parameters given before
    problem.solve(argc, argv);
    
    return 0;
}
