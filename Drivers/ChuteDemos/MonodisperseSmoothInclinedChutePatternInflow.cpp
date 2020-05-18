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
#include <Particles/BaseParticle.h>

/**
 * This code simulates monodisperse flow over smooth inclined channel.
 * \todo remove/add commented out code before the release
*/
class SmoothChute : public Chute
{
public:

    /*!
     * In the time loop, before every time step, check if we are still in the
     * period that particles should be inserted. If so, insert 11 x 26
     * particles at the various positions in the chute.
     */
    void actionsBeforeTimeStep() override {
        Chute::cleanChute();
        if (getTime() < insertTime && getTime() + getTimeStep() > insertTime)
        {
            ///\todo make a more sensible logger message
            logger(VERBOSE, "hello % : %", getTime(), getTime() + getTimeStep());
            SphericalParticle p0;
            p0.setSpecies(speciesHandler.getObject(0));
            p0.setRadius(0.0015);
            Vec3D pos;
            for (unsigned int i = 1; i < 12; i++)
            {
                for (unsigned int j = 1; j < 27; j++)
                {
                    pos.X = p0.getRadius();
                    pos.Y = p0.getRadius() * (2.0 * j - 1.0);
                    pos.Z = p0.getRadius() * (2.0 * i - 1.0);
                    p0.setPosition(pos);
                    p0.setVelocity(Vec3D(0.5, 0.0, 0.0));
                    particleHandler.copyAndAddObject(p0);
                }
            }
            logger(VERBOSE, "Number of particles before time step: %",
                            particleHandler.getNumberOfObjects());
            /*for (int i=particleHandler.getNumberOfObjects()-numberOfParticlesToBeInserted; i<particleHandler.getNumberOfObjects(); i++)
            {
                    particleHandler.getObject(i)->setForce(Vec3D(0,0,0));
                    particleHandler.getObject(i)->setTorque(Vec3D(0,0,0));
            }*/
            ///\todo IFCD: what happens here?
            insertTime = insertTime + insertTimeInterval;
        }
    }

    /*!
     * Setup the initial conditions of the chute: set the boundaries of the
     * computational domain, setup the side walls and add some initial
     * particles.
     */
    void setupInitialConditions() override {
        //Rectangular Chute geometry
        setXMax(1.0);
        setXMin(0.0);
        setYMax(0.08);
        setYMin(0.0);
        setZMax(20. * 0.003);
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
        for (unsigned int i = 1; i < 12; i++)
        {
            for (unsigned int j = 1; j < 27; j++)
            {
                pos.X = p0.getRadius();
                pos.Y = p0.getRadius() * (2.0 * j - 1.0);
                pos.Z = p0.getRadius() * (2.0 * i - 1.0);
                p0.setPosition(pos);
                p0.setVelocity(Vec3D(0.5, 0.0, 0.0));
                particleHandler.copyAndAddObject(p0);
            }
        }

    }

    ///\brief the amount of time that new particles should be inserted.
    double insertTime;
    ///\todo what does this do?
    double insertTimeInterval;
};

int main(int argc UNUSED, char *argv[] UNUSED)
{
    SmoothChute problem;
    problem.setName("FlowOverSmoothInclinedChannelPatternInflow");
    problem.setTimeMax(0.001);
    problem.setChuteAngleAndMagnitudeOfGravity(30.0,9.81);
    
    problem.insertTime = 0.01;
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

    particleSpecies->setDensity(2550.0); // sets the species type-0 density
    particleSpecies->setStiffness(1000.0);// sets the spring stiffness
    particleSpecies->setSlidingStiffness(2.0/7.0*particleSpecies->getStiffness());
    particleSpecies->setDissipation(0.0);// sets the dissipation
    particleSpecies->setSlidingFrictionCoefficient(0.12);

    wallParticleSpecies->setStiffness(979.38);
    wallParticleSpecies->setDissipation(particleSpecies->getDissipation());
    wallParticleSpecies->setSlidingStiffness(2.0/7.0*wallParticleSpecies->getStiffness());
    wallParticleSpecies->setSlidingFrictionCoefficient(0.22);

    //For this problem, we want a datafile and a restartfile as output, and
    // they should be printed every 100 time steps.
    problem.setSaveCount(100);
    problem.dataFile.setFileType(FileType::ONE_FILE);
    problem.restartFile.setFileType(FileType::ONE_FILE);
    problem.fStatFile.setFileType(FileType::NO_FILE);
    problem.eneFile.setFileType(FileType::NO_FILE);

    //problem.autoNumber();
    problem.setXBallsAdditionalArguments("-solidf -v0");

    problem.setTimeStep(1e-4);
    problem.solve(argc, argv);
    
}
