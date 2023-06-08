//Copyright (c) 2013-2023, The MercuryDPM Developers Team. All rights reserved.
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

/*!
 * \note The file contains the function definitions, so no additional linking is required. 
 */

#ifndef CONTACT_LAW_TEST_HELPERS
#define CONTACT_LAW_TEST_HELPERS

#include "DPMBase.h"
#include "Helpers/FileIOHelpers.h"
#include "Logger.h"
#include "Particles/BaseParticle.h"
#include "Species/BaseSpecies.h"
#include "Walls/InfiniteWall.h"

/**
 * Creates a DPMBase with a particles of unit size and a flat wall and loads/unloads the particle-wall contact
 * \param[in] species      particle species specifying the contact law
 * \param[in] displacement peak displacement before unloading
 * \param[in] velocity     loading/unloading velocity
 */
void loadingTest(const ParticleSpecies* species, Mdouble displacement, Mdouble velocity, Mdouble radius,
                          std::string name)
{
    class LoadingTest : public DPMBase
    {
        const ParticleSpecies* species;
        Mdouble displacement;
        Mdouble velocity;
        Mdouble radius;
    public:
        //public variables
        LoadingTest(const ParticleSpecies* species, Mdouble displacement, Mdouble velocity, Mdouble radius)
                : species(species), displacement(displacement), velocity(velocity), radius(radius)
        {}

        void setupInitialConditions() override
        {
            //setName("LoadingTest"+species->getName());
            setTimeMax(2.0 * displacement / velocity);
            setTimeStep(2e-3 * getTimeMax());
            setSaveCount(1);
            setFileType(FileType::NO_FILE);
            fStatFile.setFileType(FileType::ONE_FILE);

            setMax({radius, radius, radius + radius});
            setMin({-radius, -radius, 0});
            setSystemDimensions(3);
            setParticleDimensions(3);

            speciesHandler.copyAndAddObject(*species);

            SphericalParticle p;
            p.setSpecies(speciesHandler.getObject(0));
            p.setRadius(radius);
            p.setPosition({0, 0, radius});
            particleHandler.copyAndAddObject(p);

            InfiniteWall w;
            w.setSpecies(speciesHandler.getObject(0));
            w.set(Vec3D(0, 0, -1), Vec3D(0.0, 0.0, 0.0));
            wallHandler.copyAndAddObject(w);
        }

        void actionsBeforeTimeStep() override
        {
            BaseParticle* p = particleHandler.getLastObject();
            logger.assert_debug(p,"Empty particle handler");
            p->setAngularVelocity({0, 0, 0});

            //Moving particle normally into surface
            if (getTime() <= displacement / velocity)
            {
                p->setVelocity({0, 0, velocity});
                p->setPosition({0, 0, radius - velocity * getTime()});
            }
            else
            {
                p->setVelocity({0, 0, -velocity});
                p->setPosition({0, 0, radius - displacement - displacement + velocity * getTime()});
            }
        }
    } test(species, displacement, velocity, radius);
    test.setName(name);
    test.solve();
    helpers::writeToFile(test.getName() + ".gnu", "plot '" + test.getName() + ".fstat' u 7:9 w lp");
    logger(INFO, "finished loading test: run 'gnuplot %.gnu' to view output", test.getName());
}

/**
 * Creates a DPMBase with a particles of unit size and a flat wall and loads/unloads/reloads the particle-wall contact in tangential direction
 * \param[in] species      particle species specifying the contact law
 * \param[in] displacement peak displacement before unloading
 * \param[in] velocity     loading/unloading velocity
 */
void normalAndTangentialLoadingTest(const ParticleSpecies* species, Mdouble displacement,
                                             Mdouble tangentialDisplacement, Mdouble velocity, Mdouble radius,
                                             std::string name)
{
    class LoadingTest : public DPMBase
    {
        const ParticleSpecies* species;
        Mdouble displacement;
        Mdouble tangentialDisplacement;
        Mdouble velocity;
        Mdouble radius;
    public:
        //public variables
        LoadingTest(const ParticleSpecies* species, Mdouble displacement, Mdouble tangentialDisplacement,
                    Mdouble velocity, Mdouble radius)
                : species(species), displacement(displacement), tangentialDisplacement(tangentialDisplacement),
                  velocity(velocity), radius(radius)
        {}

        void setupInitialConditions() override
        {
            //setName("TangentialLoadingTest"+species->getName());
            setTimeMax(4.0 * tangentialDisplacement / velocity);
            setTimeStep(4e-4 * getTimeMax());
            setSaveCount(1);
            setFileType(FileType::NO_FILE);
            fStatFile.setFileType(FileType::ONE_FILE);

            setMax({radius, radius, radius + radius});
            setMin({-radius, -radius, 0});
            setSystemDimensions(3);
            setParticleDimensions(3);

            speciesHandler.copyAndAddObject(*species);

            SphericalParticle p;
            p.setSpecies(speciesHandler.getObject(0));
            p.setRadius(radius);
            p.setPosition({0, 0, radius - displacement});
            particleHandler.copyAndAddObject(p);

            InfiniteWall w;
            w.setSpecies(speciesHandler.getObject(0));
            w.set(Vec3D(0, 0, -1), Vec3D(0.0, 0.0, 0.0));
            wallHandler.copyAndAddObject(w);
        }

        void actionsBeforeTimeStep() override
        {
            BaseParticle* p = particleHandler.getLastObject();
            logger.assert_debug(p,"Empty particle handler");
            p->setAngularVelocity({0, 0, 0});
    
            //Moving particle cyclically right and left between +-tangentialDisplacement
            bool moveRight = static_cast<int>(getTime() / (2.0*tangentialDisplacement / velocity) +0.5)%2==0;
            if (moveRight)
            {
                p->setVelocity({-velocity, 0, 0});
                p->setPosition({tangentialDisplacement - velocity * getTime(), 0, radius - displacement});
            }
            else
            {
                p->setVelocity({velocity, 0, 0});
                p->setPosition({-2*tangentialDisplacement + velocity * getTime(), 0, radius - displacement});
            }
        }

    } test(species, displacement, tangentialDisplacement, velocity, radius);
    test.setName(name);
    test.solve();
    helpers::writeToFile(test.getName() + ".gnu", "plot '" + test.getName() + ".fstat' u 8:($10*$14) w lp");
    logger(INFO, "finished tangential loading test: run 'gnuplot %.gnu' to view output", test.getName());
}

/**
 * Creates a DPMBase with a particles of unit size and a flat wall, loads the particle-wall contact in normal and tangential direction, then rotates.
 * \param[in] species      particle species specifying the contact law
 * \param[in] displacement peak displacement before unloading
 * \param[in] velocity     loading/unloading velocity
 */
void objectivenessTest(const ParticleSpecies* species, Mdouble displacement, Mdouble tangentialDisplacement,
                                Mdouble velocity, Mdouble radius, std::string name)
{
    class ObjectivenessTest : public DPMBase
    {
        const ParticleSpecies* species;
        Mdouble displacement;
        Mdouble tangentialDisplacement;
        Mdouble velocity;
        Mdouble radius;
    public:
        //public variables
        ObjectivenessTest(const ParticleSpecies* species, Mdouble displacement, Mdouble tangentialDisplacement,
                          Mdouble velocity, Mdouble radius)
                : species(species), displacement(displacement), tangentialDisplacement(tangentialDisplacement),
                  velocity(velocity), radius(radius)
        {}

        void setupInitialConditions() override
        {
            //setName("ObjectivenessTest"+species->getName());
            setTimeMax((tangentialDisplacement + 0.5 * constants::pi * radius) / velocity);
            setTimeStep(1e-4 * getTimeMax());
            setSaveCount(20);
            setFileType(FileType::NO_FILE);
            dataFile.setFileType(FileType::ONE_FILE);
            fStatFile.setFileType(FileType::ONE_FILE);

            setMax(radius * Vec3D(2, 2, 2));
            setMin(radius * Vec3D(-2, -2, -2));
            setSystemDimensions(2);
            setParticleDimensions(3);

            speciesHandler.copyAndAddObject(*species);
    
            SphericalParticle p;
            p.setSpecies(speciesHandler.getObject(0));
            p.setRadius(radius);
            p.setPosition({0, radius - displacement, 0});
            particleHandler.copyAndAddObject(p);
            p.setPosition({0, -radius + displacement, 0});
            particleHandler.copyAndAddObject(p);
        }

        void actionsBeforeTimeStep() override
        {
            BaseParticle* p = particleHandler.getObject(0);
            BaseParticle* q = particleHandler.getLastObject();
            logger.assert_debug(p,"Empty particle handler");
            logger.assert_debug(q,"Empty particle handler");

            //Moving particle normally into surface
            if (getTime() <= tangentialDisplacement / velocity)
            {
                p->setAngularVelocity({0, 0, 0});
                p->setVelocity({velocity, 0, 0});
                p->setPosition({-tangentialDisplacement + velocity * getTime(), radius - displacement, 0});
                q->setAngularVelocity({0, 0, 0});
                q->setVelocity({-velocity, 0, 0});
                q->setPosition({tangentialDisplacement - velocity * getTime(), -radius + displacement, 0});
            }
            else
            {
                Mdouble angle = velocity / (radius - displacement) * (getTime() - tangentialDisplacement / velocity);
                Mdouble s = sin(angle);
                Mdouble c = cos(angle);
                p->setAngularVelocity(velocity / (radius - displacement) * Vec3D(0, 0, -1));
                //p->setAngularVelocity(Vec3D(0,0,0));
                p->setOrientation({1, 0, 0, 0});
                p->setVelocity(velocity * Vec3D(c, -s, 0));
                //p->setVelocity(Vec3D(0,0,0));
                p->setPosition((radius - displacement) * Vec3D(s, c, 0));
                q->setAngularVelocity(-p->getAngularVelocity());
                q->setOrientation(-p->getOrientation());
                q->setVelocity(-p->getVelocity());
                q->setPosition(-p->getPosition());
            }
        }

    } test(species, displacement, tangentialDisplacement, velocity, radius);
    test.setName(name);
    test.solve();
    helpers::writeToFile(test.getName() + ".gnu", "set size ratio -1; plot '" + test.getName() + ".fstat' u 14:15 every 2 w lp");
    logger(INFO, "finished objectiveness test: run 'gnuplot %.gnu' to view output", test.getName());
}

#endif // CONTACT_LAW_TEST_HELPERS
