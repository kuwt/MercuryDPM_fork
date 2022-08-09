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

#include "Chute.h"
#include "Species/LinearViscoelasticFrictionSpecies.h"
#include "Boundaries/PeriodicBoundary.h"

class ChutePeriodic : public Chute
{
public:
    ChutePeriodic()
    {
        setName("ChutePeriodic");

        //Set and add the particle-species.
        LinearViscoelasticFrictionSpecies species;
        species.setDensity(6 / constants::pi);
        species.setCollisionTimeAndRestitutionCoefficient(5e-3, 0.88, 1);
        species.setSlidingDissipation(species.getDissipation());
        species.setSlidingStiffness(2. / 7. * species.getStiffness());
        species.setSlidingFrictionCoefficient(0.5);
        speciesHandler.copyAndAddObject(species);
        
        // Chute properties
        setChuteAngleAndMagnitudeOfGravity(24,9.81);
        setChuteLength(30);
        setChuteWidth(10);
        setInflowHeight(10);
        setInflowParticleRadius(0.5);
        setFixedParticleRadius(getInflowParticleRadius());
        setRoughBottomType(MULTILAYER);

        //make chute periodic in y-direction:
        makeChutePeriodic();

        //make the chute periodic in x-direction:
        PeriodicBoundary b0;
        b0.set(Vec3D(1.0, 0.0, 0.0), getXMin(), getXMax());
        boundaryHandler.copyAndAddObject(b0);
        
        //Set time-related properties
        setTimeMax(10);
        setTimeStep(1.0e-4); //(1/50th of the collision time)
        setSaveCount(10000);
    }

    ///add particles
    void setupInitialConditions() override
    {
        Chute::setupInitialConditions();
        addFlowParticles();
        // todo: understand VTK and Xball output
//        setParticlesWriteVTK(true);
    }
    
    ///Add initial flow particles.
    ///\todo replace addFlowParticles and createInflowParticle with insertion boundary
    void addFlowParticles()
    {
        logger(INFO, "Adding flowing particles");
        // *Real* is added below by Marnix to avoid some warnings when running in parallel
        unsigned int N = particleHandler.getNumberOfRealObjects() + getChuteLength() * getChuteWidth() * getInflowHeight() /
                                                                mathsFunc::cubic(getInflowParticleRadius()) / 8;
        particleHandler.setStorageCapacity((N));
        setZMax(1.2 * getInflowHeight());
        while (particleHandler.getSize() < N)
        {
            SphericalParticle p = createInflowParticle();
            if (checkParticleForInteraction(p))
            {
                particleHandler.copyAndAddObject(p);
            }
            else
            {
                setInflowHeight(getInflowHeight() + .0001 * getMaxInflowParticleRadius());
            }
        }
        logger(INFO, "InflowHeight= %", getInflowHeight());
        setHGridDistribution(HGridDistribution::EXPONENTIAL);
    }
    
    ///creates flow particles in the whole chute
    SphericalParticle createInflowParticle()
    {
        SphericalParticle p0;
        p0.setSpecies(speciesHandler.getObject(0));
        p0.setRadius(random.getRandomNumber(getMinInflowParticleRadius(), getMaxInflowParticleRadius()));
        p0.setPosition(Vec3D(random.getRandomNumber(getXMin() + 2.0 * p0.getRadius(), getXMax()),
                             random.getRandomNumber(getYMin() + 2.0 * p0.getRadius(), getYMax()),
                             random.getRandomNumber(getZMin() + 2.0 * p0.getRadius(), getInflowHeight())));
        p0.setVelocity(Vec3D(0.0, 0.0, 0.0));
        return p0;
    }
};

