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

#include <Species/LinearViscoelasticSlidingFrictionSpecies.h>
#include "BidispersedChute.h"

BidispersedChute::BidispersedChute(const BidispersedChuteParameters& bidispersedChuteParameters) :
        parameters(bidispersedChuteParameters)
{
    setName("BidispersedChute");
    setSpeciesProperties();
    setChuteProperties();
    setTimeAndSaveCount();
    distribution = std::normal_distribution<Mdouble>(1.0, 0.025);
}

/*!
 * Overrides DPMBase::setupInitialConditions to setup the initial configuration of the system.
 * Calls subfunctions that define the content of the wall, boundary and particle handlers.
 * Finally, sets the maser to be closed initially.
 */
void BidispersedChute::setupInitialConditions()
{
    //hack to be able to override construction of the species in derived classes
    speciesHandler.clear();
    this->setSpeciesProperties();
    //make sure to create the bottom of the chute after all species etc. are set.
    createBottom();
    //Set the species for the bottom particles. 1 for large particles, 0 for particles with d=1.
    for (BaseParticle* const p : particleHandler)
    {
        p->setSpecies(speciesHandler.getObject(1));
    }
    //Make sure that there are no gaps that particles can fall through
    InfiniteWall w0;
    w0.setSpecies(speciesHandler.getObject(0));
    w0.set(Vec3D(0.0, 0.0, -1.0), Vec3D(0, 0, -1));
    wallHandler.copyAndAddObject(w0);
    
    //after the bottom is made, the boundaries can be put in place and the particles inserted
    setBoundaries();
    
    unsigned int numberOfLargeToGo = parameters.getNumberOfLargeParticles(getChuteLength() * getChuteWidth());
    unsigned int numberOfSmallToGo = parameters.getNumberOfSmallParticles(getChuteLength() * getChuteWidth());
    insertParticles(numberOfLargeToGo, numberOfSmallToGo);
}

void BidispersedChute::setStandardDeviation(Mdouble sigma)
{
    distribution = std::normal_distribution<Mdouble>(1.0, sigma);
}

/*!
 * A subfunction of setupInitialConditions; determines the initial content of the particleHandler.
 * It inserts the particles randomly into the domain in order to achieve a homogeneous distribution.
 * Note, large particles are assigned to species 1, small particles are assigned species 0,
 * fixed particles are species 0 (by default).
 */
void BidispersedChute::insertParticles(unsigned int numberOfLargeToGo, unsigned int numberOfSmallToGo)
{
    logger(DEBUG, "need % small and % large particles", numberOfLargeToGo, numberOfSmallToGo);
    
    while ((numberOfLargeToGo > 0 || numberOfSmallToGo > 0))
    {
        if (numberOfSmallToGo > 0)
        {
            insertSmallParticle(numberOfSmallToGo);
        }
        else
        {
            insertLargeParticle(numberOfLargeToGo);
        }
    }
}

/*!
 * A subfunction of insertParticles; inserts a large particle (of species 1) into the particleHandler.
 */
void BidispersedChute::insertLargeParticle(unsigned int& numberOfLargeToGo)
{
    SphericalParticle p0;
    p0.setSpecies(speciesHandler.getObject(1));
    p0.setRadius(parameters.getLargeParticleRadius() * distribution(generator));
    insertOneParticle(p0);
    numberOfLargeToGo--;
    logger(DEBUG, "Inserted large particle");
}

/*!
 * A subfunction of insertParticles; inserts a small particle (of species 2) into the particleHandler.
 */
void BidispersedChute::insertSmallParticle(unsigned int& numberOfSmallToGo)
{
    SphericalParticle p0;
    p0.setSpecies(speciesHandler.getObject(2));
    p0.setRadius(parameters.getSmallParticleRadius()* distribution(generator));
    insertOneParticle(p0);
    numberOfSmallToGo--;
    logger(DEBUG, "Inserted small particle");
}

/*!
 * A subfunction of insertParticles; inserts the given particle into the domain and the particleHandler.
 */
void BidispersedChute::insertOneParticle(BaseParticle& p0)
{
    Vec3D pos;
    p0.setVelocity(Vec3D(0.0, 0.0, 0.0));
    unsigned int toBeFailed = 1000;
    do
    {
        //if we failed to insert the particle 1000 times, increase domain in z direction
        if (toBeFailed == 0)
        {
            setZMax(getZMax() + 2 * parameters.getLargeParticleRadius());
            toBeFailed = 1000;
            logger(DEBUG, "resetting failed and increasing zmax, new zmax: %", getZMax());
        }
        pos.X = random.getRandomNumber(getXMin() + p0.getRadius(), getXMax() - p0.getRadius());
        pos.Y = random.getRandomNumber(getYMin() + p0.getRadius(), getYMax() - p0.getRadius());
        pos.Z = random.getRandomNumber(getZMin() + p0.getRadius(), getZMax());
        p0.setPosition(pos);
        toBeFailed--;
    } while (!checkParticleForInteraction(p0));
    particleHandler.copyAndAddObject(p0);
}

/*!
 * A subfunction of setupInitialConditions; sets periodic boundaries around the domain in x- and y-direction.
 */
void BidispersedChute::setBoundaries()
{
    PeriodicBoundary b0;
    b0.set(Vec3D(0.0, 1.0, 0.0), getYMin(), getYMax());
    boundaryHandler.copyAndAddObject(b0);
    if (isPeriodicInX)
    {
        b0.set(Vec3D(1.0, 0.0, 0.0), getXMin(), getXMax());
        boundaryHandler.copyAndAddObject(b0);
    }
}

/*!
 * A subfunction of the constructor; sets gravity vector, domain size, and the rough bottom,
 * which is made of particles to simulate a rough surface.
 */
void BidispersedChute::setChuteProperties()
{
    setChuteAngleAndMagnitudeOfGravity(parameters.getAngleInDegrees(), 1);
    setChuteLength(20);
    setChuteWidth(10);
    setInflowHeight(parameters.getInflowHeight());
    setZMax(4);
    setRoughBottomType(RoughBottomType::MULTILAYER);
    setFixedParticleRadius(parameters.getFixedParticleRadius());
}

/*!
 * A subfunction of the constructor; sets an appropriate time step (such that collisions are well-resolved),
 * and the frequency of writing data to output-files. Note, that these values can be overwritten in any application after constructing this class.
 */
void BidispersedChute::setTimeAndSaveCount()
{
    logger.assert_always(dynamic_cast<LinearViscoelasticSlidingFrictionSpecies*>(speciesHandler.getObject(2)) !=
                         nullptr, "cannot cast % to LinearViscoelasticSlidingFrictionSpecies",
                         speciesHandler.getObject(2)->getName());
    auto smallSpecies = *(dynamic_cast<LinearViscoelasticSlidingFrictionSpecies*>(speciesHandler.getObject(2)));
    const Mdouble massSmall = 4.0 / 3 * constants::pi * pow(parameters.getSmallParticleRadius(), 3.0)
                              * smallSpecies.getDensity();
    setTimeStep(smallSpecies.getCollisionTime(massSmall) / 50);
    logger(INFO, "time step %", getTimeStep());
    setTimeMax(20);
    setSaveCount(1e4);
}

/*!
 * A subfunction of the constructor; sets species properties such that all particles have the same density (6/pi),
 * and all collisions have the same contact time tc=sqrt(<d>/g)/200
 */
void BidispersedChute::setSpeciesProperties()
{
    const Mdouble density = 6.0 / constants::pi;
    
    const Mdouble massSmall = 4.0 / 3 * constants::pi * pow(parameters.getSmallParticleRadius(), 3.0) * density;
    const Mdouble massLarge = 4.0 / 3 * constants::pi * pow(parameters.getLargeParticleRadius(), 3.0) * density;
    logger(INFO, "mass large: %", massLarge);
    
    auto sReference = LinearViscoelasticSlidingFrictionSpecies();
    //for the reference particles (d=1, m=1, see silbert):
    sReference.setDensity(density);
    sReference.setDissipation(25); //gamma^n
    sReference.setSlidingDissipation(2.0 / 7 * sReference.getDissipation()); //  gamma^t
    sReference.setStiffness(2e5); // k^n
    sReference.setSlidingStiffness(2.0 / 7 * sReference.getStiffness()); // k^t
    sReference.setSlidingFrictionCoefficient(0.5); //mu
    
    const Mdouble r_c = sReference.getRestitutionCoefficient(1);
    const Mdouble tc_1 = sReference.getCollisionTime(1);
    
    //for the large particles
    auto sLarge = sReference;
    const Mdouble tc_l = std::sqrt(2 * parameters.getLargeParticleRadius()) * tc_1;
    sLarge.setCollisionTimeAndRestitutionCoefficient(tc_l, r_c, massLarge);
    sLarge.setSlidingDissipation(2.0 / 7 * sLarge.getDissipation());
    sLarge.setSlidingStiffness(2.0 / 7 * sLarge.getStiffness());
    sLarge.setSlidingFrictionCoefficient(0.5);
    logger(INFO, "restitution coefficient large %, collision time large %", r_c, tc_l);
    
    //for the small particles
    auto sSmall = sReference;
    const Mdouble tc_s = tc_1 * std::sqrt(parameters.getSmallParticleRadius());
    sSmall.setDensity(density);
    sSmall.setCollisionTimeAndRestitutionCoefficient(tc_s, r_c, massSmall);
    sSmall.setSlidingDissipation(2.0 / 7 * sSmall.getDissipation()); //  gamma^t
    sSmall.setSlidingStiffness(2.0 / 7 * sSmall.getStiffness()); // k^t
    sSmall.setSlidingFrictionCoefficient(0.5); //mu
    
    //add all species to handler. Note that this must be done after setting all properties, since otherwise the mixed
    //species don't make sense
    speciesHandler.copyAndAddObject(sReference);
    speciesHandler.copyAndAddObject(sLarge);
    speciesHandler.copyAndAddObject(sSmall);
}
