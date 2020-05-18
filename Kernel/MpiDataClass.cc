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

#include "MpiDataClass.h"
#include "DPMBase.h"
#include "Particles/LiquidFilmParticle.h"
#include "Particles/SphericalParticle.h"
#include "Particles/SuperQuadricParticle.h"
#include "Logger.h"

/*!
 * \brief Copies data from a BaseParticle to an MPIParticle class and returns this
 * \details In order to create ghost particles on other processors, data of particles
 * have to be transmitted to other processors. Only the required data is
 * sent. The data is sent in an MPIParticle data class and this function copies
 * the data from a particle into that class.
 * \param[in] p Pointer to a base particle from which data is copied
 * \return MPIParticle class is returned filled with data from BaseParticle p
 */

/*!
 * \brief Copies data from an MPIParticle class to a BaseParticle
 * \details When processors recieve data of ghost particles they have to add,
 * they recieve them in the MPIParticle class. This function turns the MPIParticle
 * class back into a BaseParticle and initialises hGrid values to "default"
 * \param[in] bP Pointer to an MPIParticle which contains data for a ghost particle
 * \param[in,out] p Pointer to BaseParticle, a ghost particle that will be added to the domain
 */
void MPISphericalParticle::copyDataFromMPIParticleToParticle(BaseParticle* p)
{
    
    //Set important quantities
    p->setId(id);
    p->setRadius(radius);
    p->setPosition(position);
    p->setAngularVelocity(angularVelocity);
    p->setVelocity(velocity);
    p->setOrientation(orientation);
    p->setCommunicationComplexity(communicationComplexity);

    //Set HGrid values
    p->setHGridNextObject(nullptr);
    p->setHGridPrevObject(nullptr);
    p->setHGridX(99999);
    p->setHGridY(99999);
    p->setHGridZ(99999);
    p->setHGridLevel(HGridLevel);
    
    //This is not a periodic particle
    p->setPeriodicFromParticle(nullptr);
    p->setInPeriodicDomain(false);
    
    //Fix maser if it is maser
    p->setMaserParticle(isMaser);
    
    //Fixed particles need to be fixed again
    if (isFixed)
    {
        p->fixParticle();
    }
}

void MPISuperQuadric::copyDataFromMPIParticleToParticle(BaseParticle* p)
{
    MPISphericalParticle::copyDataFromMPIParticleToParticle(p);
    p->setAxes(axes);
    p->setExponents(epsilon1, epsilon2);
}

void MPILiquidFilmParticle::copyDataFromMPIParticleToParticle(BaseParticle* p)
{
    MPISphericalParticle::copyDataFromMPIParticleToParticle(p);
    static_cast<LiquidFilmParticle*>(p)->setLiquidVolume(liquidVolume);
}

/*!
 * \brief Copies data from an MPIParticle class to a BaseParticle and sets the particleHandler and species
 * \param[in] bP Pointer to an MPIParticle which contains data for a ghost particle
 * \param[in,out] p Pointer to BaseParticle, a ghost particle that will be added to the domain
 * \param[in] particleHandler Pointer to the ParticleHandler required for creating a new particle
 * \todo MX: Maybe renamet his function to setParticleSpecies() or something
 */
void copyDataFromMPIParticleToParticle(MPIParticle* bP, BaseParticle* p, ParticleHandler* particleHandler)
{
    //Set the species of the particle, but before we can do that we have to set the handler
    p->setHandler(particleHandler);
    //p->setIndSpecies(bP->indSpecies);
    const ParticleSpecies* species = p->getHandler()->getDPMBase()->speciesHandler.getObject(bP->indSpecies);
    p->setSpecies(species);
    bP->copyDataFromMPIParticleToParticle(p);
}

/*!
 * \brief Copies data from a SuperQuadricParticle to an MPIParticle class and returns this
 * \details In order to create ghost particles on other processors, data of particles
 * have to be transmitted to other processors. Only the required data is
 * sent. The data is sent in an MPIParticle data class and this function copies
 * the data from a particle into that class.
 * \param[in] p Pointer to a SuperQuadricParticle particle from which data is copied
 * \return MPIParticle class is returned filled with data from BaseParticle p
 */
MPIParticle copyDataFromParticleToMPIParticle(BaseParticle* p)
{
    MPIParticle bP;
    bP.copyDataFromParticleToMPIParticle(p);
    return bP;
}

void MPISphericalParticle::copyDataFromParticleToMPIParticle(BaseParticle* p) {
    id = p->getId();
    indSpecies = p->getIndSpecies();
    radius = p->getRadius();
    position = p->getPosition();
    angularVelocity = p->getAngularVelocity();
    velocity = p->getVelocity();
    orientation = p->getOrientation();
    HGridLevel = p->getHGridLevel();
    communicationComplexity = p->getCommunicationComplexity();
    isMaser = p->isMaserParticle();
    isFixed = p->isFixed();
}

void MPISuperQuadric::copyDataFromParticleToMPIParticle(BaseParticle* p) {
    MPISphericalParticle::copyDataFromParticleToMPIParticle(p);
    axes = p->getAxes();
    epsilon1 = p->getExponentEps1();
    epsilon2 = p->getExponentEps2();
}

void MPILiquidFilmParticle::copyDataFromParticleToMPIParticle(BaseParticle* p) {
    MPISphericalParticle::copyDataFromParticleToMPIParticle(p);
    liquidVolume = static_cast<LiquidFilmParticle*>(p)->getLiquidVolume();
}

BaseParticle* MPISphericalParticle::newParticle () {
    return new SphericalParticle;
}

BaseParticle* MPISuperQuadric::newParticle () {
    return new SuperQuadricParticle;
}

BaseParticle* MPILiquidFilmParticle::newParticle () {
    return new LiquidFilmParticle;
}

/*!
 * \brief Copies the position from a particle to an MPIParticlePosition class
 * \param[in] particle BaseParticle which position is copied
 * \return MPIParticlePosition object containing the position of the particle
 */
MPIParticlePosition copyPositionFrom(BaseParticle* particle)
{
    MPIParticlePosition particlePosition;
    particlePosition.id = particle->getId();
    particlePosition.position = particle->getPosition();
    particlePosition.orientation = particle->getOrientation();
    if (std::is_base_of<MPILiquidFilmParticle,MPIParticle>())
        particlePosition.liquidVolume = static_cast<LiquidFilmParticle*>(particle)->getLiquidVolume();
    return particlePosition;
}

/*!
 * \brief Copies the velocity from a particle to an MPIParticleVelocity class
 * \param[in] particle BaseParticle which velocity is copied
 * \return MPIParticleVelocityn object containing the velocity of the particle
 */
MPIParticleVelocity copyVelocityFrom(BaseParticle* particle)
{
    MPIParticleVelocity particleVelocity;
    particleVelocity.velocity = particle->getVelocity();
    particleVelocity.angularVelocity = particle->getAngularVelocity();
    
    return particleVelocity;
}

Vec3D getMPISum(Vec3D& val)
{
#ifdef MERCURY_USE_MPI
    //Sum up over all domains
    Vec3D valGlobal = {0.0,0.0,0.0};
    MPIContainer& communicator = MPIContainer::Instance();
    communicator.allReduce(val.X, valGlobal.X, MPI_SUM);
    communicator.allReduce(val.Y, valGlobal.Y, MPI_SUM);
    communicator.allReduce(val.Z, valGlobal.Z, MPI_SUM);
    return valGlobal;
#else
    return val;
#endif
}

double getMPISum(double val)
{
#ifdef MERCURY_USE_MPI
    //Sum up over all domains
    double valGlobal = 0.0;
    MPIContainer& communicator = MPIContainer::Instance();
    communicator.allReduce(val, valGlobal, MPI_SUM);
    return valGlobal;
#else
    return val;
#endif
}
