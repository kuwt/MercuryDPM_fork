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

#include "MpiDataClass.h"
#include "DPMBase.h"

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
void copyDataFromMPIParticleToParticle(MPIParticle* bP, BaseParticle* p)
{
    
    //Set important quantities
    p->setId(bP->id);
    p->setRadius(bP->radius);
    p->setPosition(bP->position);
    p->setAngularVelocity(bP->angularVelocity);
    p->setVelocity(bP->velocity);
    p->setOrientation(bP->orientation);
    p->setCommunicationComplexity(bP->communicationComplexity);
    
//    p->setAxes(bP->axes);
//    p->setExponents(bP->epsilon1, bP->epsilon2);
    
    //Set HGrid values
    p->setHGridNextObject(nullptr);
    p->setHGridPrevObject(nullptr);
    p->setHGridX(99999);
    p->setHGridY(99999);
    p->setHGridZ(99999);
    p->setHGridLevel(bP->HGridLevel);
    
    //This is not a periodic particle
    p->setPeriodicFromParticle(nullptr);
    p->setInPeriodicDomain(false);
    
    //Fix maser if it is maser
    p->setMaserParticle(bP->isMaser);
    
    //Fixed particles need to be fixed again
    if (bP->isFixed)
    {
        p->fixParticle();
    }
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
    copyDataFromMPIParticleToParticle(bP, p);
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
    
    bP.id = p->getId();
    bP.indSpecies = p->getIndSpecies();
    bP.radius = p->getRadius();
//    bP.axes = p->getAxes();
//    bP.epsilon1 = p->getExponentEps1();
//    bP.epsilon2 = p->getExponentEps2();
    bP.position = p->getPosition();
    bP.angularVelocity = p->getAngularVelocity();
    bP.velocity = p->getVelocity();
    bP.orientation = p->getOrientation();
    bP.HGridLevel = p->getHGridLevel();
    bP.communicationComplexity = p->getCommunicationComplexity();
    bP.isMaser = p->isMaserParticle();
    bP.isFixed = p->isFixed();
    return bP;
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

