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


#include "ConstantMassFlowMaserBoundary.h"
#include "DPMBase.h"
#include <algorithm>

/*!
 * \details MaserBoundary constructor, sets all scalars to non-sensible values.
 */
ConstantMassFlowMaserBoundary::ConstantMassFlowMaserBoundary()
{
#ifdef MERCURY_USE_MPI
    logger(WARN,"Make sure the whole maser boundary is on one processor when using parallel code");
#endif
    distanceLeft_ = std::numeric_limits<double>::quiet_NaN();
    distanceRight_ = std::numeric_limits<double>::quiet_NaN();
    gapSize_ = std::numeric_limits<double>::quiet_NaN();
    maserIsActivated_ = false;
}

/*!
 * \details MaserBoundary constructor from a periodic boundary. It copies the properties of the periodic boundary, and
 *          then immediately opens the maser. Do not insert particles for the maser after this constructor.
 * \param[in] periodicBoundary The periodic boundary on which this maser boundary is based.
 */
ConstantMassFlowMaserBoundary::ConstantMassFlowMaserBoundary(const PeriodicBoundary& periodicBoundary) : BaseBoundary(
        periodicBoundary)
{
    logger(INFO, "Constructor ConstantMassFlowMaserBoundary(const PeriodicBoundary&) started");
    distanceLeft_ = periodicBoundary.getDistanceLeft();
    distanceRight_ = periodicBoundary.getDistanceRight();
    normal_ = periodicBoundary.getNormal();
    shift_ = periodicBoundary.getShift();
    gapSize_ = std::numeric_limits<double>::quiet_NaN();
    maserIsActivated_ = false;
    
    logger(INFO, "going to activate the maser");
    //assume that the maser should be activated immediately when it gets constructed from a periodic boundary
    activateMaser();
}

/*!
 * \details Copy method, creates a copy of the object on the heap and returns a
 * pointer to it.
 * \return      pointer to the copy
 */
ConstantMassFlowMaserBoundary* ConstantMassFlowMaserBoundary::copy() const
{
    return new ConstantMassFlowMaserBoundary(*this);
}

/*!
 * \details Set all the properties of the boundary at once.
 * \param[in] normal            Normal unit vector of the (parallel) boundary walls
 * \param[in] distanceLeft      The distance of the left wall to the origin
 * \param[in] distanceRight     The distance of the right wall to the origin
 */
void ConstantMassFlowMaserBoundary::set(Vec3D normal, Mdouble distanceLeft, Mdouble distanceRight)
{
    // factor is used to set normal to unit length
    Mdouble scaleFactor_ = 1. / std::sqrt(Vec3D::dot(normal, normal));
    normal_ = normal * scaleFactor_;
    distanceLeft_ = distanceLeft * scaleFactor_;
    distanceRight_ = distanceRight * scaleFactor_;
    shift_ = normal_ * (distanceRight_ - distanceLeft_);
    maserIsActivated_ = false;
}

/*!
 * \details Reads the boundary properties from an istream
 * \param[in,out] is        the istream from which the boundary must be read
 */
void ConstantMassFlowMaserBoundary::read(std::istream& is)
{
    BaseBoundary::read(is);
    std::string dummy;
    is >> dummy >> normal_
       >> dummy >> distanceLeft_
       >> dummy >> distanceRight_
       >> dummy >> shift_
       >> dummy >> gapSize_
       >> dummy >> maserIsActivated_;
    unsigned int n;
    is >> dummy >> n;
    const SpeciesHandler& speciesHandler = getHandler()->getDPMBase()->speciesHandler;
    for (unsigned int i = 0; i < n; ++i)
    {
        unsigned int key;
        unsigned int value;
        is >> dummy >> key >> dummy >> value;
        speciesConversionNormalToMaser_[speciesHandler.getObject(key)] = speciesHandler.getObject(value);
        speciesConversionMaserToNormal_[speciesHandler.getObject(value)] = speciesHandler.getObject(key);
    }
    logger(DEBUG, "Finished reading ConstantMassFlowMaserBoundary. \nNormal: % \nDistanceLeft: % \nDistanceRight: % "
                  "\nGapsize: % \nMaserIsActivated: %", normal_, distanceLeft_, distanceRight_, gapSize_,
           maserIsActivated_);
}

/*!
 * \details Writes boundary's properties to an ostream
 * \param[in] os    the ostream to which the boundary must be written
 */
void ConstantMassFlowMaserBoundary::write(std::ostream& os) const
{
    BaseBoundary::write(os);
    os << " normal " << normal_
       << " distanceLeft " << distanceLeft_
       << " distanceRight " << distanceRight_
       << " shift " << shift_
       << " gapSize " << gapSize_
       << " maserIsActivated " << maserIsActivated_
       << " numberOfMaserSpecies " << speciesConversionMaserToNormal_.size();
    for (auto p : speciesConversionNormalToMaser_)
    {
        os << " outflowSpeciesIndex " << p.first->getIndex() << " maserSpeciesIndex " << p.second->getIndex();
    }
}

/*!
 * \details Returns the name of the object class
 * \return      the object's class' name, i.e. 'MaserBoundary'
 */
std::string ConstantMassFlowMaserBoundary::getName() const
{
    return "ConstantMassFlowMaserBoundary";
}

/*!
 * \details Shifts the particle (using the closestToLeftBoundary_ value)
 * \param[in] p         A pointer to the particle which will be shifted.
 */
void ConstantMassFlowMaserBoundary::shiftPosition(BaseParticle* const p) const
{
    if (isClosestToRightBoundary(p))
    {
        p->move(-shift_);
    }
    else // if closest to right boundary
    {
        p->move(shift_);
    }
}

/*!
 * \details Checks the distance of given particle to the closest of both
 * walls, and creates a periodic copy of the particle if needed (i.e. if the particle
 * is closer to the periodic wall than the radius of the largest particle in the
 * system). Note that this is the same function as in PeriodicBoundary, but with the
 * extra check to make sure that only maser particles are periodic particles
 * \param[in] p         Particle to be checked and possibly periodically copied
 * \param[in,out] pH    System's ParticleHandler, (1) from which the interaction radius
 *                      of its largest particle is retrieved to determine the maximum
 *                      distance from the wall at which a particle should still have
 *                      a periodic copy created, and (2) to which a possible periodic
 *                      copy of the particle will be added
 */
void ConstantMassFlowMaserBoundary::createPeriodicParticle(BaseParticle* p, ParticleHandler& pH)
{
    if (isMaserParticle(p))
    {
        // check if particle is near the boundaries of the maser domain
        if (getDistance(p) < p->getMaxInteractionRadius() + pH.getLargestParticleLocal()->getMaxInteractionRadius())
        {
            BaseParticle* pGhost = createGhostCopy(p);
            
            // shift to the periodic location
            shiftPosition(pGhost);
            
            // add the periodic particle to the handler
            pH.addObject(pGhost);
            
            //furthermore, if the particle is on the right it has to be copied over to the outflow domain
            if (isClosestToRightBoundary(p) && maserIsActivated_)
            {
                // Copy the particle and its interactions if it is a maser particle
                BaseParticle* pGhostOut = createGhostCopy(p);
                
                // shift to the outflow location
                pGhostOut->move(gapSize_ * normal_);
                // add the periodic particle to the handler
                pH.addObject(pGhostOut);
            }
        }
    }
}

void ConstantMassFlowMaserBoundary::createPeriodicParticles(ParticleHandler& pH)
{
    unsigned numberOfParticles = pH.getSize();
    for (unsigned i = 0; i < numberOfParticles; i++)
    {
        createPeriodicParticle(pH.getObject(i), pH);
    }
}

/*!
 * \details Helper function for createPeriodicParticles. It makes a copy of the particle p, and labels the copy as a
 * ghost particle.
 * \param[in] p Particle that needs to be copied
 * \return      Copy of particle p, labeled as a ghost particle
 */
BaseParticle* ConstantMassFlowMaserBoundary::createGhostCopy(BaseParticle* const p) const
{
    // Copy the particle and its interactions
    BaseParticle* pGhost = p->copy();
    pGhost->copyInteractionsForPeriodicParticles(*p);
    
    //Set the 'last' particle. If Particle is multiply shifted, get correct original particle
    BaseParticle* last = p;
    while (last->getPeriodicFromParticle() != nullptr)
    {
        last = last->getPeriodicFromParticle();
    }
    pGhost->setPeriodicFromParticle(last);
    return pGhost;
}

/*!
 * \details Checks whether a given particle (a) is in the Maser and (b) has
 * crossed the closest wall. If so, shifts its position so as to have it appear
 * at the other wall, and creates a 'real' equivalent in the outflow domain.
 * \param[in] p         The particle to be checked and possibly shifted and copied
 * \param pH            The ParticleHandler, which is unused in this implementation
 * \return              Returns true if particle p has interacted with the
 *                      boundary, false otherwise
 */
bool ConstantMassFlowMaserBoundary::checkBoundaryAfterParticleMoved(BaseParticle* p, ParticleHandler& pH)
{
    // check if particle passed either of the boundary walls
    if (isMaserParticle(p) && getDistance(p) < 0)
    {
        // Checks if the particle is closest to the right boundary.
        // If so, and if the Maser is turned on, then create a 'real'
        // equivalent and move it over to the outflow domain
        if (isClosestToRightBoundary(p) && isCopying())
        {
            BaseParticle* pCopy = p->copy();
            pCopy->setSpecies(speciesConversionMaserToNormal_.find(p->getSpecies())->second);
            pCopy->move(gapSize_ * normal_);
            pH.addObject(pCopy);
        }
        
        // If the (original) particle has crossed a boundary wall (be it left or right),
        // then shift that particle periodically.
        shiftPosition(p);
    }
    else //remove normal particles trying to enter the gap between the maser and the outflow domain
    {
        const Mdouble distance = Vec3D::dot(p->getPosition(), normal_);
        if (distance < gapSize_ + distanceRight_ && isNormalParticle(p))
        {
            pH.removeObject(p->getIndex());
            return true;
        }
    }
    return false;
}

void ConstantMassFlowMaserBoundary::checkBoundaryAfterParticlesMove(ParticleHandler& pH)
{
    for (auto p = pH.begin(); p != pH.end(); ++p)
    {
        //If the particle has interacted with the boundary, change the iterator
        if (checkBoundaryAfterParticleMoved(*p, pH))
        {
            p--;
        }
    }
}

/*!
 * \details Turns given particle into a 'maser particle' by changing its species into
 * a 'maser particle' copy species. If the particle species is not yet in the std::map
 * speciesConversionNormalToMaser_, it and its maser copy species are added.
 * This function should be called at the beginning of the simulation, right after actually
 * filling the maser with particles, flagging each particle as one belonging to the
 * maser.
 * \param[in,out] p     The particle which is added to the maser. Its species is 'changed'
 *                      to the maser copy species.
 */
void ConstantMassFlowMaserBoundary::addParticleToMaser(BaseParticle* p)
{
    // check if particle species already known by the maser
    auto conversion = speciesConversionNormalToMaser_.find(p->getSpecies());
    if (conversion != speciesConversionNormalToMaser_.end())
    {
        //species known and flagged (i.e. 'converted')
        logger(VERBOSE, "[ConstantMassFlowMaserBoundary::addParticleToMaser()] Species conversion already present");
        p->setSpecies(conversion->second);
    }
    else
    {
        SpeciesHandler& speciesHandler = (getHandler()->getDPMBase()->speciesHandler);
        // species is not yet known by the maser, it has to be added to both maps
        ParticleSpecies* newSpecies = speciesHandler.copyAndAddObject(*(p->getSpecies()));
        speciesConversionNormalToMaser_.insert(
                std::pair<const ParticleSpecies*, const ParticleSpecies*>(p->getSpecies(), newSpecies));
        speciesConversionMaserToNormal_.insert(
                std::pair<const ParticleSpecies*, const ParticleSpecies*>(newSpecies, p->getSpecies()));
        logger(INFO, "[ConstantMassFlowMaserBoundary::addParticleToMaser()] New species conversion created");
        logger(INFO, "Original species ID: %, new species ID: %", p->getSpecies()->getId(), newSpecies->getId());
        
        //Copy over the mixed species. The delete is necessary here since it is overwritten with a copy of the old mixed
        //species, and otherwise the properties are not copied over correctly.
        //
        //setId and setIndex refer to the two different species which are present in this mixed species.
        //The highest species-index should appear first andis therefore the ID, while the second species-index is the
        //"index" of the mixed species.
        for (const BaseSpecies* const s : speciesHandler)
        {
            if (s->getId() != newSpecies->getId() && s->getId() != p->getSpecies()->getId())
            {
                BaseSpecies* newMixed = speciesHandler.getMixedObject(s->getId(), newSpecies->getId());
                const BaseSpecies* const oldMixed = speciesHandler.getMixedObject(s->getId(), p->getSpecies()->getId());
                oldMixed->copyInto(newMixed);
                newMixed->setId(newSpecies->getId());
                newMixed->setIndex(s->getId());
                logger(DEBUG, "mixed species of % and % is now \n %, should be \n %", s->getId(), newSpecies->getId(),
                       *newMixed, *oldMixed);
            }
        }
        
        // now the species IS added, so flag (convert) it!
        p->setSpecies(newSpecies);
    }
    p->move(-normal_ * gapSize_);
}

void ConstantMassFlowMaserBoundary::removeParticleFromMaser(BaseParticle* p)
{
    auto conversion = speciesConversionMaserToNormal_.find(p->getSpecies());
    if (conversion != speciesConversionMaserToNormal_.end())
    {
        p->setSpecies(conversion->second);
        p->move(normal_ * gapSize_);
    }
}

/*!
 * \details One checks if the particle is in the Maser, by checking whether its
 * species is found in the list of maser particle species
 */
bool ConstantMassFlowMaserBoundary::isMaserParticle(BaseParticle* p) const
{
    // Check if the particle is in the Maser, by checking whether its
    // species is found in the list of maser particle species
    auto positionMaserSpeciesInMap = speciesConversionMaserToNormal_.find(p->getSpecies());
    return (positionMaserSpeciesInMap != speciesConversionMaserToNormal_.end()); //Test if it is a maser particle
}

/*!
 * \details One checks if the particle is a normal particle, by checking whether
 * the species changes when we do speciesConversionNormalToMaser_().
 */
bool ConstantMassFlowMaserBoundary::isNormalParticle(BaseParticle* p) const
{
    auto toFindOutflowSpecies = speciesConversionNormalToMaser_.find(p->getSpecies());
    return (toFindOutflowSpecies != speciesConversionNormalToMaser_.end());
}

/*!
 * \details The maser boundary should be activated before the time loop starts, but after all particles are added.
 */
void ConstantMassFlowMaserBoundary::actionsBeforeTimeLoop()
{
    activateMaser();
}

/*!
 * \details Activate the maser by creating a gap betweeen the periodic domain and the outflow domain, and transforming
 * all the particles in the domain to maser particles.
 * For an explanation of the magic number 6, see the detailed documentation for gapSize_. Note that gapSize_ is set here
 * and not in set, since it is possible that the user first adds boundaries to the domain before the particles. The
 * maser boundary should only be activated after all particles are in.
 * \todo JMFT: Would it be helpful to be able to set a custom gapSize ?
 */
void ConstantMassFlowMaserBoundary::activateMaser()
{
    if (!maserIsActivated_)
    {
        logger(INFO, "Going to add particles to the maser and shift the periodic maser boundaries");
        ParticleHandler& pH = getHandler()->getDPMBase()->particleHandler;
        gapSize_ = 6 * pH.getLargestParticleLocal()->getRadius();
        logger(INFO, "just before particle loop");
        for (BaseParticle* const p : pH)
        {
            if (getDistance(p) > 0)
            {
                addParticleToMaser(p);
            }
        }
        
        
        distanceLeft_ -= gapSize_;
        distanceRight_ -= gapSize_;
        maserIsActivated_ = true;

        turnOnCopying();

        logger(INFO, "The Maser has been activated, and is now copying particles.");
    }
    else
    {
        logger(WARN, "Cannot activate the maser boundary twice!");
    }
}

///\details closing the maser: do exactly the opposite of activating the maser: close the gap by moving all particles,
///and assign the "normal" species to all maser particles again.
void ConstantMassFlowMaserBoundary::closeMaser()
{
    if (maserIsActivated_)
    {
        for (BaseParticle* const p : getHandler()->getDPMBase()->particleHandler)
        {
            if (getDistance(p) > 0)
            {
                removeParticleFromMaser(p);
            }
        }
        distanceLeft_ += gapSize_;
        distanceRight_ += gapSize_;
        maserIsActivated_ = false;

        turnOffCopying();

        logger(INFO, "The Maser has been deactivated.");
    }
    else
    {
        logger(WARN, "Cannot close the maser if it is not active");
    }
}

bool ConstantMassFlowMaserBoundary::isActivated() const
{
    return maserIsActivated_;
}

void ConstantMassFlowMaserBoundary::turnOnCopying() 
{
    maserIsCopying_ = true;
}

void ConstantMassFlowMaserBoundary::turnOffCopying()
{
    maserIsCopying_ = false;
}


bool ConstantMassFlowMaserBoundary::isCopying() const
{
    if (!isActivated())
        logger(WARN, "The Maser is closed, so isCopying() is irrelevant.");

    return maserIsCopying_;
}

Mdouble ConstantMassFlowMaserBoundary::getDistanceLeft() const
{
    return distanceLeft_;
}

Mdouble ConstantMassFlowMaserBoundary::getDistanceRight() const
{
    return distanceRight_;
}

Mdouble ConstantMassFlowMaserBoundary::getGapSize() const
{
    return gapSize_;
}
