//Copyright (c) 2013-2020, The MercuryDPM Developers Team. All rights reserved.
//For the list of developers, see <http://www.MercuryDPM.org/Team>.
//
//Redistribution and use in source and binary forms, with or without
//modification, are permitted provid->d that the following conditions are met:
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provid->d with the distribution.
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

#include "PeriodicBoundary.h"
#include "ParticleHandler.h"
#include "Particles/BaseParticle.h"
#include "MpiDataClass.h"
#include "MpiContainer.h"
#include "DPMBase.h"
#include "MercuryBase.h"

/*!
 * \details constructor
 */
PeriodicBoundary::PeriodicBoundary()
        : BasePeriodicBoundary()
{
    distanceLeft_ = std::numeric_limits<double>::quiet_NaN();
    distanceRight_ = std::numeric_limits<double>::quiet_NaN();
    scaleFactor_ = std::numeric_limits<double>::quiet_NaN();

    logger(DEBUG, "PeriodicBoundary::PeriodicBoundary() finished");
}

/*!
 * \details destructor
 */
PeriodicBoundary::~PeriodicBoundary()
{
    logger(DEBUG, "PeriodicBoundary::~PeriodicBoundary() finished");
}

/*!
 * \details Copy method; creates a copy on the heap and returns its pointer.
 */
PeriodicBoundary* PeriodicBoundary::copy() const
{
    return new PeriodicBoundary(*this);
}

/*!
 * \details Copy constructor
 */
PeriodicBoundary::PeriodicBoundary(const PeriodicBoundary& other)
    : BasePeriodicBoundary(other)
{
    normal_ = other.normal_;
    scaleFactor_ = other.scaleFactor_;
    distanceLeft_ = other.distanceLeft_;
    distanceRight_ = other.distanceRight_;
    shift_ = other.shift_;
}

/*!
 * \details Defines the boundary, given a normal vector such that all particles 
 * are within {x: position_left<=normal*x<position_right}. The shift vector is set 
 * assuming that the domain is rectangular (shift parallel to normal).
 * \param[in] normal        The normal vector pointing from the left wall into the domain
 * \param[in] distanceLeft  The (signed) distance between the left wall and the origin
 * \param[in] distanceRight The (signed) distance between the right wall and the origin
 */
void PeriodicBoundary::set(Vec3D normal, Mdouble distanceLeft, Mdouble distanceRight)
{
    // factor is used to set normal to unit length
    scaleFactor_ = 1. / std::sqrt(Vec3D::dot(normal, normal));
    normal_ = normal * scaleFactor_;
    distanceLeft_ = distanceLeft * scaleFactor_;
    distanceRight_ = distanceRight * scaleFactor_;
    logger.assert_always(distanceRight_ > distanceLeft_,
            "PeriodicBoundary::set: left distance needs to be smaller than right distance");
    shift_ = normal_ * (distanceRight_ - distanceLeft_);
}

void PeriodicBoundary::set(Vec3D normal, Vec3D positionLeft, Vec3D positionRight)
{
    set(normal, Vec3D::dot(positionLeft, normal), Vec3D::dot(positionRight, normal));
}

/*!
 * \details like PeriodicBoundary::set(normal, distanceLeft, distanceRight), but 
 * including the possibility of setting the shift direction vector.
 * \param[in]   normal The normal vector pointing from the left wall into the domain
 * \param[in]   distanceLeft The (signed) distance between the left wall and the origin
 * \param[in]   distanceRight The (signed) distance between the right wall and the origin
 * \param[in]   shiftDirection The vector over which particles will be shifted when 
 *              moving through the PeriodicBoundary
 */
void PeriodicBoundary::set(Vec3D normal, Mdouble distanceLeft, Mdouble distanceRight, Vec3D shiftDirection)
{
    // factor is used to set normal to unit length
    scaleFactor_ = 1. / std::sqrt(Vec3D::dot(normal, normal));
    normal_ = normal * scaleFactor_;
    distanceLeft_ = distanceLeft * scaleFactor_;
    distanceRight_ = distanceRight * scaleFactor_;
    // factor is used to set shift vector to correct length
    scaleFactor_ = (distanceRight_ - distanceLeft_) * Vec3D::dot(shiftDirection, normal_);
    shift_ = shiftDirection * scaleFactor_;
}

/*!
 * \details Sets the shift_ vector through setting the planewise shift.
 * We delete the component of planewiseShift that is parallel to normal_. 
 */
void PeriodicBoundary::setPlanewiseShift(Vec3D planewiseShift)
{
    planewiseShift -= Vec3D::dot(planewiseShift, normal_) / Vec3D::dot(normal_, normal_) * normal_;
    shift_ = normal_ * (distanceRight_ - distanceLeft_) + planewiseShift;
}

/*!
 * \return The vector perpendicular to the periodic boundary
 */
Vec3D PeriodicBoundary::getNormal() const
{
    return normal_;
}

/*!
 * \return The distance of the left wall to the origin, in normal direction
 */
Mdouble PeriodicBoundary::getDistanceLeft() const
{
    return distanceLeft_;
}

/*!
 * \return The distance of the left wall to the origin, in normal direction
 */
Mdouble PeriodicBoundary::getDistanceRight() const
{
    return distanceRight_;
}

/*!
 * \return The vector going from the left to the right sid-> of the periodic boundary
 */
Vec3D PeriodicBoundary::getShift() const
{
    return shift_;
}

/*!
 * \details Allows the left periodic boundary to be moved to a new position and 
 * automatically changes its shift value
 * \param[in] distanceLeft  The distance (from the origin) to which the left 
 *                          boundary is moved
 */
void PeriodicBoundary::moveLeft(Mdouble distanceLeft)
{
    distanceLeft_ = distanceLeft * scaleFactor_;
    shift_ = normal_ * (distanceRight_ - distanceLeft_);
}

/*!
 * \details Allows the right periodic wall to be moved to a new position and 
 * automatically changes its shift value
 * \param[in] distanceRight     The distance (from the origin) to which the right 
 *                              boundary is moved
 */
void PeriodicBoundary::moveRight(Mdouble distanceRight)
{
    distanceRight_ = distanceRight * scaleFactor_;
    shift_ = normal_ * (distanceRight_ - distanceLeft_);
}

/*!
 * \details Returns the distance to the closest edge of the boundary to the particle.
 * Since this function should be called before calculating any Particle-Wall 
 * interactions, it can also be used to set the shift vector in case of curved walls.
 * Positive means that the particle is insid-> the periodic domain, negative means that it is
 * outsid-> the periodic domain.
 * \param[in] p     A reference to the particle which distance to the periodic 
 *                  boundary is calculated
 */
Mdouble PeriodicBoundary::getDistance(const BaseParticle& p) const
{
    return getDistance(p.getPosition());
}

/*!
 * \details Returns the distance to the edge closest to the position
 * \param[in] position  A reference to the position which distance to the periodic 
 *                      boundary is to be calculated
 */
Mdouble PeriodicBoundary::getDistance(const Vec3D& position) const
{
    Mdouble distanceFromPlaneThroughOrigin = Vec3D::dot(position, normal_);
    return std::min(distanceFromPlaneThroughOrigin - distanceLeft_, 
                    distanceRight_ - distanceFromPlaneThroughOrigin);
}

/*!
 * \details Shifts the particle either to the left or right, using the method isClosestToLeftBoundary to determine which
 * sid-> it should be shifted to.
 * \param[in] p         A pointer to the particle which will be shifted.
 */
void PeriodicBoundary::shiftPosition(BaseParticle* p) const
{
    if (isClosestToLeftBoundary(*p))
    {
        p->move(shift_);
    }
    else
    {
        p->move(-shift_);
    }
}

/*!
 * \details Shifts the particle either to the left or right, using the method isClosestToLeftBoundary to determine which
 * sid-> it should be shifted to.
 * \param[in] p         A pointer to the particle which will be shifted.
 */
void PeriodicBoundary::shiftPosition(Vec3D& p) const
{
    if (isClosestToLeftBoundary(p))
    {
        p += shift_;
    }
    else
    {
        p -= shift_;
    }
}


/*!
 * \details Shifts two given positions by the shift_ vector. 
 * \param[in] position1     The first position to be shifted
 * \param[in] position2     The second position to be shifted
 * \todo (AT) see toDo of PeriodicBoundary::shiftPosition().
 */
void PeriodicBoundary::shiftPositions(Vec3D& position1, Vec3D& position2) const
{
    if (isClosestToLeftBoundary(position1))
    {
        position1 += shift_;
        position2 += shift_;
    }
    else
    {
        position1 -= shift_;
        position2 -= shift_;
    }
}

/*
 * \details Returns TRUE if particle checked is closest to the 'left'
 * wall, and FALSE if it is closest to the 'right' wall. 
 * \param[in] p A point to a BaseParticle that is being checked.
 * \return      true if it is closest to the left boundary, false otherwise
 */
bool PeriodicBoundary::isClosestToLeftBoundary(const BaseParticle& p) const
{
    return isClosestToLeftBoundary(p.getPosition());
}

/*
 * \details Returns TRUE if position checked is closest to the 'left'
 * wall, and FALSE if it is closest to the 'right' wall.
 * \param[in] p A position vector p that is checked.
 * \return      true if it is closest to the left boundary, false otherwise
 */
bool PeriodicBoundary::isClosestToLeftBoundary(const Vec3D& p) const
{
    const Mdouble distance = Vec3D::dot(p, normal_);
    return (distanceRight_ - distance > distance - distanceLeft_);
}

/*!
 * \details Reads the boundary properties from an istream
 * \param[in] is        the istream
 */
void PeriodicBoundary::read(std::istream& is)
{
    BasePeriodicBoundary::read(is);
    std::string dummy;
    is >> dummy >> normal_
       >> dummy >> scaleFactor_
       >> dummy >> distanceLeft_
       >> dummy >> distanceRight_
       >> dummy >> shift_;
}

/*!
 * \details Deprecated version of read().
 * \deprecated Should be gone by Mercury 2.0. Instead, use CubeInsertionBoundary::read().
 */
void PeriodicBoundary::oldRead(std::istream& is)
{
    std::string dummy;
    is >> dummy >> normal_
       >> dummy >> scaleFactor_
       >> dummy >> distanceLeft_
       >> dummy >> distanceRight_
       >> dummy >> shift_;
}

/*!
 * \details Writes boundary's properties to an ostream
 * \param[in] os    the ostream
 */
void PeriodicBoundary::write(std::ostream& os) const
{
    BasePeriodicBoundary::write(os);
    os << " normal " << normal_
       << " scaleFactor " << scaleFactor_
       << " distanceLeft " << distanceLeft_
       << " distanceRight " << distanceRight_
       << " shift " << shift_;
}

/*!
 * \details Returns the name of the object class
 * \return      the object's class' name, i.e. 'CubeInsertionBoundary'
 */
std::string PeriodicBoundary::getName() const
{
    return "PeriodicBoundary";
}

/*!
 * \details Checks the distance of given particle to the closest of both periodic 
 * walls, and creates a periodic copy of the particle if needed (i.e. if the particle
 * is closer to the periodic wall than the radius of the largest particle in the
 * system).
 * \param[in] p         Particle to be checked and possibly periodically copied
 * \param[in,out] pH    System's ParticleHandler, (1) from which the interaction radius
 *                      of its largest particle is retrieved to determine the maximum 
 *                      distance from the wall at which a particle should still have
 *                      a periodic copy created, and (2) to which a possible periodic
 *                      copy of the particle will be added
 */
void PeriodicBoundary::createPeriodicParticle(BaseParticle* p, ParticleHandler& pH)
{
    //note that getDistance sets closestToLeftBoundary_ to true or false depending on which side is closest
    const Mdouble maxDistance = p->getMaxInteractionRadius() + pH.getLargestParticle()->getMaxInteractionRadius();
    if (getDistance(*p) < maxDistance)
    {
        createGhostParticle(p);
    }
}

void PeriodicBoundary::createGhostParticle(BaseParticle* pReal)
{
    ParticleHandler& pH = getHandler()->getDPMBase()->particleHandler;

    //Step 1: Copy the particle to new ghost particle.
    BaseParticle* pGhost = pReal->copy();

    //Step 2: Copy the interactions of the ghost particle.
    pGhost->copyInteractionsForPeriodicParticles(*pReal);

    //Step 3: Shift the ghost to the 'reflected' location.
    shiftPosition(pGhost);

    //Step 4: If Particle is double shifted, get correct original particle
    BaseParticle* from = pReal;
    while (from->getPeriodicFromParticle() != nullptr)
        from = from->getPeriodicFromParticle();
    pGhost->setPeriodicFromParticle(from);
    pGhost->setPeriodicGhostParticle(true);
    
    pH.addObject(pGhost);
}

/*!
 * \details Checks the distance of given particle to the closest of both periodic 
 * walls, and creates a periodic copy of the particle if needed (i.e. if the particle
 * is closer to the periodic wall than the radius of the largest particle in the
 * system).
 * \param[in] p         Particle to be checked and possibly periodically copied
 * \param[in,out] pH    System's ParticleHandler, (1) from which the interaction radius
 *                      of its largest particle is retrieved to determine the maximum 
 *                      distance from the wall at which a particle should still have
 *                      a periodic copy created, and (2) to which a possible periodic
 *                      copy of the particle will be added
 */
void PeriodicBoundary::createPeriodicParticles(ParticleHandler& pH)
{
#ifdef MERCURY_USE_MPI
    if (NUMBER_OF_PROCESSORS == 1)
    {
#endif
    unsigned numberOfParticles = pH.getSize();
    
    for (unsigned i = 0; i < numberOfParticles; i++)
    {
        createPeriodicParticle(pH.getObject(i), pH);
    }
#ifdef MERCURY_USE_MPI
    }
#endif
}

/*!
 * \details Loops through all particles to see if they have become ghosts. If that
 * is the case their position is shifted. 
 * Note: This is only for a serial build - periodic particles work different in paralle
 * \param[in]
 * \param[out] pH  the particle handler that contains all particles that
 * need to be checked 
 */
void PeriodicBoundary::checkBoundaryAfterParticlesMove(ParticleHandler& pH)
{
#ifdef MERCURY_USE_MPI
    if (NUMBER_OF_PROCESSORS == 1)
    {
#endif
    for (auto p = pH.begin(); p != pH.end(); ++p)
    {
        if (getDistance((*p)->getPosition()) < 0)
        {
            shiftPosition(*p);
            getHandler()->getDPMBase()->hGridUpdateMove(*p, shift_.getLengthSquared());
        }
    }
#ifdef MERCURY_USE_MPI
    }
#endif
}
















