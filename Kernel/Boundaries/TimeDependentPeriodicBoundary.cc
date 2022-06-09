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

#include "TimeDependentPeriodicBoundary.h"
#include "ParticleHandler.h"
#include "Particles/BaseParticle.h"
#include "MpiDataClass.h"
#include "MpiContainer.h"
#include "DPMBase.h"
#include <functional>

/*!
 * \details constructor
 */
TimeDependentPeriodicBoundary::TimeDependentPeriodicBoundary()
        : BasePeriodicBoundary()
{
    distanceLeft_ = std::numeric_limits<double>::quiet_NaN();
    distanceRight_ = std::numeric_limits<double>::quiet_NaN();
    maxShift_ = 0;


#ifdef MERCURY_USE_MPI
    MPIContainer& communicator = MPIContainer::Instance();
    if (communicator.getNumberOfProcessors() > 1)
    {
	    logger(WARN,"LeesEdwardsBoundaries are currently not implemented in parallel MercuryDPM");
    }
#endif
    logger(DEBUG, "TimeDependentPeriodicBoundary::TimeDependentPeriodicBoundary() finished");
}

/*!
 * \details destructor
 */
TimeDependentPeriodicBoundary::~TimeDependentPeriodicBoundary()
{
    logger(DEBUG, "TimeDependentPeriodicBoundary::~TimeDependentPeriodicBoundary() finished");
}

/*!
 * \details Copy method; creates a copy on the heap and returns its pointer.
 */
TimeDependentPeriodicBoundary* TimeDependentPeriodicBoundary::copy() const
{
    return new TimeDependentPeriodicBoundary(*this);
}

/*!
 * \details Copy constructor
 */
TimeDependentPeriodicBoundary::TimeDependentPeriodicBoundary(const TimeDependentPeriodicBoundary& other)
{
    normal_ = other.normal_;
    distanceLeft_ = other.distanceLeft_;
    distanceRight_ = other.distanceRight_;
    planewiseShift_ = other.planewiseShift_;
    boost_ = other.boost_;
    maxShift_ = other.maxShift_;

}

void TimeDependentPeriodicBoundary::set(Vec3D normal, Mdouble distanceLeft, Mdouble distanceRight,
             std::function<Vec3D(Mdouble)> planewiseShift, std::function<Vec3D(Mdouble)> boost)
{
    normal_ = normal;
    distanceLeft_ = distanceLeft, 
    distanceRight_ = distanceRight, 
    planewiseShift_ = planewiseShift;
    boost_ = boost;
    maxShift_ = maxShift_;
}

void TimeDependentPeriodicBoundary::set(Vec3D normal, Vec3D positionLeft, Vec3D positionRight,
             std::function<Vec3D(Mdouble)> planewiseShift, std::function<Vec3D(Mdouble)> boost)
{
    set(normal, Vec3D::dot(positionLeft,normal), Vec3D::dot(positionRight,normal), 
            planewiseShift, boost);
}

void TimeDependentPeriodicBoundary::setPlanewiseShiftAndBoost(
            std::function<Vec3D(Mdouble)> planewiseShift, std::function<Vec3D(Mdouble)> boost)
{
    planewiseShift_ = planewiseShift;
    boost_ = boost;
}

void TimeDependentPeriodicBoundary::setMaxShift(Mdouble maxShift)
{
    maxShift_ = maxShift;
}

/*!
 * \return The vector perpendicular to the periodic boundary
 */
Vec3D TimeDependentPeriodicBoundary::getNormal() const
{
    return normal_;
}

/*!
 * \return The distance of the left wall to the origin, in normal direction
 */
Mdouble TimeDependentPeriodicBoundary::getDistanceLeft() const
{
    return distanceLeft_;
}

/*!
 * \return The distance of the left wall to the origin, in normal direction
 */
Mdouble TimeDependentPeriodicBoundary::getDistanceRight() const
{
    return distanceRight_;
}

Vec3D TimeDependentPeriodicBoundary::getShift(Mdouble time) const
{
    return getPlanewiseShift(time) + normal_ * (distanceRight_ - distanceLeft_);
}

Vec3D TimeDependentPeriodicBoundary::getPlanewiseShift(Mdouble time) const
{
    if (maxShift_ == 0)
        return planewiseShift_(time);
    if (maxShift_ > 0)
    {
        Vec3D p = planewiseShift_(time);
        Mdouble m = p.getLength();
        Vec3D n = p / m;
        return fmod(m, maxShift_) * n;
    }
    if (maxShift_ < 0)
        logger(ERROR, "[TimeDependentPeriodicBoundary::getPlanewiseShift] maxShift_ = % is negative", maxShift_);
}

Vec3D TimeDependentPeriodicBoundary::getBoost(Mdouble time) const
{
    return boost_(time);
}


/*!
 * \details Allows the left periodic boundary to be moved to a new position and 
 * automatically changes its shift value
 * \param[in] distanceLeft  The distance (from the origin) to which the left 
 *                          boundary is moved
 */
void TimeDependentPeriodicBoundary::moveLeft(Mdouble distanceLeft)
{
    distanceLeft_ = distanceLeft;
}

/*!
 * \details Allows the right periodic wall to be moved to a new position and 
 * automatically changes its shift value
 * \param[in] distanceRight     The distance (from the origin) to which the right 
 *                              boundary is moved
 */
void TimeDependentPeriodicBoundary::moveRight(Mdouble distanceRight)
{
    distanceRight_ = distanceRight;
}

Mdouble TimeDependentPeriodicBoundary::getDistance(const BaseParticle& p) const
{
    return getDistance(p.getPosition());
}

/*!
 * \details Returns the distance to the edge closest to the position
 * \param[in] position  A reference to the position which distance to the periodic 
 *                      boundary is to be calculated
 */
Mdouble TimeDependentPeriodicBoundary::getDistance(const Vec3D& position) const
{
    Mdouble distanceFromPlaneThroughOrigin = Vec3D::dot(position, normal_);
    return std::min(distanceFromPlaneThroughOrigin - distanceLeft_, 
                    distanceRight_ - distanceFromPlaneThroughOrigin);
}

/*!
 * \details Shouldn't be used for TimeDependentPeriodicBoundary.
 * Instead, use TimeDependentPeriodicBoundary::shiftandBoostParticle. 
 */
void TimeDependentPeriodicBoundary::shiftPosition(BaseParticle* p) const
{
}

void TimeDependentPeriodicBoundary::shiftAndBoostParticle(BaseParticle* p, Mdouble time) const
{
    if (isClosestToLeftBoundary(*p))
    {
        p->move(getShift(time));
        p->addVelocity(getBoost(time));
    }
    else
    {
        p->move(-getShift(time));
        p->addVelocity(-getBoost(time));
    }
}

void TimeDependentPeriodicBoundary::shiftPositions(Vec3D& position1, Vec3D& position2) const
{
    if (isClosestToLeftBoundary(position1))
    {
        position1 += getShift(0); // TODO JMFT: ?!?!?!
        position2 += getShift(0); 
    }
    else
    {
        position1 -= getShift(0);
        position2 -= getShift(0);
    }
}

/*
 * \details Returns TRUE if particle checked is closest to the 'left'
 * wall, and FALSE if it is closest to the 'right' wall. 
 * \param[in] p A point to a BaseParticle that is being checked.
 * \return      true if it is closest to the left boundary, false otherwise
 */
bool TimeDependentPeriodicBoundary::isClosestToLeftBoundary(const BaseParticle& p) const
{
    return isClosestToLeftBoundary(p.getPosition());
}

/*
 * \details Returns TRUE if position checked is closest to the 'left'
 * wall, and FALSE if it is closest to the 'right' wall.
 * \param[in] p A position vector p that is checked.
 * \return      true if it is closest to the left boundary, false otherwise
 */
bool TimeDependentPeriodicBoundary::isClosestToLeftBoundary(const Vec3D& p) const
{
    const Mdouble distance = Vec3D::dot(p, normal_);
    return (distanceRight_ - distance > distance - distanceLeft_);
}

void TimeDependentPeriodicBoundary::createPeriodicParticles(ParticleHandler& pH)
{
#ifdef MERCURY_USE_MPI
    if (NUMBER_OF_PROCESSORS == 1)
    {
#endif
        unsigned numberOfParticles = pH.getSize();

        for(unsigned i = 0; i < numberOfParticles; i++)
        {
            createPeriodicParticle(pH.getObject(i),pH);
        }
#ifdef MERCURY_USE_MPI
    }
#endif
}

void TimeDependentPeriodicBoundary::createGhostParticle(BaseParticle *pReal)
{
    ParticleHandler& pH = getHandler()->getDPMBase()->particleHandler;

    //Step 1: Copy the particle to new ghost particle.
    BaseParticle* pGhost = pReal->copy();

    //Step 2: Copy the interactions of the ghost particle.
    pGhost->copyInteractionsForPeriodicParticles(*pReal);

    //Step 3: Shift the ghost to the 'reflected' location.
    shiftAndBoostParticle(pGhost, pH.getDPMBase()->getTime());

    //Step 4: If Particle is double shifted, get correct original particle
    BaseParticle* from = pReal;
    while (from->getPeriodicFromParticle() != nullptr)
        from = from->getPeriodicFromParticle();
    pGhost->setPeriodicFromParticle(from);

    pH.addObject(pGhost);
}

void TimeDependentPeriodicBoundary::createPeriodicParticle(BaseParticle* p, ParticleHandler& pH)
{
    //note that getDistance sets closestToLeftBoundary_ to true or false depending on which side is closest
    if (getDistance(*p) < p->getMaxInteractionRadius() + pH.getLargestParticle()->getMaxInteractionRadius())
    {
        createGhostParticle(p);
    }
}

/*!
 * \details Loops through all particles to see if they have become ghosts. If that
 * is the case their position is shifted. 
 * Note: This is only for a serial build - periodic particles work different in paralle
 * \param[in]
 * \param[out] pH  the particle handler that contains all particles that
 * need to be checked 
 */
void TimeDependentPeriodicBoundary::checkBoundaryAfterParticlesMove(ParticleHandler& pH)
{
#ifdef MERCURY_USE_MPI
    if (NUMBER_OF_PROCESSORS == 1)
    {
#endif
        for (auto p = pH.begin(); p != pH.end(); ++p)
        {
            if (getDistance((*p)->getPosition()) < 0)
            {
                shiftAndBoostParticle(*p, pH.getDPMBase()->getTime());
            }
        }
#ifdef MERCURY_USE_MPI
    }
#endif
}

/*!
 * \details Reads the boundary properties from an istream
 * \param[in] is        the istream
 */
void TimeDependentPeriodicBoundary::read(std::istream& is)
{
    BasePeriodicBoundary::read(is);
    std::string dummy;
    is >> dummy >> normal_
       >> dummy >> distanceLeft_
       >> dummy >> distanceRight_
      // >> dummy >> planewiseShift_
      // >> dummy >> boost_
      ;
}

/*!
 * \details Writes boundary's properties to an ostream
 * \param[in] os    the ostream
 */
void TimeDependentPeriodicBoundary::write(std::ostream& os) const
{
    BasePeriodicBoundary::write(os);
    os << " normal " << normal_
       << " distanceLeft " << distanceLeft_
       << " distanceRight " << distanceRight_
       // << " planewiseShift " << planewiseShift_ 
       // << " boost " << boost_
       ;
}

/*!
 * \details Returns the name of the object class
 * \return      the object's class' name, i.e. 'CubeInsertionBoundary'
 */
std::string TimeDependentPeriodicBoundary::getName() const
{
    return "TimeDependentPeriodicBoundary";
}

