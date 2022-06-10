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

#ifndef TimeDependentPeriodicBoundary_H
#define TimeDependentPeriodicBoundary_H

#include "BasePeriodicBoundary.h"
#include "ParticleHandler.h"
#include "Math/Vector.h"
#include "MpiDataClass.h"
#include <functional>

/*!
 * \brief Class which creates a boundary with Lees-Edwards type periodic boundary conditions. 
 * \details A TimeDependentPeriodicBoundary is like a PeriodicBoundary, but, when a
 * particle crosses an edge,
 * <ul>
 *   <li>it is shifted as well as being copied</li>
 *   <li>it is given a boost</li>
 * </ul>
 * In general, the shift and the boost may depend on time in an arbitrary way.
 * They are to be specified as std::function<Mdouble (Mdouble)>.
 *
 * This sort of boundary is useful for studying shear flows. 
 * \details See also Lees and Edwards (J. Phys. C 1921, 
 * <a href="http://dx.doi.org/1088/0022-3719/5/15/006">doi:1088/0022-3719/5/15/006</a>). 
 * Inherits from BaseBoundary.
 * \todo Add link to paper by Lees-Edwards in the documentation of this class.
 * \todo Is implemented for 2D only now. Needs extension to 3D.
 */
class TimeDependentPeriodicBoundary : public BasePeriodicBoundary
{
public:

    /*!
     * \brief default constructor
     */
    TimeDependentPeriodicBoundary();

    /*!
     * \brief destructor
     */
    ~TimeDependentPeriodicBoundary();

    /*!
     * \brief copy method
     */
    TimeDependentPeriodicBoundary* copy() const override;

    /*!
     * \brief copy constructor
     */
    TimeDependentPeriodicBoundary(const TimeDependentPeriodicBoundary& other);
    
    /*!
     * \brief Defines a TimeDependentPeriodicBoundary by its normal and positions, and by
     * the shifting and boosting that it does (projected onto the planewise
     * direction)
     * \details 
     * \param[in] normal  Vector specifying the normal direction of the edges
     * \param[in] distanceLeft Position of the first edge
     * \param[in] distanceRight Position of the second edge
     * \param[in] shift   Vector (projected to remove normal component) by which 
     *                    particles crossing an edge are to be shifted, as a
     *                    function of time.
     * \param[in] boost   Vector (projected to remove normal component) by which
     *                    particles crossing an edge are to be boosted, as a 
     *                    function of time.
     */
    void set(Vec3D normal, Mdouble distanceLeft, Mdouble distanceRight,
             std::function<Vec3D(Mdouble)> planewiseShift, std::function<Vec3D(Mdouble)> boost);

    /*!
     * \brief As above, but by specifying two positions that the boundaries go
     * through instead of distanceLeft and distanceRight.
     */
    void set(Vec3D normal, Vec3D positionLeft, Vec3D positionRight,
             std::function<Vec3D(Mdouble)> shift, std::function<Vec3D(Mdouble)> boost);
    
    /*!
     * \brief Set the planewise shift and boost (projected onto the planewise
     * direction) as functions of time. (Boost should be the derivative of
     * shift)
     */
    void setPlanewiseShiftAndBoost(
            std::function<Vec3D(Mdouble)> shift, std::function<Vec3D(Mdouble)> boost);

    /*!
     * \brief Set the maximum shift (will take fmod w.r.t. this)
     */
    void setMaxShift(Mdouble maxShift);

    /*!
     * \brief returns the vector normal to the periodic boundary
     */
    Vec3D getNormal() const;

    /*!
     * \brief Returns the distance of the left wall to the origin, in normal direction
     */
    Mdouble getDistanceLeft() const;
    
    /*!
     * \brief Returns the distance of the right wall to the origin, in normal direction
     */
    Mdouble getDistanceRight() const;

    /*!
     * \brief Returns the vector going from the left to the right side of the periodic boundary
     */
    Vec3D getShift(Mdouble time) const;

    /*!
     * \brief Returns the planewise shift as a function of time
     */
    Vec3D getPlanewiseShift(Mdouble time) const;

    /*!
     * \brief Returns the planewise boost as a function of time
     */
    Vec3D getBoost(Mdouble time) const;

    /*!
     * \brief Sets the distance from the origin of the 'left' periodic wall
     */
    void moveLeft(Mdouble distanceLeft);

    /*!
     * \brief Sets the distance from the origin of the 'right' periodic wall
     */
    void moveRight(Mdouble distanceRight);
    
    /*!
     * \brief Returns the distance of the edge to the particle
     */
    Mdouble getDistance(const BaseParticle& p) const override;

    /*!
     * \brief Returns the distance of the edge to the position
     */
    Mdouble getDistance(const Vec3D &position) const override;

    /*!
     * \brief shifts <em>and boosts</em> the particle
     * \param[in] p A pointer to the particle which will be shifted and boosted.
     * \todo{JMFT: The time comes from p->getHandler()->getDPMBase()->getTime(), 
     * which will be undefined if p does not belong to a handler.}
     */
    virtual void shiftPosition(BaseParticle* p) const override;

    virtual void shiftAndBoostParticle(BaseParticle* p, Mdouble time) const;

    /*!
     * \brief shifts two positions (JMFT: Why, what is this for?)
     */
    virtual void shiftPositions(Vec3D &position1, Vec3D &position2) const;

    /*!
     * \brief Returns true if particle checked is closer to the 'left'
     * edge, and false if it is closer to the 'right' edge
     */
    virtual bool isClosestToLeftBoundary(const BaseParticle& p) const;

    /*!
     * \brief Returns true if position checked is closer to the 'left'
     * edge, and false if it is closer to the 'right' edge
     */
    virtual bool isClosestToLeftBoundary(const Vec3D& p) const override;

    /*!
     * \brief Checks distance of particle to closer edge and creates a periodic
     * copy if necessary
     */
    virtual void createPeriodicParticles(ParticleHandler &pH) override;

    /*!
     * \brief Creates and adds a ghost particle from a given real particle
     * \todo{JMFT: The time comes from p->getHandler()->getDPMBase()->getTime(), 
     * which will be undefined if p does not belong to a handler.}
     */
    void createGhostParticle(BaseParticle *pReal);

    /*!
     * \brief Creates a single periodic particle if required from a given particle
     */
    void createPeriodicParticle(BaseParticle* p, ParticleHandler &pH) override;

    /*!
     * \brief Loops over particles, checks if each particle has crossed either
     * boundary edge, and applies a shift if that is the case.
     */
    virtual void checkBoundaryAfterParticlesMove(ParticleHandler& pH) override;

    /*!
     * \brief reads boundary properties from istream
     */
    virtual void read(std::istream& is) override;

    /*!
     * \brief writes boundary properties to ostream
     */
    virtual void write(std::ostream& os) const override;

    /*!
     * \brief Returns the name of the object
     */
    virtual std::string getName() const override;
    
protected:
    /*!
     * \brief outward unit normal vector for right edge
     */
    Vec3D normal_;

    /*!
     * \brief position of left edge, s.t. normal*x = distanceLeft_
     */
    Mdouble distanceLeft_;

    /*!
     * \brief position of right edge, s.t. normal*x = distanceRight_
     */
    Mdouble distanceRight_; 

    /*!
     * \brief shift from left to right boundary in the planewise direction
     * (Note: total shift = planewiseShift_(time) + normal * (distanceRight_ - distanceLeft_)
     */
    std::function<Vec3D (Mdouble)> planewiseShift_;

    /*!
     * \brief boost from the left to right boundary
     */
    std::function<Vec3D (Mdouble)> boost_; 

    /*!
     * \brief Maximum shifting (will take fmod(shift, maxshift) )
     */
    Mdouble maxShift_;

};
#endif
