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


#ifndef CONSTANTMASSFLOWMASERBOUNDARY_H
#define CONSTANTMASSFLOWMASERBOUNDARY_H

#include <map>

#include "Boundaries/BaseBoundary.h"
#include "Boundaries/PeriodicBoundary.h"
#include "Interactions/BaseInteraction.h"
#include "Math/Vector.h"
#include "Particles/BaseParticle.h"

class ParticleSpecies;

/*!
 * \brief Variation on the PeriodicBoundary which also has an outflow part
 * \details Creates a boundary which divides the domain into two parts: a
 * periodic domain and an outflow domain, with a gap inbetween. 
 * Any particle flowing through the right of the periodic domain is
 * moved to both the left side of the periodic domain (as with a PeriodicBoundary), 
 * and also copied the left side of the outflow domain. Furthermore, the
 * particles near the right side of the periodic domain also exert forces on the
 * left side of the outflow domain as if there were no gap. 
 * However, the outflow domain has no influence on the periodic domain. 
 * When an outflow particle tries to
 * enter the gap between the periodic domain and the outflow domain, it gets
 * removed.
 *
 * There are two distinct properties of the state of the Maser. `Activated' or
 * `closed' refer to whether the particles within the Maser have been set to the
 * Maser species and whether they have been moved by the gapSize. Particles will
 * be affected by the Maser (i.e. moved periodically) iff the Maser is
 * activated.  On the other hand, `copying' refers to whether new particles are
 * being produced by the Maser. Thus, a Maser that is activated but not copying
 * behaves just as a periodic boundary.
 *
 * By default, activateMaser() also turns on copying. And when a Maser is
 * closed, the value of isCopying() is irrelevant, although it is set to false
 * for housekeeping.
 *
 * \todo Add functionality which allows for opening the maser boundary after a
 * certain time, being a normal periodic boundary until then \todo Consider
 * writing a destructor that closes the gap again \todo Consider splitting it in
 * 2 DPMBase instances, one for the periodic domain and one for the outflow
 * domain \todo Consider re-using the PeriodicBoundary by adding it as a data
 * member
 *
 * The difference between SubcriticalMaserBoundary (formerly known as
 * MaserBoundaryOldStyle) and ConstantMassFlowMaserBoundary (formerly known
 * simply as MaserBoundary) is that in ConstantMassFlowMaserBoundary, the left
 * hand side of the periodic box does not have any influence on the rest of the
 * flow, but the right side of the periodic box and the left side of the outflow
 * domain interact.  The ConstantMassFlowMaserBoundary is most useful for fast
 * (supercritical) flows, and for flows for which the flux across the boundary
 * needs to be controlled. The SubcriticalMaserBoundary is more useful for slow
 * flows, as the ConstantMassFlowMaserBoundary might generate "pulse-waves" in
 * those cases. 
 *
 * For a compact overview of the behaviour of ConstantMassFlowMaserBoundary,
 * please look at the output of ConstantMassFlowMaserSelfTest.
 *
 * \todo Which Maser is it used in Denissen2019?  To cite the Maser: I. F. C.
 * Denissen, T. Weinhart, A. Te Voortwis, S. Luding, J. M. N. T. Gray and A. R.
 * Thornton, Bulbous head formation in bidisperse shallow granular flow over an
 * inclined plane. Journal of Fluid Mechanics, 866:263--297, mar 2019.
 */
class ConstantMassFlowMaserBoundary : public BaseBoundary { public:
    /*!  \brief MaserBoundary constructor
     */
    ConstantMassFlowMaserBoundary();
    
    /*!  \brief Maserboundary constructor that takes a periodic boundary, and
     * converts it to a maser boundary
     */
    explicit ConstantMassFlowMaserBoundary(const PeriodicBoundary&
            periodicBoundary);
    
    /*!  \brief Creates a copy of this maser on the heap.
     */
    ConstantMassFlowMaserBoundary* copy() const override;
    
    /*!  \brief Sets all boundary properties at once and adds particles of the
     * handler to the maser.
     */
    void set(Vec3D normal, Mdouble distanceLeft, Mdouble distanceRight);
    
    /*!  \brief reads boundary properties from istream
     */
    void read(std::istream& is) override;
    
    /*!  \brief writes boundary properties to ostream
     */
    void write(std::ostream& os) const override;
    
    /*!  \brief Returns the name of the object
     */
    std::string getName() const override;
    
    /*!  \brief Creates periodic particles when the particle is a maser particle
     * and is sufficiently close to one of the boundary walls.
     */
    void createPeriodicParticle(BaseParticle* p, ParticleHandler& pH) override;
    
    void createPeriodicParticles(ParticleHandler& pH) override;
    
    /*!  \brief Shifts the particle to its 'periodic' position if it is a maser
     * particle and has crossed either of the walls. Creates a 'normal' particle
     * at its current position if it is a maser particle which crossed the RIGHT
     * boundary wall.
     */
    bool checkBoundaryAfterParticleMoved(BaseParticle* p, ParticleHandler& pH);
    
    /*!  \brief Evaluates what the particles have to do after they have changed
     * position
     */
    void checkBoundaryAfterParticlesMove(ParticleHandler& pH) override;
    
    /*!  \brief Converts a 'normal' particle into a maser particle.
     */
    void addParticleToMaser(BaseParticle* p);
    
    /*!  \brief Convert a maser particle into a 'normal' particle
     */
    void removeParticleFromMaser(BaseParticle* p);
    
    /*!  \brief Returns true if the particle is a Maser particle, and false
     * otherwise.
     */
    bool isMaserParticle(BaseParticle* p) const;
    
    /*!  \brief Returns true if the particle is a Normal particle, and false
     * otherwise.
     */
    bool isNormalParticle(BaseParticle* p) const;
    
    /*!  \brief Does everything that needs to be done for this boundary between
     * setupInitialConditions and the time loop, in this case, it activates the
     * maser.
     */
    void actionsBeforeTimeLoop() override;
    
    /*!  \brief Opens the gap, and transforms particles to maser particles. Also
     * calls turnOnCopying().
     */
    void activateMaser();
    
    /*!  \brief Stops copying particles (and act merely as a chute)
     */
    void closeMaser();

    /*! 
     * \brief Returns whether the Maser is activated or not. 
     */
    bool isActivated() const;

    /*!
     * \brief Start copying particles.
     */
    void turnOnCopying();

    /*! 
     * \brief Stop copying particles.
     */
    void turnOffCopying();

    /*!
     * \brief Returns whether the Maser is copying particles or not.
     */
    bool isCopying() const;
    
    Mdouble getDistanceLeft() const;
    
    Mdouble getDistanceRight() const;
    
    Mdouble getGapSize() const;

private:
    
    /*!
     * \brief Shifts the particle to its 'periodic' position
     */
    void shiftPosition(BaseParticle* p) const;
    
    /*!
     * \brief Creates a copy of the input particle, that gets removed again in DPMBase::removeDuplicatePeriodicParticles
     */
    BaseParticle* createGhostCopy(BaseParticle* p) const;
    
    /*!
     * \brief Returns whether the given particle is closer to the right boundary of the periodic part.
     * \param[in] p Particle for which we would like to know whether it is closest to the right boundary
     * \return      True if p is closer to the right boundary, false otherwise
     */
    bool isClosestToRightBoundary(const BaseParticle* const p) const
    {
        const Mdouble distance = Vec3D::dot(p->getPosition(), normal_);
        return (distanceRight_ - distance < distance - distanceLeft_);
    }
    
    /*!
     * \brief Returns the distance of the wall to the particle
     * \param[in] p     Pointer to the particle of which we want to know the distance to the wall to
     * \return          Distance of the particle to the boundary: positive means the particle is inside the periodic
     *                  part of the boundary, negative means it's outside.
     */
    Mdouble getDistance(BaseParticle* p) const
    {
        const Mdouble distance = Vec3D::dot(p->getPosition(), normal_);
        return std::min(distance - distanceLeft_, distanceRight_ - distance);
    }
    
    /*!
     * \brief Normal unit vector of both maser walls. Points in the flowing direction.
     */
    Vec3D normal_;
    
    /*!
     * \brief position of left boundary wall, s.t. normal*x=position_left
     */
    Mdouble distanceLeft_;
    
    /*!
     * \brief position of right boundary wall, s.t. normal*x=position_right
     */
    Mdouble distanceRight_;
    
    /*!
     * \brief distance between the right side of the periodic domain and the start of the outflow domain.
     * \details I.e., each particle in the maser is moved
     * -distanceToOutflowDomain_ * normal when it becomes part of the maser, and
     *  moved distanceToOutflowDomain_ * normal when it is inserted in the
     *  outflow domain. Generally this is 6 times the radius of the largest
     *  particle, so that the ghost particles do not touch each other: the
     *  centre of the ghostparticle is at most 1 diameter (2 radii) away from
     *  the boundary, so a ghost particle can extend at most 3 particle radii
     *  away from the domain. Do this on both sides, and it follows that the gap
     *  should be at least 6 diameters wide.
     * \todo JMFT: Do you mean 6 radii?
     */
    Mdouble gapSize_;
    
    /*!
     * \brief Direction in which particles are to be shifted when they cross the boundary.
     * \details I.e., the vector pointing from a point the left boundary wall to the equivalent point
     * on the right one.
     */
    Vec3D shift_;
    
    /*!
     * \brief List of 'normal' particles' species, and their maser counterparts
     */
    std::map<const ParticleSpecies*, const ParticleSpecies*> speciesConversionNormalToMaser_;
    
    /*!
     * \brief List of 'maser' particles' species, and their normal counterparts
     */
    std::map<const ParticleSpecies*, const ParticleSpecies*> speciesConversionMaserToNormal_;
    
    /*!
     * \brief Flag whether or not the gap is created and particles transformed already.
     */
    bool maserIsActivated_;

    /*!
     * \brief Flag whether or not the Maser is copying particles.
     */
    bool maserIsCopying_;
    
};

#endif
