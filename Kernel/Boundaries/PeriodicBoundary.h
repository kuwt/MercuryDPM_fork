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

#ifndef PeriodicBoundary_H
#define PeriodicBoundary_H

#include "BasePeriodicBoundary.h"
#include "Math/Vector.h"
#include "MpiDataClass.h"

class ParticleHandler;

    /*!
     * \brief Defines a pair of periodic walls. Inherits from BaseBoundary.
     * \details The particles are in {x: position_left<=normal*x <position_right}, 
     * with normal being the outward unit normal vector of the right wall. If a 
     * particle moves outside these boundaries, it will be shifted.
     */

class PeriodicBoundary : public BasePeriodicBoundary
{
public:
    /*!
     * \brief default constructor
     */
    PeriodicBoundary();
    
    /*!
     * \brief destructor
     */
    ~PeriodicBoundary();
    
    /*!
     * \brief copy method
     */
    PeriodicBoundary* copy() const override;
    
    /*!
     * \brief Defines a periodic wall
     */
    void set(Vec3D normal, Mdouble distanceLeft, Mdouble distanceRight);

    void set(Vec3D normal, Vec3D positionLeft, Vec3D positionRight);

    /*!
   * \brief For general parallelogramic domains, the direction of the shift vector can to be set manually.
   */
    void set(Vec3D normal, Mdouble distanceLeft, Mdouble distanceRight, Vec3D shiftDirection);
    
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
    Vec3D getShift() const;
    
    /*!
     * \brief Sets the distance from the origin of the 'left' periodic wall
     */
    void moveLeft(Mdouble distanceLeft);
    
    /*!
     * \brief Sets the distance from the origin of the 'right' periodic wall
     */
    void moveRight(Mdouble distanceRight);
    
    /*!
     * \brief Returns the distance of the wall to the particle
     */
    Mdouble getDistance(const BaseParticle& p) const;
    
    /*!
     * \brief Returns the distance of the wall to the position
     */
    Mdouble getDistance(const Vec3D &position) const;
    
    /*!
     * \brief shifts the particle 
     * \param[in] p A pointer to the particle which will be shifted.
     */
    virtual void shiftPosition(BaseParticle* p) const;
    
    /*!
     * \brief shifts two positions
     */
    virtual void shiftPositions(Vec3D &postition1, Vec3D &postion2) const;

    /*!
     * \brief Returns TRUE if particle checked is closest to the 'left' 
     * wall, and FALSE if it is closest to the 'right' wall
     */
    virtual bool isClosestToLeftBoundary(const BaseParticle& p) const;

    /*
     * \details Returns TRUE if position checked is closest to the 'left'
     * wall, and FALSE if it is closest to the 'right' wall.
     */
    virtual bool isClosestToLeftBoundary(const Vec3D& p) const;

    /*!
     * \brief reads boundary properties from istream
     */
    virtual void read(std::istream& is) override;
   
    /*!
     * \brief deprecated version of CubeInsertionBoundary::read().
     */
    MERCURY_DEPRECATED
    void oldRead(std::istream& is);
    
    /*!
     * \brief writes boundary properties to ostream
     */
    virtual void write(std::ostream& os) const override;
    
    /*!
     * \brief Returns the name of the object
     */
    virtual std::string getName() const override;
    
    /*!   
     * \brief Checks distance of particle to closest wall and creates periodic 
     * copy if necessary
     */
    virtual void createPeriodicParticles(ParticleHandler &pH) override;

    /*!
     * \brief Creates and adds a ghost particle from a give real particle
     */
    void createGhostParticle(BaseParticle *pReal);

    /*!
     * \brief Creates a single periodic particle if required from a given particle
     */
    void createPeriodicParticle(BaseParticle* p, ParticleHandler &pH) override;

    /*!
     * \brief Checks if particle has crossed either boundary wall, and applies a shift
     * if that is the case.
     */
    virtual void checkBoundaryAfterParticlesMove(ParticleHandler& pH);
 
protected:
  /*!
   * \brief outward unit normal vector for right wall
   */
    Vec3D normal_;
  /*!
   * \brief position of left wall, s.t. normal*x=position_left
   */
    Mdouble distanceLeft_;
  /*!
   * \brief position of right wall, s.t. normal*x=position_right
   */
    Mdouble distanceRight_; 
  /*!
   * \brief This is the normal to rescale the normal vector to a unit vectors.
   */
    Mdouble scaleFactor_;  
  /*!
   * \brief shift from left to right boundary
   */
    Vec3D shift_;
    
};
#endif
