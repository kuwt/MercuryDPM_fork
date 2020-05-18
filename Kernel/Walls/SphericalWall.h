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

///This is a class defining walls. It defines the 
///interaction of regular walls and periodic walls
///with particles as defined in Particle
///Modifications:

#ifndef SphericalWall_H
#define SphericalWall_H

#include "BaseWall.h"
#include "Math/Vector.h"

/*!
 * \brief A infinite wall fills the half-space {point: (position_-point)*normal_<=0}. 
 * \details Thus, the surface of the wall is a plane through position position_
 * with normal_ the outward unit normal vector of the wall 
 * (pointing away from the particles, into the wall).
 * Please note that this wall is infinite and straight.
 * 
 * A particle touches an infinite wall if (position_-point)*normal_<=radius.
 */

class SphericalWall : public BaseWall
{
public:
    
    /*!
     * \brief Default constructor, the normal is infinitely long.
     */
    SphericalWall();
    
    /*!
     * \brief Copy constructor, copy the given wall.
     */
    SphericalWall(const SphericalWall& w);
    
    
    /*!
     * \brief Constructor setting values.
     */
    SphericalWall(Mdouble radius, const ParticleSpecies* species);
    
    /*!
     * \brief Default destructor.
     */
    ~SphericalWall() override;
    
    /*!
     * \brief Wall copy method. It calls the copy constructor of this Wall, useful for polymorphism
     */
    SphericalWall* copy() const override;
    
    /*!
     * \brief Defines a spherical wall with radius r.
     */
    void setRadius(Mdouble radius);
    
    using BaseWall::move;
    
    /*!
     * \brief Returns the distance of the wall to the particle.
     */
    Mdouble getDistance(const Vec3D& otherPosition) const;
    
    /*!
     * \brief Returns the distance of the wall to the particle.
     */
    Mdouble getRadius() const;
    
    /*!
     * \brief Compute the distance from the wall for a given BaseParticle and return if there is a collision. If there is a collision, also return the normal vector.
     */
    bool getDistanceAndNormal(const BaseParticle& p, Mdouble& distance, Vec3D& normal_return) const override;
    
    /*!
     * \brief Reads SphericalWall from a restart file.
     */
    void read(std::istream& is) override;
    
    /*!
     * \brief Writes the SphericalWall to an output stream, usually a restart file.
     */
    void write(std::ostream& os) const override;
    
    /*!
     * \brief Returns the name of the object, in this case the string "SphericalWall".
     */
    std::string getName() const override;

private:
    /*!
     * This is the factor used to rescale the normal given by the user to a unit
     * vector. It is only used by the deprecated function move(Mdouble).
     */
    Mdouble radius_;
};

#endif
