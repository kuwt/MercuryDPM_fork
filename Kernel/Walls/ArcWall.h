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

#ifndef ARCWALL_H
#define ARCWALL_H

#include "Walls/BaseWall.h"
#include <cmath>

class WallHandler;

class BaseParticle;

class BaseWall;

/*!
 * \brief A wall that is the inside of an arc of a cylinder.
 * \details The ArcWall is specified by the cylinder's axis (in turn by a
 * position vector and a direction vector), its radius, a 'centreline direction'
 * (a unit vector normal to the direction vector) for the centreline of the arc,
 * and a semiangle. 
 * For a wall that is the outside of an arc, use AxisymmetricIntersectionOfWalls
 * or Combtooth instead.
 */
class ArcWall : public BaseWall
{
public:
    /*!
     * \brief Default constructor
     */
    ArcWall();
    
    /*!
     * \brief Copy constructor
     */
    ArcWall(const ArcWall& aw);
    
    /*!
     * \brief Default destructor
     */
    ~ArcWall() override = default;
    
    /*!
     * \brief Set
     */
    void set(Vec3D axis, Vec3D pos, Mdouble radius, Vec3D centreline, Mdouble semiangle);
    
    ArcWall* copy() const override;
    
    bool getDistanceAndNormal(const BaseParticle& p,
                              Mdouble& distance, Vec3D& normal_return) const override;
    
    BaseInteraction* getInteractionWith(BaseParticle* p,
                                                     unsigned timeStamp,
                                                     InteractionHandler* interactionHandler) override;
    
    void read(std::istream& is) override;
    
    void write(std::ostream& os) const override;
    
    std::string getName() const override;

private:
    Vec3D axis_;
    Vec3D pos_;
    Mdouble radius_;
    Vec3D centreline_;
    Mdouble semiangle_;
    
};

#endif
