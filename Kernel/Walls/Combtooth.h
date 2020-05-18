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

#ifndef COMBTOOTH_H
#define COMBTOOTH_H

#include "BaseWall.h"
#include "Math/Vector.h"
#include <cmath>

class Combtooth : public BaseWall
{
public:
    /*!
     * \brief Default constructor
     */
    Combtooth();
    
    /*!
     * \brief Copy constructor
     */
    Combtooth(const Combtooth& ct);
    
    /*!
     * \brief Default destructor
     */
    ~Combtooth() override;
    
    /*!
     * \brief Set
     */
    void set(Vec3D axis, Vec3D position, Mdouble radius);
    
    /*!
     * \brief Copy
     */
    Combtooth* copy() const override;
    
    bool getDistanceAndNormal(const BaseParticle& p,
                              Mdouble& distance, Vec3D& normal_return) const override;
    
    BaseInteraction* getInteractionWith(BaseParticle* p,
                                                     unsigned timeStamp,
                                                     InteractionHandler* interactionHandler) override;
    
    void read(std::istream& is) override;
    
    void write(std::ostream& os) const override;
    
    std::string getName() const override;

private:
    Vec3D axis_; // unit vector pointing in direction of axis
    Vec3D position_; // position vector of a point that the axis goes through
    Mdouble radius_;
    
};

#endif
