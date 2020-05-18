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

#ifndef RESTRICTEDWALL_H
#define RESTRICTEDWALL_H

#include "BaseWall.h"
#include "Math/Vector.h"
#include "InfiniteWall.h"

/*!
 * \brief Restriction of a wall to the intersection with another wall
 */

class RestrictedWall : public BaseWall
{
public:
    
    /*!
     * \brief Default constructor, the normal is infinitely long.
     */
    RestrictedWall();
    
    /*!
     * \brief Copy constructor, copy the given wall.
     */
    RestrictedWall(const RestrictedWall& w);
    
    
    /*!
     * \brief Constructor setting values.
     */
    RestrictedWall(BaseWall* wall, InfiniteWall* restriction);
    
    /*!
     * \brief Default destructor.
     */
    ~RestrictedWall() override;
    
    /*!
     * \brief Wall copy method. It calls the copy constructor of this Wall, useful for polymorphism
     */
    RestrictedWall* copy() const override;
    
    /*!
     * \brief Defines a standard wall, given an outward normal vector s.t. normal*x=normal*point for all x of the wall.
     */
    void set(BaseWall* wall, InfiniteWall* restriction);
    
    /*!
     * \brief Compute the distance from the wall for a given BaseParticle and return if there is a collision. If there is a collision, also return the normal vector.
     */
    bool getDistanceAndNormal(const BaseParticle& p, Mdouble& distance, Vec3D& normal_return) const override;
    
    /*!
     * \brief Reads RestrictedWall from a restart file.
     */
    void read(std::istream& is) override;
    
    /*!
     * \brief Reads RestrictedWall from an old-style restart file.
     */
    void oldRead(std::istream& is);
    
    /*!
     * \brief Writes the RestrictedWall to an output stream, usually a restart file.
     */
    void write(std::ostream& os) const override;
    
    /*!
     * \brief Returns the name of the object, in this case the string "RestrictedWall".
     */
    std::string getName() const override;
    
    /*!
     * \brief Look up the interaction between this wall and a BaseParticle at a certain timeStamp.
     */
    BaseInteraction*
    getInteractionWith(BaseParticle* p, unsigned timeStamp, InteractionHandler* interactionHandler) override;
    
    void writeVTK(VTKContainer& vtk) const override;

private:
    /*!
     * Outward normal vector. This does not have to be a unit vector.
     */
    BaseWall* wall_;
    InfiniteWall* restriction_;
};

#endif
