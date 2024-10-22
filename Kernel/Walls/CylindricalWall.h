//Copyright (c) 2013-2023, The MercuryDPM Developers Team. All rights reserved.
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

#ifndef CYLINDRICALWALL_H
#define CYLINDRICALWALL_H

#include "BaseWall.h"

/*!
 * \class CylindricalWall
 * \brief
 * \deprecated This class will be gone in Mercury 2.0, use 
 * AxisymmetricIntersectionOfWalls instead.
 */
class CylindricalWall : public BaseWall
{
public:
    /*!
     * \brief
     */
    CylindricalWall();
    
    /*!
     * \brief
     */
    CylindricalWall(const CylindricalWall& p);
    
    /*!
     * \brief
     */
    explicit CylindricalWall(double radius);
    
    /*!
     * \brief Wall copy method. It calls the copy contrustor of this Wall, usefull for polymorfism
     */
    CylindricalWall* copy() const override;
    
    /*!
     * \brief Defines a standard wall, given an outward normal vector s. t. normal*x=position
     */
    void set(Mdouble radius);
    
    /*!
     * \brief Since this function should be called before calculating any Particle-Wall interactions, it can also be used to set the normal vector in case of curved walls.
     */
    bool getDistanceAndNormal(const BaseParticle& P, Mdouble& distance, Vec3D& normal_return) const override;
    
    /*!
     * \brief reads wall
     */
    void read(std::istream& is) override;
    
    /*!
     * \brief
     */
    void oldRead(std::istream& is);
    
    /*!
     * \brief outputs wall
     */
    void write(std::ostream& os) const override;
    
    /*!
     * \brief Returns the name of the object
     */
    std::string getName() const override;
    
    /*!
     * \brief access function for radius
     */
    double getRadius() const;

private:
    Mdouble radius_;///
};

#endif
