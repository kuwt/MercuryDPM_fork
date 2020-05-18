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

#include "SimpleDrumSuperquadrics.h"
#include "Particles/BaseParticle.h"
#include "WallHandler.h"
#include "DPMBase.h"

SimpleDrumSuperquadrics::SimpleDrumSuperquadrics()
{
    logger(DEBUG, "SimpleDrumSuperquadrics() finished");
    radius_ = 1;
    wall.set({1,0,0}, {1,0,0});
}

/*!
 * \param[in] other The AxisymmetricIntersectionOfWalls that must be copied.
 */
SimpleDrumSuperquadrics::SimpleDrumSuperquadrics(const SimpleDrumSuperquadrics& other)
        : BaseWall(other)
{
    wall = other.wall;
    radius_ = other.radius_;
    logger(DEBUG, "AxisymmetricIntersectionOfWalls(const AxisymmetricIntersectionOfWalls &p) finished");
}

SimpleDrumSuperquadrics::~SimpleDrumSuperquadrics()
{
    logger(DEBUG, "SimpleDrumSuperquadricsuadrics() finished.");
}

/*!
 * \param[in] other The AxisymmetricIntersectionOfWalls that must be copied.
 */
SimpleDrumSuperquadrics&
SimpleDrumSuperquadrics::operator=(const SimpleDrumSuperquadrics& other)
{
    if (this == &other)
    {
        return *this;
    }
    else
    {
        return *(other.copy());
    }
}

/*!
 * \return pointer to a IntersectionOfWalls object allocated using new.
 */
SimpleDrumSuperquadrics* SimpleDrumSuperquadrics::copy() const
{
    return new SimpleDrumSuperquadrics(*this);
}

/*!
 * \details First, the particle is translated by the vector position_, then the 
 * distance normal and tangential to the orientation is computed. This normal 
 * and tangential direction is interpreted as the x and z coordinate. With the 
 * particle shifted into the XZ plane, the distance and normal is computed, as 
 * if the AxisymmetricIntersectionOfWalls would be a simple IntersectionOfWalls.
 * Finally, the object and the normal is rotated back to the original position.
 * 
 * See also AxisymmetricIntersectionOfWalls for details.
 */
bool SimpleDrumSuperquadrics::getDistanceAndNormal(const BaseParticle& p, Mdouble& distance,
                                                           Vec3D& normalReturn) const
{
    Vec3D normalDirection = p.getPosition();
    normalDirection.Y = 0;
    normalDirection.normalise();
    wall.set(normalDirection, normalDirection * radius_);
    logger(DEBUG, "setting wall to %, %", normalDirection, normalDirection * radius_);
    logger(DEBUG, "angular velocity of wall: %", getAngularVelocity());
    //determine wall distance, normal and contact in axissymmetric coordinates
    //and transform from axisymmetric coordinates
    return  (wall.getDistanceAndNormal(p, distance, normalReturn));
}

bool SimpleDrumSuperquadrics::getDistanceNormalOverlapSuperquadric(const SuperQuadricParticle& p, Mdouble& distance, Vec3D& normal_return,
                                          Mdouble& overlap) const
{
    Vec3D normalDirection = p.getPosition();
    normalDirection.Y = 0;
    if (normalDirection.getLengthSquared() < 1e-2)
        return false;
    normalDirection.normalise();
    wall.set(normalDirection, normalDirection * radius_);
    logger(DEBUG, "setting wall to [%], [%]", normalDirection, normalDirection * radius_);
    return wall.getDistanceNormalOverlapSuperquadric(p, distance, normal_return, overlap);
}

/*!
 * \param[in] is The input stream from which the AxisymmetricIntersectionOfWalls
 * is read, usually a restart file.
 */
void SimpleDrumSuperquadrics::read(std::istream& is)
{
    BaseWall::read(is);
    wall.read(is);
    is >> radius_;
}

/*!
 * \param[in] os The output stream where the AxisymmetricIntersectionOfWalls must be written
 *  to, usually a restart file.
 */
void SimpleDrumSuperquadrics::write(std::ostream& os) const
{
    BaseWall::write(os);
    wall.write(os);
    os << radius_;
}

/*!
 * \return The string "AxisymmetricIntersectionOfWalls".
 */
std::string SimpleDrumSuperquadrics::getName() const
{
    return "SimpleDrumSuperquadrics";
}

void SimpleDrumSuperquadrics::setAxis(Vec3D a)
{
    setOrientationViaNormal(a);
}
