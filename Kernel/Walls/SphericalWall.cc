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

#include <limits>
#include "Logger.h"
#include "SphericalWall.h"
#include "Particles/BaseParticle.h"
#include "InteractionHandler.h"
#include "WallHandler.h"
#include "DPMBase.h"

SphericalWall::SphericalWall()
{
    radius_ = std::numeric_limits<double>::quiet_NaN();
    logger(DEBUG, "SphericalWall::SphericalWall ) finished");
}

/*!
 * \param[in] w SphericalWall that has to be copied.
 * \details First copy the attributes of the BaseWall, then copy the ones that are
 * specific for the SphericalWall.
 */
SphericalWall::SphericalWall(const SphericalWall& w)
        : BaseWall(w)
{
    radius_ = w.radius_;
    logger(DEBUG, "SphericalWall::SphericalWall(const SphericalWall &p) finished");
}

SphericalWall::SphericalWall(Mdouble radius, const ParticleSpecies* species)
{
    setRadius(radius);
}

SphericalWall::~SphericalWall()
{
    logger(DEBUG, "SphericalWall::~SphericalWall finished");
}

/*!
 * Wall copy method. It calls the copy constructor of this Wall, useful for polymorphism
 */
SphericalWall* SphericalWall::copy() const
{
    return new SphericalWall(*this);
}

/*
 * \param[in] normal A Vec3D that represents the normal to the wall.
 * \param[in] point A Vec3D which is a point on the wall.
 * \details Sets the wall such that for all points x on the wall it holds that 
 * normal*x=normal*point.
 */
void SphericalWall::setRadius(Mdouble radius)
{
    logger.assert(radius >= 0, "radius=% cannot be negative", radius);
    radius_ = radius;
}

/*!
 * \param[in] otherPosition The position to which the distance must be computed to.
 * \return The distance of the wall to the particle.
 */
Mdouble SphericalWall::getDistance(const Vec3D& otherPosition) const
{
    return Vec3D::getLength(getPosition() - otherPosition);
}

/*!
 * \return The radius of the wall .
 */
Mdouble SphericalWall::getRadius() const
{
    return radius_;
}

/*!
 * \param[in] p BaseParticle for which the distance to the wall must be computed.
 * \param[out] distance Distance between the particle and the wall.
 * \param[out] normal_return The normal of this wall, will only be set if there is a collision.
 * \return A boolean value for whether or not there is a collision.
 * \details First the distance is checked. If there is no collision, this
 * function will return false and give the distance. If there is a collision, the
 * function will return true and give the distance and the normal vector of this wall.
 * Since this function should be called before calculating any 
 * Particle-Wall interactions, it can also be used to set the normal vector in 
 * case of curved walls.
 */
bool SphericalWall::getDistanceAndNormal(const BaseParticle& p, Mdouble& distance, Vec3D& normal_return) const
{
    normal_return = p.getPosition() - getPosition();
    distance = Vec3D::getLength(normal_return) - radius_;
    if (distance >= p.getWallInteractionRadius(this))
        return false;
    normal_return /= distance + radius_;
    //logger(WARN,"p% q% q% q% q%", getPosition(), p.getPosition(), normal_return, distance, p.getWallInteractionRadius(this));
    return true;
}

/*!
 * \param[in] is The input stream from which the SphericalWall is read.
 */
void SphericalWall::read(std::istream& is)
{
    BaseWall::read(is);
    std::string dummy;
    is >> dummy >> radius_;
}

/*!
 * \param[in] os The output stream the SphericalWall is written to.
 */
void SphericalWall::write(std::ostream& os) const
{
    BaseWall::write(os);
    os << " radius " << radius_;
}

/*!
 * \return The string "SphericalWall", which is the name of this class.
 */
std::string SphericalWall::getName() const
{
    return "SphericalWall";
}
