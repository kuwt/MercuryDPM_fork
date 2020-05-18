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

#include "InfiniteWall.h"
#include "Particles/BaseParticle.h"
#include "InteractionHandler.h"
#include "WallHandler.h"
#include "DPMBase.h"

InfiniteWall::InfiniteWall()
{
    logger(DEBUG, "InfiniteWall::InfiniteWall ) finished");
}

/*!
 * \param[in] w InfiniteWall that has to be copied.
 * \details First copy the attributes of the BaseWall, then copy the ones that are
 * specific for the InfiniteWall.
 */
InfiniteWall::InfiniteWall(const InfiniteWall& w)
        : BaseWall(w)
{
    logger(DEBUG, "InfiniteWall::InfiniteWall(const InfiniteWall &p) finished");
}

InfiniteWall::InfiniteWall(const ParticleSpecies* s)
        : BaseWall()
{
    setSpecies(s);
    logger(DEBUG, "InfiniteWall::InfiniteWall(const ParticleSpecies* s) finished");
}

InfiniteWall::InfiniteWall(Vec3D normal, Vec3D point, const ParticleSpecies* species)
{
    setNormal(normal);
    setPosition(point);
    setSpecies(species);
}

/*!
 *
 * @param PointA first coordinate
 * @param PointB second coordinate
 * @param PointC third coordinate
 * @param species
 * \details Builds an infinite wall through 3 points, normal is defined with Right hand rule following the three input points.
 */
InfiniteWall::InfiniteWall(Vec3D PointA, Vec3D PointB, Vec3D PointC, const ParticleSpecies* species)
{
//function returns infinite wall passing through 3 points
//subtract second and third from the first
    Vec3D SubtB = PointA - PointB;
    Vec3D SubtC = PointA - PointC;
    
    Vec3D WallNormal = Vec3D::cross(SubtB, SubtC);
    //Check if walls coordinates  inline, if true Do not build wall and give error message.
    
    if (WallNormal.getLengthSquared() == 0.0)
    {
        logger(ERROR,
               "Error Building InfiniteWall out of 3 coordinates. Coordinates are in line, Wall not constructed.");
    }
    else
    {
        setNormal(WallNormal);
        setPosition(PointA);
        setSpecies(species);
    }
    
    
}

InfiniteWall::~InfiniteWall()
{
    logger(DEBUG, "InfiniteWall::~InfiniteWall finished");
}

/*!
 * Wall copy method. It calls the copy constructor of this Wall, useful for polymorphism
 */
InfiniteWall* InfiniteWall::copy() const
{
    return new InfiniteWall(*this);
}

/*
 * \details Defines a standard wall, given a normal vector pointing into the
 * wall (i.e. out of the flow domain) and a point that the wall goes through,
 * to give a plane defined by
 * normal*x=normal*point.
 * \param[in] normal A Vec3D that represents the normal to the wall.
 * \param[in] point A Vec3D which is a point on the wall.
 */
void InfiniteWall::set(Vec3D normal, Vec3D point)
{
    setNormal(normal);
    setPosition(point);
}

/*!
 * \param[in] normal The vector normal to the wall.
 */
void InfiniteWall::setNormal(const Vec3D normal)
{
    setOrientationViaNormal(normal);
}

/*!
 * \details Defines a standard wall, given an normal vector pointing into the wall (i.e. out of the flow domain), 
 * to give a plane defined by normal*x=position
 * \param[in] normal A Vec3D that represents the normal vector to the wall.
 * \param[in] positionInNormalDirection The position of the wall in the direction
 *  of the normal vector.
 */
void InfiniteWall::set(Vec3D normal, Mdouble positionInNormalDirection)
{
    logger(WARN, "InfiniteWall::set(Vec3D, Mdouble) is deprecated. Use set(Vec3D, Vec3D) instead.");
    set(normal, positionInNormalDirection * normal);
}

/*!
 * \param[in] otherPosition The position to which the distance must be computed to.
 * \return The distance of the wall to the particle.
 */
Mdouble InfiniteWall::getDistance(Vec3D otherPosition) const
{
    return getOrientation().getDistance(otherPosition, getPosition());
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
bool InfiniteWall::getDistanceAndNormal(const BaseParticle& p, Mdouble& distance, Vec3D& normal_return) const
{
    distance = getDistance(p.getPosition());
    if (distance >= p.getWallInteractionRadius(this))
        return false;
    normal_return = getOrientation().getAxis();
    return true;
}

/*!
 * \param[in] is The input stream from which the InfiniteWall is read. Only needed for backward compatibility.
 */
void InfiniteWall::read(std::istream& is)
{
    BaseWall::read(is);
    Vec3D normal;
    if (helpers::readOptionalVariable(is,"normal",normal)) setNormal(normal);
}

/*!
 * \param[in] is The input stream from which the InfiniteWall old style is read.
 */
void InfiniteWall::oldRead(std::istream& is)
{
    std::string dummy;
    Vec3D velocity;
    Vec3D position;
    Vec3D normal;
    is >> dummy >> normal >> dummy >> position >> dummy >> velocity;
    setPosition(position);
    setVelocity(velocity);
    setOrientation(Quaternion(normal));
}

/*!
 * \return The string "InfiniteWall", which is the name of this class.
 */
std::string InfiniteWall::getName() const
{
    return "InfiniteWall";
}

/*!
 * \return The 3D vector that represents the normal to the wall.
 */
Vec3D InfiniteWall::getNormal() const
{
    return getOrientation().getAxis();
}

/*!
 * Used to create an array of points for representing the wall in the VTK viewer.
 * E.g., for a InfiniteWall through p=(0,0,0) with normal n=(0,0,-1) in a domain (-1,1)^3, createVTK returns
 *   {(0,-1,-1), (0,-1,1), (0,1,1), (0,1,-1)}
 * Calling addToVTK will then create a triangle strip connecting these points with triangle faces.
 */
void InfiniteWall::createVTK(std::vector<Vec3D>& myPoints) const
{
    Vec3D max = getHandler()->getDPMBase()->getMax();
    Vec3D min = getHandler()->getDPMBase()->getMin();
    createVTK(myPoints, min, max);
}

void InfiniteWall::createVTK(std::vector<Vec3D>& myPoints, const Vec3D min, const Vec3D max) const
{
    const Vec3D& n = getOrientation().getAxis();
    const Vec3D& p = getPosition();
    
    if (fabs(n.X) > 0.5)
    {
        // If the wall normal has a nonzero x-component,
        // We first find four intersection points with the four domain edges pointing in x-direction.
        // Because these points might be outside the domain in x-direction, we use the intersection function have points in the domain
        myPoints.emplace_back(p.X - ((min.Y - p.Y) * n.Y + (min.Z - p.Z) * n.Z) / n.X, min.Y, min.Z);
        myPoints.emplace_back(p.X - ((min.Y - p.Y) * n.Y + (max.Z - p.Z) * n.Z) / n.X, min.Y, max.Z);
        myPoints.emplace_back(p.X - ((max.Y - p.Y) * n.Y + (max.Z - p.Z) * n.Z) / n.X, max.Y, max.Z);
        myPoints.emplace_back(p.X - ((max.Y - p.Y) * n.Y + (min.Z - p.Z) * n.Z) / n.X, max.Y, min.Z);
        intersectVTK(myPoints, Vec3D(+1, 0, 0), Vec3D(max.X, 0, 0));
        intersectVTK(myPoints, Vec3D(-1, 0, 0), Vec3D(min.X, 0, 0));
    }
    else if (fabs(n.Y) > 0.5)
    {
        // Else, if the wall normal has a nonzero y-component ...
        myPoints.emplace_back(min.X, p.Y - ((min.X - p.X) * n.X + (min.Z - p.Z) * n.Z) / n.Y, min.Z);
        myPoints.emplace_back(min.X, p.Y - ((min.X - p.X) * n.X + (max.Z - p.Z) * n.Z) / n.Y, max.Z);
        myPoints.emplace_back(max.X, p.Y - ((max.X - p.X) * n.X + (max.Z - p.Z) * n.Z) / n.Y, max.Z);
        myPoints.emplace_back(max.X, p.Y - ((max.X - p.X) * n.X + (min.Z - p.Z) * n.Z) / n.Y, min.Z);
        intersectVTK(myPoints, Vec3D(0, +1, 0), Vec3D(0, max.Y, 0));
        intersectVTK(myPoints, Vec3D(0, -1, 0), Vec3D(0, min.Y, 0));
    }
    else
    {
        // Else, the wall normal has to have a a nonzero z-component ...
        myPoints.emplace_back(min.X, min.Y, p.Z - ((min.Y - p.Y) * n.Y + (min.X - p.X) * n.X) / n.Z);
        myPoints.emplace_back(min.X, max.Y, p.Z - ((max.Y - p.Y) * n.Y + (min.X - p.X) * n.X) / n.Z);
        myPoints.emplace_back(max.X, max.Y, p.Z - ((max.Y - p.Y) * n.Y + (max.X - p.X) * n.X) / n.Z);
        myPoints.emplace_back(max.X, min.Y, p.Z - ((min.Y - p.Y) * n.Y + (max.X - p.X) * n.X) / n.Z);
        intersectVTK(myPoints, Vec3D(0, 0, +1), Vec3D(0, 0, max.Z));
        intersectVTK(myPoints, Vec3D(0, 0, -1), Vec3D(0, 0, min.Z));
    }
}

void InfiniteWall::writeVTK(VTKContainer& vtk) const
{
    std::vector<Vec3D> points;
    createVTK(points);
    addToVTK(points, vtk);
}


bool
InfiniteWall::getDistanceNormalOverlapSuperquadric(const SuperQuadricParticle& p, Mdouble& distance, Vec3D& normal_return,
                                                   Mdouble& overlap) const
{
    //first check: if the bounding sphere does not touch the wall, there is no contact.
    if (getDistance(p.getPosition()) >= p.getWallInteractionRadius(this))
    {
        return false;
    }
    Vec3D normalBodyFixed = getOrientation().getAxis();
    p.getOrientation().rotateBack(normalBodyFixed);
    Vec3D xWallBodyFixed = getPosition() - p.getPosition();
    p.getOrientation().rotateBack(xWallBodyFixed);
    Vec3D axes = p.getAxes();
    Mdouble eps1 = p.getExponentEps1();
    Mdouble eps2 = p.getExponentEps2();
    
    Vec3D furthestPoint = getFurthestPointSuperQuadric(normalBodyFixed, axes, eps1, eps2);
    overlap = Vec3D::dot(xWallBodyFixed - furthestPoint, -normalBodyFixed);
    if (overlap > 0)
    {
        Vec3D overlapBody = overlap * normalBodyFixed;
        Vec3D contactPoint = furthestPoint - overlapBody / 2;
        p.getOrientation().rotate(contactPoint);
        contactPoint += p.getPosition();
        distance = (contactPoint - overlapBody / 2 - p.getPosition()).getLength();
        normal_return = getOrientation().getAxis();
        return true;
    }
    return false;
}


///Largely untested, use at your own risk for anything other than ellipsoids.
Vec3D InfiniteWall::getFurthestPointSuperQuadric(const Vec3D& normalBodyFixed, const Vec3D& axes, Mdouble eps1,
                                                 Mdouble eps2) const
{
    Vec3D furthestPoint;
    if (std::abs(normalBodyFixed.X) > 1e-10)
    {
        Mdouble alpha = std::abs(normalBodyFixed.Y * axes.Y / normalBodyFixed.X / axes.X);
        Mdouble gamma = std::pow(1 + std::pow(alpha, 2 / eps2), eps2 / eps1 - 1);
        Mdouble beta = std::pow(gamma * std::abs(normalBodyFixed.Z * axes.Z / normalBodyFixed.X / axes.X),
                                eps1 / (2 - eps1));
        furthestPoint.X = axes.X * mathsFunc::sign(normalBodyFixed.X) /
                          (std::pow(std::pow(1 + std::pow(alpha, 2 / eps2), eps2 / eps1) + std::pow(beta, 2 / eps1),
                                    eps1 / 2));
        furthestPoint.Y = axes.Y * alpha / axes.X * std::abs(furthestPoint.X) * mathsFunc::sign(normalBodyFixed.Y);
        furthestPoint.Z = axes.Z * beta / axes.X * std::abs(furthestPoint.X) * mathsFunc::sign(normalBodyFixed.Z);
    }
    else if (std::abs(normalBodyFixed.Y) > 1e-10)
    {
        Mdouble beta = std::pow(std::abs(normalBodyFixed.Z * axes.Z / normalBodyFixed.Y / axes.Y), eps1 / (2 - eps1));
        furthestPoint.Y = axes.Y / std::pow((1 + std::pow(beta, 2/eps1)), eps1 / 2) * mathsFunc::sign(normalBodyFixed.Y);
        furthestPoint.Z = axes.Z / axes.Y * std::abs(furthestPoint.Y) * beta * mathsFunc::sign(normalBodyFixed.Z);
    }
    else
    {
        furthestPoint.Z = axes.Z * mathsFunc::sign(normalBodyFixed.Z);
    }
    return furthestPoint;
}
