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

#include "BasicIntersectionOfWalls.h"
#include "Particles/BaseParticle.h"
#include "Particles/SphericalParticle.h"
#include "InteractionHandler.h"
#include "WallHandler.h"
#include "DPMBase.h"

BasicIntersectionOfWalls::BasicIntersectionOfWalls()
{
    logger(DEBUG, "BasicIntersectionOfWalls::BasicIntersectionOfWalls ) finished");
}

/*!
 * \param[in] w BasicIntersectionOfWalls that has to be copied.
 * \details First copy the attributes of the BaseWall, then copy the ones that are
 * specific for the BasicIntersectionOfWalls.
 */
BasicIntersectionOfWalls::BasicIntersectionOfWalls(const BasicIntersectionOfWalls& b)
        : BaseWall(b)
{
    for (auto& w : b.walls_)
    {
        walls_.push_back(w->copy());
    }
    logger(DEBUG, "BasicIntersectionOfWalls::BasicIntersectionOfWalls(const BasicIntersectionOfWalls &p) finished");
}

BasicIntersectionOfWalls::~BasicIntersectionOfWalls()
{
    for (auto& w : walls_)
    {
        delete w;
    }
    logger(DEBUG, "BasicIntersectionOfWalls::~BasicIntersectionOfWalls finished");
}

/*!
 * Wall copy method. It calls the copy constructor of this Wall, useful for polymorphism
 */
BasicIntersectionOfWalls* BasicIntersectionOfWalls::copy() const
{
    return new BasicIntersectionOfWalls(*this);
}

/*! 
 * \detail Suppose your simulation adds to an BasicIntersectionOfWalls after a
 * certain time or condition is met. Checking the number of objects is useful
 * for checking if this has happened yet, when restarting. 
 */
unsigned long BasicIntersectionOfWalls::getNumberOfObjects()
{
    return walls_.size();
}


/*
 * \param[in] normal A Vec3D that represents the normal to the wall.
 * \param[in] point A Vec3D which is a point on the wall.
 * \details Sets the wall such that for all points x on the wall it holds that 
 * normal*x=normal*point.
 */
///\todo TW maybe the Restricted wall should be templated with the wall type such that we don't need to use new and delete.
void BasicIntersectionOfWalls::add(BaseWall& wall)
{
    walls_.push_back(wall.copy());
    walls_.back()->setId(walls_.size());
}

/*!
 * \param[in] p BaseParticle for which the distance to the wall must be computed.
 * \param[out] distance Distance between the particle and the wall.
 * \param[out] normal The normal of this wall, will only be set if there is a collision.
 * \return A boolean value for whether or not there is a collision.
 * \details Return the smallest distance and normal
 */
bool BasicIntersectionOfWalls::getDistanceAndNormal(const BaseParticle& p, Mdouble& distance, Vec3D& normal) const
{
    if (walls_.empty())
    {
        logger(DEBUG, "Empty BasicIntersectionOfWalls");
        return false;
    }
    
    distance = -constants::inf; //distance of the closest wall
    Mdouble distance2 = -constants::inf; //distance of the second-closest wall
    Mdouble distance3 = -constants::inf; //distance of the third-closest wall
    Mdouble distanceCurrent; //distance of teh current wall
    Vec3D normal2; //normal of the second-closest wall
    Vec3D normal3; //normal of the third-closest wall
    Vec3D normalCurrent; //normal of the current wall
    unsigned int id = static_cast<unsigned int>(-1); //id of the closest wall
    unsigned int id2 = static_cast<unsigned int>(-1); //id of the second-closest wall
    unsigned int id3 = static_cast<unsigned int>(-1); //id of the third-closest wall
    Mdouble wallInteractionRadius = p.getWallInteractionRadius(this);
    Vec3D position = p.getPosition() - getPosition();
    getOrientation().rotateBack(position);
    SphericalParticle shifted;
    shifted.setSpecies(p.getSpecies());
    shifted.setPosition(position);
    shifted.setRadius(p.getRadius());
    
    //The object has to touch each wall each wall (distanceCurrent) and keep the minimum distance (distance) and wall index (id)
    for (unsigned int i = 0; i < walls_.size(); i++)
    {
        //immediately stop if one of teh walls is not in interaction distance
        if (!walls_[i]->getDistanceAndNormal(shifted, distanceCurrent, normalCurrent))
            return false;
        // Find out which of the walls is interacting with the particle.
        // Keep the minimum distance (distance) and wall index (id)
        // and store up to two walls (id2, id3) and their distances (distance2, distance3),
        // if the possible contact point is near the intersection between id and id2 (and id3)
        if (distanceCurrent > distance)
        {
            if (distance > -wallInteractionRadius)
            {
                if (distance2 > -wallInteractionRadius)
                {
                    distance3 = distance;
                    normal3 = normal;
                    id3 = id;
                }
                else
                {
                    distance2 = distance;
                    normal2 = normal;
                    id2 = id;
                }
            }
            distance = distanceCurrent;
            normal = normalCurrent;
            id = i;
        }
        else if (distanceCurrent > -wallInteractionRadius)
        {
            if (distance2 > -wallInteractionRadius)
            {
                distance3 = distanceCurrent;
                normal3 = normalCurrent;
                id3 = i;
            }
            else
            {
                distance2 = distanceCurrent;
                normal2 = normalCurrent;
                id2 = i;
            }
        }
    }
    //logger(INFO, "particle %: contact with % % %, n % % %", p.getId(), distance, distance2, distance3, normal, normal2, normal3);
    
    
    //If we are here, the closest wall is id;
    //if distance2>-P.Radius, it's a possible face or vertex contact
    //if distance3>-P.Radius, it's a possible vertex contact
    if (distance2 > -wallInteractionRadius)
    {
        //contact point on wall id
        SphericalParticle contact;
        contact.setSpecies(p.getSpecies());
        contact.setPosition(position + distance * normal);
        contact.setRadius(0);
        
        //If the distance of D to id2 is positive, the contact is with the intersection
        bool contactPointOutsideID2 = !walls_[id2]->getDistanceAndNormal(contact, distanceCurrent, normalCurrent);
        
        if (distance3 > -wallInteractionRadius &&
            !walls_[id3]->getDistanceAndNormal(contact, distanceCurrent, normalCurrent))
        {
            if (contactPointOutsideID2)
            {
                //possible contact is with intersection of id,id2,id3
                //we know id2<id3
                unsigned int index =
                        (id < id2) ? ((id3 - 2) * (id3 - 1) * id3 / 6 + (id2 - 1) * id2 / 2 + id) :
                        (id < id3) ? ((id3 - 2) * (id3 - 1) * id3 / 6 + (id - 1) * id / 2 + id2) :
                        ((id - 2) * (id - 1) * id / 6 + (id3 - 1) * id3 / 2 + id2);
                //find vertex C
                Matrix3D N(normal.X, normal.Y, normal.Z, normal2.X, normal2.Y, normal2.Z, normal3.X, normal3.Y,
                           normal3.Z);
                Vec3D P(Vec3D::dot(position, normal) + distance, Vec3D::dot(position, normal2) + distance2,
                        Vec3D::dot(position, normal3) + distance3);
                Vec3D vertex = N.ldivide(P);
                normal = position - vertex;
                distance = sqrt(normal.getLengthSquared());
                if (distance >= wallInteractionRadius)
                    return false; //no contact
                normal /= -distance;
                //logger(INFO, "vertex contact with % % %: pos %, n %, d %", id, id2, id3, p.getPosition(), normal, distance);
                getOrientation().rotate(normal);
                return true;
            }
            else
            {
                contactPointOutsideID2 = true;
                distance2 = distance3;
                id2 = id3;
            }
        }
        
        if (contactPointOutsideID2)
        {
            //find contact point C
            normal3 = Vec3D::cross(normal, normal2);
            Matrix3D N(normal.X, normal.Y, normal.Z, normal2.X, normal2.Y, normal2.Z, normal3.X, normal3.Y,
                       normal3.Z);
            Vec3D P(Vec3D::dot(position, normal) + distance, Vec3D::dot(position, normal2) + distance2,
                    Vec3D::dot(position, normal3));
            Vec3D C = N.ldivide(P);
            normal = position - C;
            distance = sqrt(normal.getLengthSquared());
            if (distance >= wallInteractionRadius)
                return false; //no contact
            normal /= -distance;
            //logger(INFO, "edge contact with % %: normal %", id, id2, normal);
            getOrientation().rotate(normal);
            return true;
        }
    }
    //logger(INFO, "face contact with %", id);
    getOrientation().rotate(normal);
    return true;
}

/*!
 * \param[in] is The input stream from which the BasicIntersectionOfWalls is read.
 */
void BasicIntersectionOfWalls::read(std::istream& is)
{
    BaseWall::read(is);
    std::string dummy;
    unsigned size;
    // read e.g. "numIntersectionOfWalls 2"
    is >> dummy >> size;
    for (unsigned i = 0; i < size; i++)
    {
        std::string type;
        // read e.g. "IntersectionOfWalls"
        is >> type;
        BaseWall* wall = getHandler()->createObject(type);
        wall->setHandler(getHandler());
        wall->read(is);
        walls_.push_back(wall);
        walls_.back()->setId(walls_.size());
    }
}

/*!
 * \param[in] os The output stream the BasicIntersectionOfWalls is written to.
 */
void BasicIntersectionOfWalls::write(std::ostream& os) const
{
    BaseWall::write(os);
    os << " numIntersectionOfWalls " << walls_.size();
    for (auto w : walls_)
    {
        os << " ";
        w->write(os);
    }
}

/*!
 * \return The string "BasicIntersectionOfWalls", which is the name of this class.
 */
std::string BasicIntersectionOfWalls::getName() const
{
    return "BasicIntersectionOfWalls";
}

void BasicIntersectionOfWalls::getVTK(std::vector<Vec3D>& points, std::vector<std::vector<double>>& triangleStrips)
{
//    for (auto w : walls_) {
//        w->writeVTK (points, triangleStrips);
//    }
    //writes points and strips for all walls; points are added to the global point vector, but the strips are held back
    std::vector<std::vector<double>> myTriangleStrips;
    unsigned long n = points.size();
    for (auto w : walls_)
    {
        w->getVTK(points, myTriangleStrips);
    }
    //add position of the BasicIntersectionOfWalls to the point
    for (std::vector<Vec3D>::iterator p = points.begin() + n; p != points.end(); p++)
    {
        getOrientation().rotate(*p);
    }
    //create a vector which points are in the wall (actually, only the new points are necessary)
    std::vector<bool> pointInWall;
    pointInWall.reserve(points.size());
    SphericalParticle particle;
    particle.setSpecies(getSpecies());
    particle.setRadius(1e-10); //points within that distance are declared part of the wall
    Mdouble distance;
    Vec3D normal;
    for (auto p : points)
    {
        particle.setPosition(p);
        pointInWall.push_back(getDistanceAndNormal(particle, distance, normal));
    }
    //now loop through myTriangleStrips to find the strip parts that are fully inside the wall
    std::vector<double> strip;
    for (auto t : myTriangleStrips)
    {
        for (unsigned i : t)
        {
            if (pointInWall[i] == true)
            {
                strip.push_back(i);
            }
            else
            {
                if (strip.size() > 2)
                    triangleStrips.push_back(strip);
                strip.clear();
            }
        }
        if (strip.size() > 2)
            triangleStrips.push_back(strip);
        strip.clear();
    }
    ///\todo this function could be improved; might not plot full wall
}

BaseWall* BasicIntersectionOfWalls::getObject(unsigned i)
{
    logger.assert_always(walls_.size() > i, "Index % exceeds number of walls %", i, walls_.size());
    return walls_[i];
}
