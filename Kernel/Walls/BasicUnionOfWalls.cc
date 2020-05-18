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

#include "BasicUnionOfWalls.h"
#include "Particles/BaseParticle.h"
#include "Particles/SphericalParticle.h"
#include "InteractionHandler.h"
#include "WallHandler.h"
#include "DPMBase.h"

BasicUnionOfWalls::BasicUnionOfWalls()
{
    logger(DEBUG, "BasicUnionOfWalls::BasicUnionOfWalls ) finished");
}

/*!
 * \param[in] w BasicUnionOfWalls that has to be copied.
 * \details First copy the attributes of the BaseWall, then copy the ones that are
 * specific for the BasicUnionOfWalls.
 */
BasicUnionOfWalls::BasicUnionOfWalls(const BasicUnionOfWalls& b)
        : BaseWall(b)
{
    for (auto& w : b.walls_)
    {
        walls_.push_back(w->copy());
    }
    logger(DEBUG, "BasicUnionOfWalls::BasicUnionOfWalls(const BasicUnionOfWalls &p) finished");
}

BasicUnionOfWalls::~BasicUnionOfWalls()
{
    for (auto& w : walls_)
    {
        delete w;
    }
    logger(DEBUG, "BasicUnionOfWalls::~BasicUnionOfWalls finished");
}

/*!
 * Wall copy method. It calls the copy constructor of this Wall, useful for polymorphism
 */
BasicUnionOfWalls* BasicUnionOfWalls::copy() const
{
    return new BasicUnionOfWalls(*this);
}

/*!
 * \detail Suppose your simulation adds to an BasicUnionOfWalls after a
 * certain time or condition is met. Checking the number of objects is useful
 * for checking if this has happened yet, when restarting.
 */
unsigned long BasicUnionOfWalls::getNumberOfObjects()
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
void BasicUnionOfWalls::add(BaseWall& wall)
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
bool BasicUnionOfWalls::getDistanceAndNormal(const BaseParticle& p, Mdouble& distance, Vec3D& normal) const
{
    if (walls_.empty()) return false;

    Vec3D position = p.getPosition() - getPosition();
    getOrientation().rotateBack(position);
    SphericalParticle shifted;
    shifted.setSpecies(p.getSpecies());
    shifted.setPosition(position);
    shifted.setRadius(p.getRadius());

    //check wall after wall; the first wall that returns an interaction is chosen
    for (auto w : walls_) {
        if (w->getDistanceAndNormal(shifted, distance, normal) == true) {
            getOrientation().rotate(normal);
            return true;
        }
    }

    return false;
}

/*!
 * \param[in] is The input stream from which the BasicUnionOfWalls is read.
 */
void BasicUnionOfWalls::read(std::istream& is)
{
    BaseWall::read(is);
    std::string dummy;
    unsigned size;
    is >> dummy >> size;
    for (unsigned i = 0; i < size; i++)
    {
        std::string type;
        is >> type;
        BaseWall* wall = getHandler()->createObject(type);
        wall->setHandler(getHandler());
        wall->read(is);
        walls_.push_back(wall);
        walls_.back()->setId(walls_.size());
    }
}

/*!
 * \param[in] os The output stream the BasicUnionOfWalls is written to.
 */
void BasicUnionOfWalls::write(std::ostream& os) const
{
    BaseWall::write(os);
    os << " numWalls " << walls_.size();
    for (auto w : walls_)
    {
        os << " ";
        w->write(os);
    }
}

/*!
 * \return The string "BasicUnionOfWalls", which is the name of this class.
 */
std::string BasicUnionOfWalls::getName() const
{
    return "BasicUnionOfWalls";
}

void BasicUnionOfWalls::getVTK(std::vector<Vec3D>& points, std::vector<std::vector<double>>& triangleStrips)
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
    //add position of the BasicUnionOfWalls to the point
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

BaseWall* BasicUnionOfWalls::getObject(unsigned i)
{
    logger.assert_always(walls_.size() > i, "Index % exceeds number of walls %", i, walls_.size());
    return walls_[i];
}
