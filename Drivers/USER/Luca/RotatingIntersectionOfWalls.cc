//Copyright (c) 2013-2014, The MercuryDPM Developers Team. All rights reserved.
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

#include "RotatingIntersectionOfWalls.h"
#include <InteractionHandler.h>
#include <Particles/BaseParticle.h>

RotatingIntersectionOfWalls::RotatingIntersectionOfWalls()
{
    logger(DEBUG, "RotatingIntersectionOfWalls() constructed.");
}

/*!
 * \param[in] other The RotatingIntersectionOfWalls that must be copied.
 */
RotatingIntersectionOfWalls::RotatingIntersectionOfWalls(const RotatingIntersectionOfWalls& other)
    : BaseWall(other)
{
    wallObjects_ = other.wallObjects_;
    for (auto& wall : wallObjects_)
    {
        if (getHandler() != nullptr)
            wall.setHandler(getHandler());
    }
    A_ = other.A_;
    AB_ = other.AB_;
    C_ = other.C_;
    logger(DEBUG, "RotatingIntersectionOfWalls(RotatingIntersectionOfWalls&) constructed.");
}

RotatingIntersectionOfWalls::RotatingIntersectionOfWalls(std::vector<RotatingIntersectionOfWalls::normalAndPosition> walls, const ParticleSpecies* species)
{
    setSpecies(species);
    for (auto wall : walls)
    {
        addObject(wall.normal,wall.position);
    }
}


RotatingIntersectionOfWalls::~RotatingIntersectionOfWalls()
{
    logger(DEBUG, "~RotatingIntersectionOfWalls() has been called.");
}

void RotatingIntersectionOfWalls::setSpecies(const ParticleSpecies* species)
{
    BaseWall::setSpecies(species);
    for (auto wall : wallObjects_)
    {
        wall.setSpecies(species);
    }
}

/*!
 * \param[in] other The RotatingIntersectionOfWalls that must be copied.
 */
RotatingIntersectionOfWalls& RotatingIntersectionOfWalls::operator=(const RotatingIntersectionOfWalls& other)
{
    logger(DEBUG, "RotatingIntersectionOfWalls::operator= called.");
    if (this == &other)
    {
        return *this;
    }
    return *(other.copy());
}

/*!
 * \return pointer to a RotatingIntersectionOfWalls object allocated using new.
 */
RotatingIntersectionOfWalls* RotatingIntersectionOfWalls::copy() const
{
    return new RotatingIntersectionOfWalls(*this);
}

void RotatingIntersectionOfWalls::setHandler(WallHandler* wallHandler)
{
    BaseWall::setHandler(wallHandler);
    for (InfiniteWall w : wallObjects_)
    {
        BaseWall::setHandler(wallHandler);
    }
}

/*!
 * \param[in] normal The normal to this wallObject.
 * \param[in] point One of the points of the wallObject.
 * \details Adds a wall to the set of finite walls, given an outward unit normal
 *  vector s.t. normal*x=normal*point for all x of the wallObject. First make the
 * InfiniteWall, then compute all intersections, which are then stored in A_, AB_
 * and C_.
 */
void RotatingIntersectionOfWalls::addObject(Vec3D normal, Vec3D point)
{
    normal.normalise();

    //n is the index of the new wall
    std::size_t n = wallObjects_.size();
    InfiniteWall w;
    if (getSpecies()) w.setSpecies(getSpecies());
    w.set(normal, point);
    wallObjects_.push_back(w);

    // AB[n*(n-1)/2+m] is the direction of the intersecting line between walls m and n, m<n
    // A[n*(n-1)/2+m] is a point on the intersecting line between walls m and n, m<n
    // See http://www.netcomuk.co.uk/~jenolive/vect18d.html for finding the line where two planes meet
    AB_.resize(n * (n + 1) / 2); //0 + 1 + 2 + ... + indexNew, total number of walls you need
    A_.resize(n * (n + 1) / 2);
    for (std::size_t m = 0; m < n; m++)
    {
        std::size_t id = (n - 1) * n / 2 + m;
        //first we cross the wall normals and normalize to obtain AB
        AB_[id] = Vec3D::cross(wallObjects_[m].getNormal(), wallObjects_[n].getNormal());
        AB_[id] /= sqrt(AB_[id].getLengthSquared());
        //then we find a point A (using AB*x=0 as a third plane)
        Mdouble invdet = 1.0 / (+wallObjects_[n].getNormal().X * (wallObjects_[m].getNormal().Y * AB_[id].Z - AB_[id].Y * wallObjects_[m].getNormal().Z)
                - wallObjects_[n].getNormal().Y * (wallObjects_[m].getNormal().X * AB_[id].Z - wallObjects_[m].getNormal().Z * AB_[id].X)
                + wallObjects_[n].getNormal().Z * (wallObjects_[m].getNormal().X * AB_[id].Y - wallObjects_[m].getNormal().Y * AB_[id].X));

        A_[id] = Vec3D( +(wallObjects_[m].getNormal().Y * AB_[id].Z - AB_[id].Y * wallObjects_[m].getNormal().Z) * Vec3D::dot(wallObjects_[n].getPosition(),wallObjects_[n].getNormal())
                        -(wallObjects_[n].getNormal().Y * AB_[id].Z - wallObjects_[n].getNormal().Z * AB_[id].Y) * Vec3D::dot(wallObjects_[m].getPosition(),wallObjects_[m].getNormal())
                        +(wallObjects_[n].getNormal().Y * wallObjects_[m].getNormal().Z - wallObjects_[n].getNormal().Z * wallObjects_[m].getNormal().Y) * 0.0,
                        -(wallObjects_[m].getNormal().X * AB_[id].Z - wallObjects_[m].getNormal().Z * AB_[id].X) * Vec3D::dot(wallObjects_[n].getPosition(),wallObjects_[n].getNormal())
                        +(wallObjects_[n].getNormal().X * AB_[id].Z - wallObjects_[n].getNormal().Z * AB_[id].X) * Vec3D::dot(wallObjects_[m].getPosition(),wallObjects_[m].getNormal())
                        -(wallObjects_[n].getNormal().X * wallObjects_[m].getNormal().Z - wallObjects_[m].getNormal().X * wallObjects_[n].getNormal().Z) * 0.0,
                        +(wallObjects_[m].getNormal().X * AB_[id].Y - AB_[id].X * wallObjects_[m].getNormal().Y) * Vec3D::dot(wallObjects_[n].getPosition(),wallObjects_[n].getNormal())
                        -(wallObjects_[n].getNormal().X * AB_[id].Y - AB_[id].X * wallObjects_[n].getNormal().Y) * Vec3D::dot(wallObjects_[m].getPosition(),wallObjects_[m].getNormal())
                        +(wallObjects_[n].getNormal().X * wallObjects_[m].getNormal().Y - wallObjects_[m].getNormal().X * wallObjects_[n].getNormal().Y) * 0.0) * invdet;
    }

    // C[(n-2)*(n-1)*n/6+(m-1)*m/2+l] is a point intersecting walls l, m and n, l<m<n
    C_.resize((n - 1) * n * (n + 1) / 6);
    for (std::size_t m = 0; m < n; m++)
    {
        for (std::size_t l = 0; l < m; l++)
        {
            std::size_t id = (n - 2) * (n - 1) * n / 6 + (m - 1) * m / 2 + l;
            Mdouble invdet = 1.0 / (+wallObjects_[n].getNormal().X * (wallObjects_[m].getNormal().Y * wallObjects_[l].getNormal().Z - wallObjects_[l].getNormal().Y * wallObjects_[m].getNormal().Z)
                    - wallObjects_[n].getNormal().Y * (wallObjects_[m].getNormal().X * wallObjects_[l].getNormal().Z - wallObjects_[m].getNormal().Z * wallObjects_[l].getNormal().X)
                    + wallObjects_[n].getNormal().Z * (wallObjects_[m].getNormal().X * wallObjects_[l].getNormal().Y - wallObjects_[m].getNormal().Y * wallObjects_[l].getNormal().X));
            C_[id] = Vec3D(+(wallObjects_[m].getNormal().Y * wallObjects_[l].getNormal().Z - wallObjects_[l].getNormal().Y * wallObjects_[m].getNormal().Z) * Vec3D::dot(wallObjects_[n].getPosition(),wallObjects_[n].getNormal())
                    - (wallObjects_[n].getNormal().Y * wallObjects_[l].getNormal().Z - wallObjects_[n].getNormal().Z * wallObjects_[l].getNormal().Y) * Vec3D::dot(wallObjects_[m].getPosition(),wallObjects_[m].getNormal())
                    + (wallObjects_[n].getNormal().Y * wallObjects_[m].getNormal().Z - wallObjects_[n].getNormal().Z * wallObjects_[m].getNormal().Y) * Vec3D::dot(wallObjects_[l].getPosition(),wallObjects_[l].getNormal()),
                    -(wallObjects_[m].getNormal().X * wallObjects_[l].getNormal().Z - wallObjects_[m].getNormal().Z * wallObjects_[l].getNormal().X) * Vec3D::dot(wallObjects_[n].getPosition(),wallObjects_[n].getNormal())
                            + (wallObjects_[n].getNormal().X * wallObjects_[l].getNormal().Z - wallObjects_[n].getNormal().Z * wallObjects_[l].getNormal().X) * Vec3D::dot(wallObjects_[m].getPosition(),wallObjects_[m].getNormal())
                            - (wallObjects_[n].getNormal().X * wallObjects_[m].getNormal().Z - wallObjects_[m].getNormal().X * wallObjects_[n].getNormal().Z) * Vec3D::dot(wallObjects_[l].getPosition(),wallObjects_[l].getNormal()),
                    +(wallObjects_[m].getNormal().X * wallObjects_[l].getNormal().Y - wallObjects_[l].getNormal().X * wallObjects_[m].getNormal().Y) * Vec3D::dot(wallObjects_[n].getPosition(),wallObjects_[n].getNormal())
                            - (wallObjects_[n].getNormal().X * wallObjects_[l].getNormal().Y - wallObjects_[l].getNormal().X * wallObjects_[n].getNormal().Y) * Vec3D::dot(wallObjects_[m].getPosition(),wallObjects_[m].getNormal())
                            + (wallObjects_[n].getNormal().X * wallObjects_[m].getNormal().Y - wallObjects_[m].getNormal().X * wallObjects_[n].getNormal().Y) * Vec3D::dot(wallObjects_[l].getPosition(),wallObjects_[l].getNormal())) * invdet;
        }
    }

    logger(VERBOSE, "%", *this);
    for (InfiniteWall w : wallObjects_)
        logger(VERBOSE, "wallObject %, %", w.getNormal(), w.getPosition());
    for (Vec3D v : A_)
        logger(VERBOSE, "A %", v);
    for (Vec3D v : AB_)
        logger(VERBOSE, "AB %", v);
    for (Vec3D v : C_)
        logger(VERBOSE, "C %", v);
}

/*!
 * \param[in] normal The normal to the wallObject.
 * \param[in] position The position of the wallObject in the direction of the normal vector.
 */
void RotatingIntersectionOfWalls::addObject(Vec3D normal, Mdouble position)
{
    logger(WARN, "This function is deprecated, use RotatingIntersectionOfWalls::addObject(Vec3D, Vec3D) instead.");
    addObject(normal, position*normal);
}

/*!
 * \param[in] points A vector of 3D-vectors which contains the points between which the polygon is drawn.
 * \param[in] prismAxis A 3D-vector which represents the direction in which the prism is extended infinitely.
 * \details Create an open prism which is a polygon with no connection between the
 * first and last point, and extending infinitely in the other direction, which is
 * defined as PrismAxis. Do this by adding the walls between the consecutive points
 * one by one.
 */
void RotatingIntersectionOfWalls::createOpenPrism(std::vector<Vec3D> points, Vec3D prismAxis)
{
   // clear();
    //note: use i+1 < points.size() instead of i < points.size()-1, otherwise it creates havoc if point has zero entries.
    for (unsigned int i = 0; i+1 < points.size(); i++)
        addObject(Vec3D::cross(points[i] - points[i + 1], prismAxis), points[i]);
}

/*!
 * \param[in] points A vector of 3D-vectors which contains the points between which the polygon is drawn.
 * \param[in] prismAxis A 3D-vector which represents the direction in which the prism is extended infinitely.
 * \details Create an open prism which is a polygon and extending infinitely in the other direction, which is
 * defined as PrismAxis. Do this by first creating an open prism and then connect
 * the last and the first point.
 */
void RotatingIntersectionOfWalls::createPrism(std::vector<Vec3D> points, Vec3D prismAxis)
{
    createOpenPrism(points, prismAxis);
    addObject(Vec3D::cross(points.back() - points.front(), prismAxis), points.front());
}

/*!
 * \param[in] points A vector of 3D-vectors which contains the points between which the polygon is drawn.
 * \details Create an open prism which is a polygon with no connection between the
 * first and last point, and extending infinitely in the direction perpendicular
 * to the first and second wall. Do this by first computing in which direction the
 * wall must be extended infinitely, then call createOpenPrism(points, prismAxis).
 */
void RotatingIntersectionOfWalls::createOpenPrism(std::vector<Vec3D> points)
{
    Vec3D prismAxis = Vec3D::cross(
            Vec3D::getUnitVector(points[1] - points[0]),
            Vec3D::getUnitVector(points[2] - points[0]));
    createOpenPrism(points, prismAxis);
}

/*!
 * \param[in] points A vector of 3D-vectors which contains the points between which the polygon is drawn.
 * \details Create an open prism which is a polygon, and extending infinitely in
 *  the direction perpendicular to the first and second wall. Do this by first
 * computing in which direction the wall must be extended infinitely, then call
 * createOpenPrism(points, prismAxis).
 */
void RotatingIntersectionOfWalls::createPrism(std::vector<Vec3D> points)
{
    Vec3D prismAxis = Vec3D::cross(
            Vec3D::getUnitVector(points[1] - points[0]),
            Vec3D::getUnitVector(points[2] - points[0]));
    createPrism(points, prismAxis);
}

/*!
 * \param[in] p BaseParticle we want to calculate the distance and whether it collided of.
 * \param[out] distance The distance of the BaseParticle to this wall.
 * \param[out] normal_return If there was a collision, the normal vector to this wall will be placed here.
 * \return A boolean which says whether or not there was a collision.
 * \details This function computes whether or not there is a collision between
 * a given BaseParticle and this RotatingIntersectionOfWalls. If there is a collision, this
 * function also computes the distance between the BaseParticle and RotatingIntersectionOfWalls
 * and the normal of the RotatingIntersectionOfWalls at the intersection point. It does
 * this by calling RotatingIntersectionOfWalls::getDistanceAndNormal(const Vec3D& , Mdouble , Mdouble&, Vec3D&) const.
 * Since this function should be called before calculating any
 * Particle-Wall interactions, it can also be used to set the normal vector in
 * case of curved walls.
 */
bool RotatingIntersectionOfWalls::getDistanceAndNormal(const BaseParticle& p, Mdouble& distance, Vec3D& normal_return) const
{
    return getDistanceAndNormal(p.getPosition(), p.getInteractionRadius(), distance, normal_return);
}

/*!
 * \param[in] position The position of the object there is possible an interaction with.
 * \param[in] wallInteractionRadius The maximum distance between the RotatingIntersectionOfWalls
 *  and the input argument position for which there is an interaction.
 * \param[out] distance The distance of the object at position to this wall.
 * \param[out] normal_return If there was an interaction, the normal vector to
 * this wall will be placed here.
 * \return A boolean which says whether or not there was an interaction.
 * \details This function computes whether or not there is an interaction between
 * an object at the given distance and this RotatingIntersectionOfWalls. If there is an interaction, this
 * function also computes the distance between the BaseParticle and RotatingIntersectionOfWalls
 * and the normal of the RotatingIntersectionOfWalls at the intersection point.
 * First check if the distance between the object at position and the
 * RotatingIntersectionOfWalls is smaller or greater than the wallInteractionRadius. If
 * there is no interaction, return false, the output parameters then have no meaning.
 * If there is an interaction, find out which (one or more) of the InfiniteWall
 * there is an interaction with. Then compute the distance between the particle
 * and InfiniteWall and the normal to the interaction point.
 * Since this function should be called before calculating any
 * Particle-Wall interactions, it can also be used to set the normal vector in
 * case of curved walls.
 */
bool RotatingIntersectionOfWalls::getDistanceAndNormal(const Vec3D& position, Mdouble wallInteractionRadius, Mdouble &distance, Vec3D &normal_return) const
{
    if (wallObjects_.size()==0)
    {
        logger(DEBUG,"Empty RotatingIntersectionOfWalls");
        return false;
    }

    distance = -1e20;
    Mdouble distance2 = -1e20;
    Mdouble distance3 = -1e20;
    Mdouble distanceCurrent;
    unsigned int id=0;
    unsigned int id2=0;
    unsigned int id3=0;

    //The object has to touch each wall  each wall (distanceCurrent) and keep the minimum distance (distance) and wall index (id)
    for (unsigned int i = 0; i < wallObjects_.size(); i++)
    {
        // Calculate distance to each wall (distanceCurrent);
        distanceCurrent = wallObjects_[i].getDistance(position);
        // The object has to touch each wall (distanceCurrent >= wallInteractionRadius), otherwise return false (i.e. no contact)
        // This means that for each InfiniteWall in wallObjects_, the particle is either "inside"
        // the wall or touching it. If not, there is no interaction.
        if (distanceCurrent >= wallInteractionRadius)
            return false;
        // Find out which of the InfiniteWalls is interacting with the particle.
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
                    id3 = id;
                }
                else
                {
                    distance2 = distance;
                    id2 = id;
                }
            }
            distance = distanceCurrent;
            id = i;
        }
        else if (distanceCurrent > -wallInteractionRadius)
        {
            if (distance2 > -wallInteractionRadius)
            {
                distance3 = distanceCurrent;
                id3 = i;
            }
            else
            {
                distance2 = distanceCurrent;
                id2 = i;
            }
        }
    }

    //If we are here, the closest wall is id;
    //if distance2>-P.Radius (and distance3>-P.Radius), the possible contact point
    // is near the intersection between id and id2 (and id3)
    if (distance2 > -wallInteractionRadius)
    {
        //D is the point on wall id closest to P
        Vec3D D = position + wallObjects_[id].getNormal() * distance;
        //If the distance of D to id2 is positive, the contact is with the intersection
        bool intersection_with_id2 = (wallObjects_[id2].getDistance(D) > 0.0);

        if (distance3 > -wallInteractionRadius && (wallObjects_[id3].getDistance(D) > 0.0))
        {
            if (intersection_with_id2)
            {
                //possible contact is with intersection of id,id2,id3
                //we know id2<id3
                unsigned int index =
                        (id < id2) ? ((id3 - 2) * (id3 - 1) * id3 / 6 + (id2 - 1) * id2 / 2 + id) :
                        (id < id3) ? ((id3 - 2) * (id3 - 1) * id3 / 6 + (id - 1) * id / 2 + id2) :
                                     ((id - 2) * (id - 1) * id / 6 + (id3 - 1) * id3 / 2 + id2);
                normal_return = position - C_[index];
                distance = sqrt(normal_return.getLengthSquared());
                if (distance >= wallInteractionRadius)
                    return false; //no contact
                normal_return /= -distance;
                return true; //contact with id,id2,id3
            }
            else
            {
                intersection_with_id2 = true;
                distance2 = distance3;
                id2 = id3;
            }
        }

        if (intersection_with_id2)
        { //possible contact is with intersection of id,id2
            unsigned int index = (id > id2) ? ((id - 1) * id / 2 + id2) : ((id2 - 1) * id2 / 2 + id);
            Vec3D AC = position - A_[index];
            normal_return = AC - AB_[index] * Vec3D::dot(AC, AB_[index]);
            distance = sqrt(normal_return.getLengthSquared());
            if (distance >= wallInteractionRadius)
                return false; //no contact
            normal_return /= -distance;
            return true; //contact with two walls
        }
    }
    //contact is with id
    normal_return = wallObjects_[id].getNormal();
    return true;
}

/*!
 * \param[in] move A reference to a Vec3D that denotes the direction and length
 * it should be moved with.
 * \details A function that moves the InterSectionOfWalls in a certain direction
 * by both moving the walls and all intersections. Note that the directions of the
 * intersections are not moved since they don't change when moving the RotatingIntersectionOfWalls
 * as a whole.
 * \todo We should use the position_ and orientation_ of the RotatingIntersectionOfWalls;
 * that way, RotatingIntersectionOfWalls can be moved with the standard BaseInteractable::move function,
 * getting rid of an anomaly in the code and removing the virtual from the move function. \author weinhartt
 */
void RotatingIntersectionOfWalls::move(const Vec3D& move)
{
    BaseInteractable::move(move);
    for(Vec3D& a : A_)
    {
		a += move;
	}
	for(Vec3D& c : C_)
	{
		c += move;
	}
	for(InfiniteWall& o : wallObjects_)
	{
		o.move(move);
	}
}

// HACK - BEGIN
void RotatingIntersectionOfWalls::moveAlongZ(double displacement)
{
    // number of infinite walls composing the intersection of walls
    int nOfWalls = wallObjects_.size();
    
    Vec3D newPoint;
    
    // loops through all walls
    for (int ii = 1; ii < nOfWalls + 1; ii++)
    {
        // if the normal has components along Z then move, otherwise skip the traslation
        if (wallObjects_[ii - 1].getNormal().Z)
        {
            // compute the new position
            newPoint.X = wallObjects_[ii - 1].getPosition().X;
            newPoint.Y = wallObjects_[ii - 1].getPosition().Y;
            newPoint.Z = wallObjects_[ii - 1].getPosition().Z + displacement;
            
            wallObjects_[ii - 1].set(wallObjects_[ii - 1].getNormal(), newPoint);
        }
        else
        {
            wallObjects_[ii - 1].set(wallObjects_[ii - 1].getNormal(), wallObjects_[ii - 1].getPosition());
        }
        
        // re-do everything done in the RotatingIntersectionOfWalls::addObject(Vec3D normal, Vec3D point) function
        AB_.resize(ii * (ii + 1) / 2); //0 + 1 + 2 + ... + indexNew, total number of walls you need
        A_.resize(ii * (ii + 1) / 2);
        for (std::size_t m = 0; m < ii; m++)
        {
            std::size_t id = (ii - 1) * ii / 2 + m;
            //first we cross the wall normals and normalize to obtain AB
            AB_[id] = Vec3D::cross(wallObjects_[m].getNormal(), wallObjects_[ii].getNormal());
            AB_[id] /= sqrt(AB_[id].getLengthSquared());
            //then we find a point A (using AB*x=0 as a third plane)
            Mdouble invdet = 1.0 / (+wallObjects_[ii].getNormal().X * (wallObjects_[m].getNormal().Y * AB_[id].Z - AB_[id].Y * wallObjects_[m].getNormal().Z)
                                    - wallObjects_[ii].getNormal().Y * (wallObjects_[m].getNormal().X * AB_[id].Z - wallObjects_[m].getNormal().Z * AB_[id].X)
                                    + wallObjects_[ii].getNormal().Z * (wallObjects_[m].getNormal().X * AB_[id].Y - wallObjects_[m].getNormal().Y * AB_[id].X));
            
            A_[id] = Vec3D( +(wallObjects_[m].getNormal().Y * AB_[id].Z - AB_[id].Y * wallObjects_[m].getNormal().Z) * Vec3D::dot(wallObjects_[ii].getPosition(),wallObjects_[ii].getNormal())
                           -(wallObjects_[ii].getNormal().Y * AB_[id].Z - wallObjects_[ii].getNormal().Z * AB_[id].Y) * Vec3D::dot(wallObjects_[m].getPosition(),wallObjects_[m].getNormal())
                           +(wallObjects_[ii].getNormal().Y * wallObjects_[m].getNormal().Z - wallObjects_[ii].getNormal().Z * wallObjects_[m].getNormal().Y) * 0.0,
                           -(wallObjects_[m].getNormal().X * AB_[id].Z - wallObjects_[m].getNormal().Z * AB_[id].X) * Vec3D::dot(wallObjects_[ii].getPosition(),wallObjects_[ii].getNormal())
                           +(wallObjects_[ii].getNormal().X * AB_[id].Z - wallObjects_[ii].getNormal().Z * AB_[id].X) * Vec3D::dot(wallObjects_[m].getPosition(),wallObjects_[m].getNormal())
                           -(wallObjects_[ii].getNormal().X * wallObjects_[m].getNormal().Z - wallObjects_[m].getNormal().X * wallObjects_[ii].getNormal().Z) * 0.0,
                           +(wallObjects_[m].getNormal().X * AB_[id].Y - AB_[id].X * wallObjects_[m].getNormal().Y) * Vec3D::dot(wallObjects_[ii].getPosition(),wallObjects_[ii].getNormal())
                           -(wallObjects_[ii].getNormal().X * AB_[id].Y - AB_[id].X * wallObjects_[ii].getNormal().Y) * Vec3D::dot(wallObjects_[m].getPosition(),wallObjects_[m].getNormal())
                           +(wallObjects_[ii].getNormal().X * wallObjects_[m].getNormal().Y - wallObjects_[m].getNormal().X * wallObjects_[ii].getNormal().Y) * 0.0) * invdet;
        }
        
        C_.resize((ii - 1) * ii * (ii + 1) / 6);
        for (std::size_t m = 0; m < ii; m++)
        {
            for (std::size_t l = 0; l < m; l++)
            {
                std::size_t id = (ii - 2) * (ii - 1) * ii / 6 + (m - 1) * m / 2 + l;
                Mdouble invdet = 1.0 / (+wallObjects_[ii].getNormal().X * (wallObjects_[m].getNormal().Y * wallObjects_[l].getNormal().Z - wallObjects_[l].getNormal().Y * wallObjects_[m].getNormal().Z)
                                        - wallObjects_[ii].getNormal().Y * (wallObjects_[m].getNormal().X * wallObjects_[l].getNormal().Z - wallObjects_[m].getNormal().Z * wallObjects_[l].getNormal().X)
                                        + wallObjects_[ii].getNormal().Z * (wallObjects_[m].getNormal().X * wallObjects_[l].getNormal().Y - wallObjects_[m].getNormal().Y * wallObjects_[l].getNormal().X));
                C_[id] = Vec3D(+(wallObjects_[m].getNormal().Y * wallObjects_[l].getNormal().Z - wallObjects_[l].getNormal().Y * wallObjects_[m].getNormal().Z) * Vec3D::dot(wallObjects_[ii].getPosition(),wallObjects_[ii].getNormal())
                               - (wallObjects_[ii].getNormal().Y * wallObjects_[l].getNormal().Z - wallObjects_[ii].getNormal().Z * wallObjects_[l].getNormal().Y) * Vec3D::dot(wallObjects_[m].getPosition(),wallObjects_[m].getNormal())
                               + (wallObjects_[ii].getNormal().Y * wallObjects_[m].getNormal().Z - wallObjects_[ii].getNormal().Z * wallObjects_[m].getNormal().Y) * Vec3D::dot(wallObjects_[l].getPosition(),wallObjects_[l].getNormal()),
                               -(wallObjects_[m].getNormal().X * wallObjects_[l].getNormal().Z - wallObjects_[m].getNormal().Z * wallObjects_[l].getNormal().X) * Vec3D::dot(wallObjects_[ii].getPosition(),wallObjects_[ii].getNormal())
                               + (wallObjects_[ii].getNormal().X * wallObjects_[l].getNormal().Z - wallObjects_[ii].getNormal().Z * wallObjects_[l].getNormal().X) * Vec3D::dot(wallObjects_[m].getPosition(),wallObjects_[m].getNormal())
                               - (wallObjects_[ii].getNormal().X * wallObjects_[m].getNormal().Z - wallObjects_[m].getNormal().X * wallObjects_[ii].getNormal().Z) * Vec3D::dot(wallObjects_[l].getPosition(),wallObjects_[l].getNormal()),
                               +(wallObjects_[m].getNormal().X * wallObjects_[l].getNormal().Y - wallObjects_[l].getNormal().X * wallObjects_[m].getNormal().Y) * Vec3D::dot(wallObjects_[ii].getPosition(),wallObjects_[ii].getNormal())
                               - (wallObjects_[ii].getNormal().X * wallObjects_[l].getNormal().Y - wallObjects_[l].getNormal().X * wallObjects_[ii].getNormal().Y) * Vec3D::dot(wallObjects_[m].getPosition(),wallObjects_[m].getNormal())
                               + (wallObjects_[ii].getNormal().X * wallObjects_[m].getNormal().Y - wallObjects_[m].getNormal().X * wallObjects_[ii].getNormal().Y) * Vec3D::dot(wallObjects_[l].getPosition(),wallObjects_[l].getNormal())) * invdet;
            }
        }
    }
}
    
void RotatingIntersectionOfWalls::rotateAroundZ(double angle)
{
   // number of infinite walls composing the intersection of walls
   int nOfWalls = wallObjects_.size();

   Vec3D newNormal;
   Vec3D newPoint;

   // loops through all walls
   for (int ii = 1; ii < nOfWalls + 1; ii++)
   {
      // if the normal has components along X or Y then rotate, otherwise skip the rotation
      if (wallObjects_[ii - 1].getNormal().X || wallObjects_[ii - 1].getNormal().X)
      {
         // compute the new normal
         newNormal.setZero();
         newNormal.X = std::cos(angle)*wallObjects_[ii - 1].getNormal().X - std::sin(angle)*wallObjects_[ii - 1].getNormal().Y;
         newNormal.Y = std::sin(angle)*wallObjects_[ii - 1].getNormal().X + std::cos(angle)*wallObjects_[ii - 1].getNormal().Y;

         // compute the new position
         newPoint.X = std::cos(angle)*wallObjects_[ii - 1].getPosition().X - std::sin(angle)*wallObjects_[ii - 1].getPosition().Y;
         newPoint.Y = std::sin(angle)*wallObjects_[ii - 1].getPosition().X + std::cos(angle)*wallObjects_[ii - 1].getPosition().Y;
         newPoint.Z = wallObjects_[ii - 1].getPosition().Z;

         wallObjects_[ii - 1].set(newNormal, newPoint);
      }
      else
      {
         wallObjects_[ii - 1].set(wallObjects_[ii - 1].getNormal(), wallObjects_[ii - 1].getPosition());
      }

      // re-do everything done in the RotatingIntersectionOfWalls::addObject(Vec3D normal, Vec3D point) function
      AB_.resize(ii * (ii + 1) / 2); //0 + 1 + 2 + ... + indexNew, total number of walls you need
      A_.resize(ii * (ii + 1) / 2);
      for (std::size_t m = 0; m < ii; m++)
      {
         std::size_t id = (ii - 1) * ii / 2 + m;
         //first we cross the wall normals and normalize to obtain AB
         AB_[id] = Vec3D::cross(wallObjects_[m].getNormal(), wallObjects_[ii].getNormal());
         AB_[id] /= sqrt(AB_[id].getLengthSquared());
         //then we find a point A (using AB*x=0 as a third plane)
         Mdouble invdet = 1.0 / (+wallObjects_[ii].getNormal().X * (wallObjects_[m].getNormal().Y * AB_[id].Z - AB_[id].Y * wallObjects_[m].getNormal().Z)
                  - wallObjects_[ii].getNormal().Y * (wallObjects_[m].getNormal().X * AB_[id].Z - wallObjects_[m].getNormal().Z * AB_[id].X)
                  + wallObjects_[ii].getNormal().Z * (wallObjects_[m].getNormal().X * AB_[id].Y - wallObjects_[m].getNormal().Y * AB_[id].X));

         A_[id] = Vec3D( +(wallObjects_[m].getNormal().Y * AB_[id].Z - AB_[id].Y * wallObjects_[m].getNormal().Z) * Vec3D::dot(wallObjects_[ii].getPosition(),wallObjects_[ii].getNormal())
                          -(wallObjects_[ii].getNormal().Y * AB_[id].Z - wallObjects_[ii].getNormal().Z * AB_[id].Y) * Vec3D::dot(wallObjects_[m].getPosition(),wallObjects_[m].getNormal())
                          +(wallObjects_[ii].getNormal().Y * wallObjects_[m].getNormal().Z - wallObjects_[ii].getNormal().Z * wallObjects_[m].getNormal().Y) * 0.0,
                          -(wallObjects_[m].getNormal().X * AB_[id].Z - wallObjects_[m].getNormal().Z * AB_[id].X) * Vec3D::dot(wallObjects_[ii].getPosition(),wallObjects_[ii].getNormal())
                          +(wallObjects_[ii].getNormal().X * AB_[id].Z - wallObjects_[ii].getNormal().Z * AB_[id].X) * Vec3D::dot(wallObjects_[m].getPosition(),wallObjects_[m].getNormal())
                          -(wallObjects_[ii].getNormal().X * wallObjects_[m].getNormal().Z - wallObjects_[m].getNormal().X * wallObjects_[ii].getNormal().Z) * 0.0,
                          +(wallObjects_[m].getNormal().X * AB_[id].Y - AB_[id].X * wallObjects_[m].getNormal().Y) * Vec3D::dot(wallObjects_[ii].getPosition(),wallObjects_[ii].getNormal())
                          -(wallObjects_[ii].getNormal().X * AB_[id].Y - AB_[id].X * wallObjects_[ii].getNormal().Y) * Vec3D::dot(wallObjects_[m].getPosition(),wallObjects_[m].getNormal())
                          +(wallObjects_[ii].getNormal().X * wallObjects_[m].getNormal().Y - wallObjects_[m].getNormal().X * wallObjects_[ii].getNormal().Y) * 0.0) * invdet;
      }

      C_.resize((ii - 1) * ii * (ii + 1) / 6);
      for (std::size_t m = 0; m < ii; m++)
      {
         for (std::size_t l = 0; l < m; l++)
         {
              std::size_t id = (ii - 2) * (ii - 1) * ii / 6 + (m - 1) * m / 2 + l;
              Mdouble invdet = 1.0 / (+wallObjects_[ii].getNormal().X * (wallObjects_[m].getNormal().Y * wallObjects_[l].getNormal().Z - wallObjects_[l].getNormal().Y * wallObjects_[m].getNormal().Z)
                      - wallObjects_[ii].getNormal().Y * (wallObjects_[m].getNormal().X * wallObjects_[l].getNormal().Z - wallObjects_[m].getNormal().Z * wallObjects_[l].getNormal().X)
                      + wallObjects_[ii].getNormal().Z * (wallObjects_[m].getNormal().X * wallObjects_[l].getNormal().Y - wallObjects_[m].getNormal().Y * wallObjects_[l].getNormal().X));
              C_[id] = Vec3D(+(wallObjects_[m].getNormal().Y * wallObjects_[l].getNormal().Z - wallObjects_[l].getNormal().Y * wallObjects_[m].getNormal().Z) * Vec3D::dot(wallObjects_[ii].getPosition(),wallObjects_[ii].getNormal())
                      - (wallObjects_[ii].getNormal().Y * wallObjects_[l].getNormal().Z - wallObjects_[ii].getNormal().Z * wallObjects_[l].getNormal().Y) * Vec3D::dot(wallObjects_[m].getPosition(),wallObjects_[m].getNormal())
                      + (wallObjects_[ii].getNormal().Y * wallObjects_[m].getNormal().Z - wallObjects_[ii].getNormal().Z * wallObjects_[m].getNormal().Y) * Vec3D::dot(wallObjects_[l].getPosition(),wallObjects_[l].getNormal()),
                      -(wallObjects_[m].getNormal().X * wallObjects_[l].getNormal().Z - wallObjects_[m].getNormal().Z * wallObjects_[l].getNormal().X) * Vec3D::dot(wallObjects_[ii].getPosition(),wallObjects_[ii].getNormal())
                              + (wallObjects_[ii].getNormal().X * wallObjects_[l].getNormal().Z - wallObjects_[ii].getNormal().Z * wallObjects_[l].getNormal().X) * Vec3D::dot(wallObjects_[m].getPosition(),wallObjects_[m].getNormal())
                              - (wallObjects_[ii].getNormal().X * wallObjects_[m].getNormal().Z - wallObjects_[m].getNormal().X * wallObjects_[ii].getNormal().Z) * Vec3D::dot(wallObjects_[l].getPosition(),wallObjects_[l].getNormal()),
                      +(wallObjects_[m].getNormal().X * wallObjects_[l].getNormal().Y - wallObjects_[l].getNormal().X * wallObjects_[m].getNormal().Y) * Vec3D::dot(wallObjects_[ii].getPosition(),wallObjects_[ii].getNormal())
                              - (wallObjects_[ii].getNormal().X * wallObjects_[l].getNormal().Y - wallObjects_[l].getNormal().X * wallObjects_[ii].getNormal().Y) * Vec3D::dot(wallObjects_[m].getPosition(),wallObjects_[m].getNormal())
                              + (wallObjects_[ii].getNormal().X * wallObjects_[m].getNormal().Y - wallObjects_[m].getNormal().X * wallObjects_[ii].getNormal().Y) * Vec3D::dot(wallObjects_[l].getPosition(),wallObjects_[l].getNormal())) * invdet;
         }
      }
   }
}
// HACK - END

/*!
 * \param[in] is The input stream from which the RotatingIntersectionOfWalls is read, usually a restart file.
 */
void RotatingIntersectionOfWalls::read(std::istream& is)
{
    BaseWall::read(is);
    std::string dummy;
    int n;
    is >> dummy >> n;

    Vec3D normal;
    Vec3D position;
    for (int i = 0; i < n; i++)
    {
        is >> dummy >> normal >> dummy >> position;
        addObject(normal, position);
    }
}

/*!
 * \param[in] os The output stream where the RotatingIntersectionOfWalls must be written
 *  to, usually a restart file.
 */
void RotatingIntersectionOfWalls::write(std::ostream& os) const
{
    BaseWall::write(os);
    os << " numIntersectionOfWalls " << wallObjects_.size();
    for (std::vector<InfiniteWall>::const_iterator it = wallObjects_.begin(); it != wallObjects_.end(); ++it)
    {
        os << " normal " << it->getNormal() << " position " << it->getPosition();
    }
}

/*!
 * \return The string "RotatingIntersectionOfWalls".
 */
std::string RotatingIntersectionOfWalls::getName() const
{
    return "RotatingIntersectionOfWalls";
}

void RotatingIntersectionOfWalls::writeVTK (VTKContainer& vtk) const
{
    for (auto wall=wallObjects_.begin(); wall!=wallObjects_.end(); wall++)
    {
        std::vector<Vec3D> points;
        wall->createVTK (points);
        for (auto other=wallObjects_.begin(); other!=wallObjects_.end(); other++)
        {
            if (other!=wall)
            {
                intersectVTK(points, -other->getNormal(), other->getPosition());
            }
        }
        wall->addToVTK (points, vtk);
    }
}
