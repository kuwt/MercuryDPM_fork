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

#include "IntersectionOfWalls.h"
#include "InteractionHandler.h"
#include "WallHandler.h"
#include "DPMBase.h"
//#include "Species/BaseSpecies.h"
#include "Particles/BaseParticle.h"

IntersectionOfWalls::IntersectionOfWalls()
{
    logger(DEBUG, "IntersectionOfWalls() constructed.");
}

/*!
 * \param[in] other The IntersectionOfWalls that must be copied.
 */
IntersectionOfWalls::IntersectionOfWalls(const IntersectionOfWalls& other)
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
    logger(DEBUG, "IntersectionOfWalls(IntersectionOfWalls&) constructed.");
}

IntersectionOfWalls::IntersectionOfWalls(const std::vector<normalAndPosition>& walls,
                                         const ParticleSpecies* species)
{
    setSpecies(species);
    for (auto wall : walls)
    {
        addObject(wall.normal, wall.position);
    }
}


IntersectionOfWalls::~IntersectionOfWalls()
{
    logger(DEBUG, "~IntersectionOfWalls() has been called.");
}

void IntersectionOfWalls::setSpecies(const ParticleSpecies* species)
{
    BaseWall::setSpecies(species);
    for (auto wall : wallObjects_)
    {
        wall.setSpecies(species);
    }
}

/*!
 * \param[in] other The IntersectionOfWalls that must be copied.
 */
IntersectionOfWalls& IntersectionOfWalls::operator=(const IntersectionOfWalls& other)
{
    logger(DEBUG, "IntersectionOfWalls::operator= called.");
    if (this == &other)
    {
        return *this;
    }
    return *(other.copy());
}

/*!
 * \return pointer to a IntersectionOfWalls object allocated using new.
 */
IntersectionOfWalls* IntersectionOfWalls::copy() const
{
    return new IntersectionOfWalls(*this);
}

void IntersectionOfWalls::clear()
{
    wallObjects_.clear();
    A_.clear();
    AB_.clear();
    C_.clear();
}

void IntersectionOfWalls::setHandler(WallHandler* wallHandler)
{
    BaseWall::setHandler(wallHandler);
    for (InfiniteWall& w : wallObjects_)
    {
        w.setHandler(wallHandler);
    }
}

/*! 
 * \detail Suppose your simulation adds to an IntersectionOfWalls after a
 * certain time or condition is met. Checking the number of objects is useful
 * for checking if this has happened yet, when restarting. 
 */
unsigned int IntersectionOfWalls::getNumberOfObjects()
{
    return wallObjects_.size();
}

/*!
 * \param[in] normal The normal to this wallObject.
 * \param[in] point One of the points of the wallObject.
 * \details Adds a wall to the set of finite walls, given an outward unit normal
 *  vector s.t. normal*x=normal*point for all x of the wallObject. First make the
 * InfiniteWall, then compute all intersections, which are then stored in A_, AB_
 * and C_.
 */
void IntersectionOfWalls::addObject(Vec3D normal, Vec3D point)
{
    normal.normalise();
    
    //n is the index of the new wall
    std::size_t n = wallObjects_.size();
    InfiniteWall w;
    if (getSpecies()) w.setSpecies(getSpecies());
    w.set(normal, point);
    wallObjects_.push_back(w);
    setPointsAndLines(n);
}

/*!
 *
 * \param PointA first coordinate plane passes through;
 * \param PointB second coordinate plane passes through;
 * \param PointC third coordinate plane passes through;
 * \details calls the IntersectionOfWalls::addObject(normal, point) by calculating the normal from these three coordinates
 *
 */
void IntersectionOfWalls::add3PointObject(Vec3D PointA, Vec3D PointB, Vec3D PointC)
{
    Vec3D SubtB = PointA - PointB;
    Vec3D SubtC = PointA - PointC;
    
    Vec3D WallNormal = Vec3D::cross(SubtB, SubtC);
    //Check if walls coordinates  inline, if true Do not build wall and give error message.
    if (WallNormal.getLengthSquared() == 0.0)
    {
        logger(ERROR,
               "Error Building IntersectionOfWalls::add3PointObject out of 3 coordinates. Coordinates are in line, Wall not constructed.");
    }
    else
    {
        addObject(WallNormal, PointA);
    }
    
}


/*!
 *
 * \param PointA first coordinate plane passes through;
 * \param PointB second coordinate plane passes through;
 * \param PointC third coordinate plane passes through;
 * \param WallNormal the normal of the wal of the plane, Please note that the Wallnormal input is defined inverse to other Mercury functions
 * Due to the usage of this function for reading in STL files
 * \param Thickness the height of the apex
 * \param wallidentifier to identify which wall does not get constructed for bug finding in STL files
 * \details constructs a tethrahedron with the apex in minus normal direction.
 *
 */
void IntersectionOfWalls::addTetraSTL(Vec3D PointA, Vec3D PointB, Vec3D PointC, Vec3D WallNormal, Mdouble Thickness,
                                      int wallidentifier)
{
    //Check if wall coordinates  inline, if true Do not build wall and give error message. But keep continuing
    if (WallNormal.getLengthSquared() == 0.0)
    {
        std::cout << "Error Building Plate number " << wallidentifier
                  << "out of 3 coordinates. Coordinates are in line, Wall not constructed." << std::endl;
    }
    else
    {
        //do a check whether the normal follows the RHR (Right Hand Rule) or not
        //todo: Bert Use this check a lot, possibly make a function of this
        Vec3D SubtB = PointA - PointB;
        Vec3D SubtC = PointA - PointC;
        
        Vec3D WallNormalRHR = Vec3D::cross(SubtB, SubtC);
        
        //normalise for easy check
        WallNormalRHR.normalise();
        WallNormal.normalise();
        
        //if RHRchecl is 1, wall normal and RHR normal are in same direction, if -1 theyre not, then point B and C need to be swapped
        Mdouble RHRcheck = Vec3D::dot(WallNormalRHR, WallNormal);
        //todo: Bert Officially need to check for other answers than 1, however, it will either be -1 or 1 in ideal case
        
        if (RHRcheck == 1)
        {
            //calculate centroid
            Vec3D mid = (PointA + PointB + PointC) / 3.0;
            //shift centroid in normal direction
            //mental note: if Same direction as STL normal it should be subtracted to go IN the wall
            Vec3D midT = mid - WallNormalRHR * Thickness;
            //generate the base wall first through the input coordinates
            add3PointObject(PointA, PointC, PointB);
            
            //add sidewalls
            add3PointObject(PointA, midT, PointC);
            add3PointObject(PointC, midT, PointB);
            add3PointObject(PointB, midT, PointA);
        }
        else
        {
            //calculate centroid
            Vec3D mid = (PointA + PointB + PointC) / 3.0;
            //shift centroid in normal direction
            //mental note: if opposite direction as STL normal it should be added to go IN the wall
            Vec3D midT = mid + WallNormalRHR * Thickness;
            //generate the base wall first through the input coordinates
            add3PointObject(PointA, PointB, PointC);
            
            //add sidewalls
            add3PointObject(PointA, midT, PointB);
            add3PointObject(PointB, midT, PointC);
            add3PointObject(PointC, midT, PointA);
        }
        
    }
}


/*!
 *
 * \param PointA first coordinate plane passes through;
 * \param PointB second coordinate plane passes through;
 * \param PointC third coordinate plane passes through;
 * \param Thickness the height of the apex
 * \details constructs a tethrahedron with the apex in the direction according to the right hand rule
 */
void IntersectionOfWalls::addTetra(const Vec3D& PointA, const Vec3D& PointB, const Vec3D& PointC, Mdouble& Thickness)
{
    
    //generate other coordinate by finding the centre and shift it the distance thickness in the normal direction
    //calculate normal
    Vec3D SubtB = PointA - PointB;
    Vec3D SubtC = PointA - PointC;
    
    Vec3D WallNormal = Vec3D::cross(SubtB, SubtC);
    
    //Check if walls coordinates  inline, if true Do not build wall and give error message.
    if (WallNormal.getLengthSquared() == 0.0)
    {
        logger(ERROR,
               "Error Building IntersectionOfWalls::addTetra out of 3 coordinates. "
               "Coordinates are in line, Wall not constructed.");
    }
    else
    {
        //calculate centroid
        Vec3D mid = (PointA + PointB + PointC) / 3.0;
        //shift centroid in normal direction
        Vec3D midT = mid + WallNormal * Thickness;
        //generate the base wall first through the input coordinates
        add3PointObject(PointC, PointB, PointA);
        
        //add sidewalls
        add3PointObject(PointA, midT, PointB);
        add3PointObject(PointB, midT, PointC);
        add3PointObject(PointC, midT, PointA);
    }
}


//funtion builds a triangle wall through the 3 points with a thickness
//note normal is defined OPPOSITE as to normal Mercury convention due to STL convention
void
IntersectionOfWalls::addPlate(const Vec3D& PointA, const Vec3D& PointB, const Vec3D& PointC, const Vec3D& WallNormal,
                              const Mdouble& Thickness, int wallidentifier)
{
    
    //Check if walls coordinates  inline, if true Do not build wall and give error message. But keep continuing
    if (WallNormal.getLengthSquared() == 0.0)
    {
        std::cout << "Error Building Plate number " << wallidentifier
                  << "out of 3 coordinates. Coordinates are in line, Wall not constructed." << std::endl;
    }
    else
    {
        Vec3D PointAT = PointA - WallNormal * Thickness;
        Vec3D PointBT = PointB - WallNormal * Thickness;
        Vec3D PointCT = PointC - WallNormal * Thickness;
        
        //generate the base wall first through the input coordinates
        add3PointObject(PointC, PointB, PointA);
        
        //add sidewalls
        add3PointObject(PointA, PointB, PointAT);
        add3PointObject(PointB, PointC, PointBT);
        add3PointObject(PointC, PointA, PointCT);
        
        
        //add opposite wall
        add3PointObject(PointAT, PointBT, PointCT);
        
        
    }
    
    
}

/*!
 * \param[in] normal The normal to this wallObject.
 * \param[in] position One of the points of the wallObject.
 * \details Adds a wall to the set of finite walls, given an outward unit normal
 *  vector s.t. normal*x=normal*point for all x of the wallObject. First make the
 * InfiniteWall, then compute all intersections, which are then stored in A_, AB_
 * and C_.
 */
void IntersectionOfWalls::addObject(Quaternion orientation, Vec3D position)
{
    //n is the index of the new wall
    std::size_t n = wallObjects_.size();
    InfiniteWall w;
    if (getSpecies()) w.setSpecies(getSpecies());
    w.setOrientation(orientation);
    w.setPosition(position);
    wallObjects_.push_back(w);
    setPointsAndLines(n);
}

void IntersectionOfWalls::setPointsAndLines(unsigned int n)
{
    // AB[n*(n-1)/2+m] is the direction of the intersecting line between walls m and n, m<n
    // A[n*(n-1)/2+m] is a point on the intersecting line between walls m and n, m<n
    // See http://www.netcomuk.co.uk/~jenolive/vect18d.html for finding the line where two planes meet
    AB_.resize(n * (n + 1) / 2); //0 + 1 + 2 + ... + indexNew, total number of walls you need
    A_.resize(n * (n + 1) / 2);
    for (std::size_t m = 0; m < n; m++)
    {
        std::size_t id = (n - 1) * n / 2 + m;
        //first we cross the wall normals and normalise to obtain AB
        AB_[id] = Vec3D::cross(wallObjects_[m].getNormal(), wallObjects_[n].getNormal());
        AB_[id] /= sqrt(AB_[id].getLengthSquared());
        //then we find a point A (using AB*x=0 as a third plane)
        Mdouble invdet = 1.0 / (+wallObjects_[n].getNormal().X *
                                (wallObjects_[m].getNormal().Y * AB_[id].Z - AB_[id].Y * wallObjects_[m].getNormal().Z)
                                - wallObjects_[n].getNormal().Y * (wallObjects_[m].getNormal().X * AB_[id].Z -
                                                                   wallObjects_[m].getNormal().Z * AB_[id].X)
                                + wallObjects_[n].getNormal().Z * (wallObjects_[m].getNormal().X * AB_[id].Y -
                                                                   wallObjects_[m].getNormal().Y * AB_[id].X));
        
        A_[id] = Vec3D(+(wallObjects_[m].getNormal().Y * AB_[id].Z - AB_[id].Y * wallObjects_[m].getNormal().Z) *
                       Vec3D::dot(wallObjects_[n].getPosition(), wallObjects_[n].getNormal())
                       - (wallObjects_[n].getNormal().Y * AB_[id].Z - wallObjects_[n].getNormal().Z * AB_[id].Y) *
                         Vec3D::dot(wallObjects_[m].getPosition(), wallObjects_[m].getNormal())
                       + (wallObjects_[n].getNormal().Y * wallObjects_[m].getNormal().Z -
                          wallObjects_[n].getNormal().Z * wallObjects_[m].getNormal().Y) * 0.0,
                       -(wallObjects_[m].getNormal().X * AB_[id].Z - wallObjects_[m].getNormal().Z * AB_[id].X) *
                       Vec3D::dot(wallObjects_[n].getPosition(), wallObjects_[n].getNormal())
                       + (wallObjects_[n].getNormal().X * AB_[id].Z - wallObjects_[n].getNormal().Z * AB_[id].X) *
                         Vec3D::dot(wallObjects_[m].getPosition(), wallObjects_[m].getNormal())
                       - (wallObjects_[n].getNormal().X * wallObjects_[m].getNormal().Z -
                          wallObjects_[m].getNormal().X * wallObjects_[n].getNormal().Z) * 0.0,
                       +(wallObjects_[m].getNormal().X * AB_[id].Y - AB_[id].X * wallObjects_[m].getNormal().Y) *
                       Vec3D::dot(wallObjects_[n].getPosition(), wallObjects_[n].getNormal())
                       - (wallObjects_[n].getNormal().X * AB_[id].Y - AB_[id].X * wallObjects_[n].getNormal().Y) *
                         Vec3D::dot(wallObjects_[m].getPosition(), wallObjects_[m].getNormal())
                       + (wallObjects_[n].getNormal().X * wallObjects_[m].getNormal().Y -
                          wallObjects_[m].getNormal().X * wallObjects_[n].getNormal().Y) * 0.0) * invdet;
    }
    
    // C[(n-2)*(n-1)*n/6+(m-1)*m/2+l] is a point intersecting walls l, m and n, l<m<n
    C_.resize((n - 1) * n * (n + 1) / 6);
    for (std::size_t m = 0; m < n; m++)
    {
        for (std::size_t l = 0; l < m; l++)
        {
            std::size_t id = (n - 2) * (n - 1) * n / 6 + (m - 1) * m / 2 + l;
            Mdouble invdet = 1.0 / (+wallObjects_[n].getNormal().X *
                                    (wallObjects_[m].getNormal().Y * wallObjects_[l].getNormal().Z -
                                     wallObjects_[l].getNormal().Y * wallObjects_[m].getNormal().Z)
                                    - wallObjects_[n].getNormal().Y *
                                      (wallObjects_[m].getNormal().X * wallObjects_[l].getNormal().Z -
                                       wallObjects_[m].getNormal().Z * wallObjects_[l].getNormal().X)
                                    + wallObjects_[n].getNormal().Z *
                                      (wallObjects_[m].getNormal().X * wallObjects_[l].getNormal().Y -
                                       wallObjects_[m].getNormal().Y * wallObjects_[l].getNormal().X));
            C_[id] = Vec3D(+(wallObjects_[m].getNormal().Y * wallObjects_[l].getNormal().Z -
                             wallObjects_[l].getNormal().Y * wallObjects_[m].getNormal().Z) *
                           Vec3D::dot(wallObjects_[n].getPosition(), wallObjects_[n].getNormal())
                           - (wallObjects_[n].getNormal().Y * wallObjects_[l].getNormal().Z -
                              wallObjects_[n].getNormal().Z * wallObjects_[l].getNormal().Y) *
                             Vec3D::dot(wallObjects_[m].getPosition(), wallObjects_[m].getNormal())
                           + (wallObjects_[n].getNormal().Y * wallObjects_[m].getNormal().Z -
                              wallObjects_[n].getNormal().Z * wallObjects_[m].getNormal().Y) *
                             Vec3D::dot(wallObjects_[l].getPosition(), wallObjects_[l].getNormal()),
                           -(wallObjects_[m].getNormal().X * wallObjects_[l].getNormal().Z -
                             wallObjects_[m].getNormal().Z * wallObjects_[l].getNormal().X) *
                           Vec3D::dot(wallObjects_[n].getPosition(), wallObjects_[n].getNormal())
                           + (wallObjects_[n].getNormal().X * wallObjects_[l].getNormal().Z -
                              wallObjects_[n].getNormal().Z * wallObjects_[l].getNormal().X) *
                             Vec3D::dot(wallObjects_[m].getPosition(), wallObjects_[m].getNormal())
                           - (wallObjects_[n].getNormal().X * wallObjects_[m].getNormal().Z -
                              wallObjects_[m].getNormal().X * wallObjects_[n].getNormal().Z) *
                             Vec3D::dot(wallObjects_[l].getPosition(), wallObjects_[l].getNormal()),
                           +(wallObjects_[m].getNormal().X * wallObjects_[l].getNormal().Y -
                             wallObjects_[l].getNormal().X * wallObjects_[m].getNormal().Y) *
                           Vec3D::dot(wallObjects_[n].getPosition(), wallObjects_[n].getNormal())
                           - (wallObjects_[n].getNormal().X * wallObjects_[l].getNormal().Y -
                              wallObjects_[l].getNormal().X * wallObjects_[n].getNormal().Y) *
                             Vec3D::dot(wallObjects_[m].getPosition(), wallObjects_[m].getNormal())
                           + (wallObjects_[n].getNormal().X * wallObjects_[m].getNormal().Y -
                              wallObjects_[m].getNormal().X * wallObjects_[n].getNormal().Y) *
                             Vec3D::dot(wallObjects_[l].getPosition(), wallObjects_[l].getNormal())) * invdet;
        }
    }
    
    logger(VERBOSE, "%", *this);
    for (const InfiniteWall& w : wallObjects_)
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
void IntersectionOfWalls::addObject(Vec3D normal, Mdouble position)
{
    logger(WARN, "This function is deprecated, use IntersectionOfWalls::addObject(Vec3D, Vec3D) instead.");
    addObject(normal, position * normal);
}

/*!
 * \param[in] points A vector of 3D-vectors which contains the points between which the polygon is drawn.
 * \param[in] prismAxis A 3D-vector which represents the direction in which the prism is extended infinitely.
 * \details Create an open prism which is a polygon with no connection between the
 * first and last point, and extending infinitely in the other direction, which is
 * defined as PrismAxis. Do this by adding the walls between the consecutive points
 * one by one.
 */
void IntersectionOfWalls::createOpenPrism(std::vector<Vec3D> points, Vec3D prismAxis)
{
    clear();
    //note: use i+1 < points.size() instead of i < points.size()-1, otherwise it creates havoc if point has zero entries.
    for (unsigned int i = 0; i + 1 < points.size(); i++)
        addObject(Vec3D::cross(points[i] - points[i + 1], prismAxis), points[i]);
}

/*!
 * \param[in] points A vector of 3D-vectors which contains the points between which the polygon is drawn.
 * \param[in] prismAxis A 3D-vector which represents the direction in which the prism is extended infinitely.
 * \details Create an open prism which is a polygon and extending infinitely in the other direction, which is
 * defined as PrismAxis. Do this by first creating an open prism and then connect
 * the last and the first point.
 */
void IntersectionOfWalls::createPrism(std::vector<Vec3D> points, Vec3D prismAxis)
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
void IntersectionOfWalls::createOpenPrism(std::vector<Vec3D> points)
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
void IntersectionOfWalls::createPrism(std::vector<Vec3D> points)
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
 * a given BaseParticle and this IntersectionOfWalls. If there is a collision, this
 * function also computes the distance between the BaseParticle and IntersectionOfWalls
 * and the normal of the IntersectionOfWalls at the intersection point. It does
 * this by calling IntersectionOfWalls::getDistanceAndNormal(const Vec3D& , Mdouble , Mdouble&, Vec3D&) const.
 * Since this function should be called before calculating any
 * Particle-Wall interactions, it can also be used to set the normal vector in
 * case of curved walls.
 */
bool IntersectionOfWalls::getDistanceAndNormal(const BaseParticle& p, Mdouble& distance, Vec3D& normal_return) const
{
    //transform coordinates into position-orientation frame
    Vec3D position = p.getPosition() - getPosition();
    getOrientation().rotateBack(position);
    ///\todo do this for all walls
    BaseSpecies* s = getHandler()->getDPMBase()->speciesHandler.getMixedObject(p.getSpecies(), getSpecies());
    if (getDistanceAndNormal(position, p.getRadius() + s->getInteractionDistance(), distance, normal_return))
    {
        getOrientation().rotate(normal_return);
        return true;
    }
    else
    {
        return false;
    }
}

/*!
 * \param[in] position The position of the object there is possible an interaction with.
 * \param[in] wallInteractionRadius The maximum distance between the IntersectionOfWalls
 *  and the input argument position for which there is an interaction.
 * \param[out] distance The distance of the object at position to this wall.
 * \param[out] normal_return If there was an interaction, the normal vector to
 * this wall will be placed here.
 * \return A boolean which says whether or not there was an interaction.
 *
 * \details This function computes whether a particle at the given position (position) and radius (wallInteractionRadius) overlaps with the IntersectionOfWalls.
 * - First, the distances to each InfiniteWall is computed, and the three walls with the largest distances (smallest overlaps) are identified (distance, distance2, distance3, id, id2, id3).
 *
 * - If the largest distance is bigger than the wallInteractionRadius, there is no contact. This, if any distance>wallInteractionRadius is detected we return false.
 * - If the second-largest distance is bigger than the wallInteractionRadius, it is a face contact (contact with a single wall.
 * - Otherwise, we need to determine if it is a edge or vertex contact.
 * In the latter two cases, the function returns true, and the distance and normal vector is returned.
 */
bool IntersectionOfWalls::getDistanceAndNormal(const Vec3D& position, Mdouble wallInteractionRadius, Mdouble& distance,
                                               Vec3D& normal_return) const
{
    if (wallObjects_.empty())
    {
        logger(DEBUG, "Empty IntersectionOfWalls");
        return false;
    }
    
    distance = -1e20;
    Mdouble distance2 = -1e20;
    Mdouble distance3 = -1e20;
    Mdouble distanceCurrent;
    unsigned int id = 0;
    unsigned int id2 = 0;
    unsigned int id3 = 0;
    
    //For each wall we calculate the distance (distanceCurrent). To compute The walls with the minimal overlap
    //The object has to touch each wall  each wall (distanceCurrent) and keep the minimum distance (distance) and wall index (id
    for (unsigned int i = 0; i < wallObjects_.size(); i++)
    {
        // Calculate distance to each wall (distanceCurrent);
        distanceCurrent = wallObjects_[i].getDistance(position);
        // The object has to touch each wall (distanceCurrent >= wallInteractionRadius), otherwise return false (i.e. no contact)
        // This means that for each InfiniteWall in wallObjects_, the particle is either "inside"
        // the wall or touching it. If not, there is no interaction.
        if (distanceCurrent >= wallInteractionRadius)
        {
            return false;
        }
        // Find out which of the InfiniteWalls is interacting with the particle.
        // Keep the minimum distance (distance) and wall index (id)
        // and store up to two walls (id2, id3) and their distances (distance2, distance3),
        // if the possible contact point is near the intersection between id and id2 (and id3)
        if (distanceCurrent > distance)
        {
            if (distance > -wallInteractionRadius) //if distance was set previously
            {
                if (distance2 > -wallInteractionRadius) //if distance2 was set previously
                {
                    distance3 = distance2;
                    id3 = id2;
                }
                distance2 = distance;
                id2 = id;
            }
            distance = distanceCurrent;
            id = i;
        }
        else if (distanceCurrent < -wallInteractionRadius)
        {
            continue;
        }
        else if (distanceCurrent > distance2)
        {
            if (distance2 > -wallInteractionRadius) //if distance2 was set previously
            {
                distance3 = distance2;
                id3 = id2;
            }
            distance2 = distanceCurrent;
            id2 = i;
        }
        else if (distanceCurrent > distance3)
        {
            distance3 = distanceCurrent;
            id3 = i;
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
            //E is the point on wall id2 closest to P
            Vec3D E = position + wallObjects_[id2].getNormal() * distance2;
            //If the distance of D to id2 is positive, the contact is with the intersection
            bool intersection_with_id3 = (wallObjects_[id3].getDistance(E) > 0.0);
            
            if (intersection_with_id2)
            {
                if (intersection_with_id3)
                {
                    //possible contact is with intersection of id,id2,id3
                    //we know id2<id3
                    if (id2 > id3)
                    {
                        auto id = id2;
                        id2 = id3;
                        id3 = id;
                    }
                    unsigned int index =
                            (id < id2) ? ((id3 - 2) * (id3 - 1) * id3 / 6 + (id2 - 1) * id2 / 2 + id) :
                            (id < id3) ? ((id3 - 2) * (id3 - 1) * id3 / 6 + (id - 1) * id / 2 + id2) :
                            ((id - 2) * (id - 1) * id / 6 + (id3 - 1) * id3 / 2 + id2);
                    normal_return = position - C_[index];
                    distance = sqrt(normal_return.getLengthSquared());
                    if (distance <= wallInteractionRadius) //note what if nan?
                    {
                        normal_return /= -distance;
                        return true; //contact with id,id2,id3
                    }
                    else
                    {
                        if (distance == distance)
                            return false;
                    }
                }
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
            if (distance <= wallInteractionRadius) //note what if nan?
            {
                normal_return /= -distance;
                return true; //contact with id,id2,id3
            }
            else
            {
                if (distance == distance)
                    return false;
            }
        }
    }
    //contact is with id
    normal_return = wallObjects_[id].getNormal();
    return true;
}

///*!
// * \param[in] move A reference to a Vec3D that denotes the direction and length
// * it should be moved with.
// * \details A function that moves the InterSectionOfWalls in a certain direction
// * by both moving the walls and all intersections. Note that the directions of the
// * intersections are not moved since they don't change when moving the IntersectionOfWalls
// * as a whole.
// * \todo We should use the position_ and orientation_ of the IntersectionOfWalls;
// * that way, IntersectionOfWalls can be moved with the standard BaseInteractable::move function,
// * getting rid of an anomaly in the code and removing the virtual from the move function. \author weinhartt
// */
//void IntersectionOfWalls::move(const Vec3D& move)
//{
//    BaseInteractable::move(move);
//    for (Vec3D& a : A_)
//    {
//        a += move;
//    }
//    for (Vec3D& c : C_)
//    {
//        c += move;
//    }
//    for (InfiniteWall& o : wallObjects_)
//    {
//        o.move(move);
//    }
//}

/*!
 * \param[in] is The input stream from which the IntersectionOfWalls is read, usually a restart file.
 */
void IntersectionOfWalls::read(std::istream& is)
{
    BaseWall::read(is);
    std::string dummy;
    int n;
    is >> dummy >> n;
    
    Vec3D normal;
    Vec3D position;
    for (int i = 0; i < n; i++)
    {
        is >> dummy;
        if (dummy != "normal")
        {
            Quaternion orientation;
            is >> position >> dummy >> orientation;
            addObject(orientation, position);
        }
        else
        {
            is >> normal >> dummy >> position;
            addObject(normal, position);
        }
    }
}

/*!
 * \param[in] os The output stream where the IntersectionOfWalls must be written
 *  to, usually a restart file.
 */
void IntersectionOfWalls::write(std::ostream& os) const
{
    BaseWall::write(os);
    os << " numIntersectionOfWalls " << wallObjects_.size();
    for (const auto& wallObject : wallObjects_)
    {
        os << " position " << wallObject.getPosition() << " orientation " << wallObject.getOrientation();
    }
}

/*!
 * \return The string "IntersectionOfWalls".
 */
std::string IntersectionOfWalls::getName() const
{
    return "IntersectionOfWalls";
}

void IntersectionOfWalls::writeVTK(VTKContainer& vtk) const
{
    Vec3D max = getHandler()->getDPMBase()->getMax()-getPosition();
    Vec3D min = getHandler()->getDPMBase()->getMin()-getPosition();
    for (auto wall = wallObjects_.begin(); wall != wallObjects_.end(); wall++)
    {
        std::vector<Vec3D> points;
        wall->createVTK(points, min, max);
        for (auto other = wallObjects_.begin(); other != wallObjects_.end(); other++)
        {
            if (other != wall)
            {
                intersectVTK(points, -other->getNormal(), other->getPosition());
            }
        }
        //rotate into real frame
        for (auto& p : points)
        {
            getOrientation().rotate(p);
            p += getPosition();
        }
        wall->addToVTK(points, vtk);
    }
}

//void IntersectionOfWalls::rotate(const Vec3D& rotate)
//{
//    size_t n = 0;
//    for (auto& wall : wallObjects_)
//    {
//        Vec3D p = wall.getPosition() - getPosition();
//        wall.getOrientation().rotateBack(p);
//        wall.setOrientation(wall.getOrientation().updateAngularDisplacement(rotate));
//        wall.getOrientation().rotate(p);
//        wall.setPosition(p + getPosition());
//        setPointsAndLines(n);
//        n++;
//    }
//}

