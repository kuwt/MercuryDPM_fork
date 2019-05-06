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

#ifndef ROTATINGINTERSECTIONOFWALLS_H
#define ROTATINGINTERSECTIONOFWALLS_H

#include <vector>
#include <Walls/BaseWall.h>
#include <Walls/InfiniteWall.h>
#include <Math/Vector.h>

/*!
 * \brief A RotatingIntersectionOfWalls is convex polygon defined as an intersection of InfiniteWall's.
 * \details It can be defined as the intersection of a set
 * of #InfiniteWalls, defined by the normal vector into the wall and a point on the wall:
 * \code
 *   RotatingIntersectionOfWalls w;
 *   //for each wall, specify a normal and position vector
 *   w.addObject(Vec3D(-1.0,0.0,0.0),Vec3D(1.0,0.0,0.0));
 *   w.addObject(Vec3D(1.0,0.0,0.0),Vec3D(0.0,0.0,0.0));
 *   w.addObject(Vec3D(0.0,-1.0,0.0),Vec3D(0.0,1.0,0.0));
 *   w.addObject(Vec3D(0.0,1.0,0.0),Vec3D(0.0,0.0,0.0));
 *   wallHandler.copyAndAddObject(w);
 * \endcode
 * A particle of radius *r* and position *x* touches an #InfiniteWall with normal *n* and position *p* if \f$p-n\cdot x\leq r\f$
 * (note 'touching particles' also includes particles that are completely enclosed inside the wall).
 * A particle touches an #RotatingIntersectionOfWalls if it touches all InfiniteWall objects (shown in the image below).
 * <img src="T8_fig2_finitewall.jpg" height="250px">
 *
 * For a demonstration on how to use this class, see \ref T8 and \ref HourGlass3DDemo (shown in the image below).
 * <img src="HourGlass2DDemo.png" height="250px">
 */
class RotatingIntersectionOfWalls : public BaseWall
{
public:
    struct normalAndPosition {
        Vec3D normal;
        Vec3D position;
    };

    /*!
     * \brief Default constructor.
     */
    RotatingIntersectionOfWalls();

    /*!
     * \brief Copy constructor.
     */
    RotatingIntersectionOfWalls(const RotatingIntersectionOfWalls& other);

    /*!
     * \brief Constructor setting values.
     */
    RotatingIntersectionOfWalls(std::vector<normalAndPosition> walls, const ParticleSpecies* species);

    /*!
     * \brief Destructor.
     */
    virtual ~RotatingIntersectionOfWalls();

    /*!
     * Copy assignment operator.
     */
    RotatingIntersectionOfWalls& operator=(const RotatingIntersectionOfWalls& other);

    /*!
     * \brief Wall copy method. It calls the copy constructor of this Wall, useful for polymorphism
     */
    RotatingIntersectionOfWalls* copy() const override;

    /*!
     * \brief Removes all parts of the walls.
     * \deprecated Please don't use any clear() anymore, it will be gone soon.
     */

    /// sets species of subwalls as well
    void setSpecies(const ParticleSpecies* species);

    void setHandler(WallHandler* wallHandler) override;

    /*!
     * \brief Adds a wall to the set of infinite walls, given an outward normal vector s.t. normal*x=normal*point
     */
    void addObject(Vec3D normal, Vec3D point);

    /*!
     * \brief Adds a wall to the set of finite walls, given an outward normal vector s. t. normal*x=position
     * \deprecated Don't use this function, instead use the function addObject(Vec3D, Vec3D).
     */
    MERCURY_DEPRECATED
    void addObject(Vec3D normal, Mdouble position);

    /*!
     * \brief Creates an open prism which is a polygon between the points, except the first and last point,  and extends infinitely in the PrismAxis direction.
     */
    void createOpenPrism(std::vector<Vec3D> points, Vec3D prismAxis);

    /*!
     * \brief Creates an open prism which is a polygon between the points and extends infinitely in the PrismAxis direction.
     */
    void createPrism(std::vector<Vec3D> points, Vec3D prismAxis);

    /*!
     * \brief Creates an open prism which is a polygon between the points, except the first and last point, and extends infinitely in the direction perpendicular to the first and second wall.
     */
    void createOpenPrism(std::vector<Vec3D> points);

    /*!
     * \brief Creates an open prism which is a polygon between the points and extends infinitely in the direction perpendicular to the first and second wall.
     */
    void createPrism(std::vector<Vec3D> points);

    /*!
     * \brief Compute the distance from the wall for a given BaseParticle and return if there is a collision. If there is a collision, also return the normal vector.
     */
    bool getDistanceAndNormal(const BaseParticle& p, Mdouble& distance, Vec3D& normal_return) const override;

    /*!
     * \brief Compute the distance from the wall for a given BaseParticle and return if there is an interaction. If there is an interaction, also return the normal vector.
     */
    bool getDistanceAndNormal(const Vec3D& postition, Mdouble wallInteractionRadius, Mdouble& distance, Vec3D& normal_return) const;

    /*!
     * \brief Move the RotatingIntersectionOfWalls to a new position, which is a Vec3D from the old position.
     */
    void move(const Vec3D& move) override;

    // HACK - BEGIN
    void moveAlongZ(double displacement);
    void rotateAroundZ(double angle);
    // HACK - END

    /*!
     * \brief Reads an RotatingIntersectionOfWalls from an input stream, for example a restart file.
     */
    void read(std::istream& is) override;

    /*!
     * \brief Writes an RotatingIntersectionOfWalls to an output stream, for example a restart file.
     */
    void write(std::ostream& os) const override;

    /*!
     * \brief Returns the name of the object, here the string "RotatingIntersectionOfWalls".
     */
    std::string getName() const override;

    void writeVTK (VTKContainer& vtk) const override;

protected:
    /*!
     * \brief The wall "segments"/directions that together make up the finite wall.
     * \details An intersection of walls exists of a number of infinite walls that
     * are cut of at the intersection points. These InfiniteWall are saved in this
     * vector called iWObjects_.
     */
    std::vector<InfiniteWall> wallObjects_;

private:
    /*!
     * \brief A vector that stores a point for each intersecting line between two different InfiniteWall.
     * \details A[n*(n-1)/2+m] is a point on the intersecting line between walls
     *  m and n, m<n.
     */
    std::vector<Vec3D> A_;

    /*!
     * \brief A vector that stores the direction of the intersecting lines between two different InfiniteWall.
     * \details AB[n*(n-1)/2+m] is the direction of the intersecting line between
     *  walls m and n, m<n.
     */
    std::vector<Vec3D> AB_;

    /*!
     * \brief A vector that stores the intersection point of three different InfiniteWall.
     * \details C[(n-2)*(n-1)*n/6+(m-1)*m/2+l] is a point intersecting walls
     * l, m and n, l<m<n
     */
    std::vector<Vec3D> C_;
};

#endif
