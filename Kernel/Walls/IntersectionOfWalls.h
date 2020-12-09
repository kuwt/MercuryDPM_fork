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

#ifndef INTERSECTIONOFWALLS_H
#define INTERSECTIONOFWALLS_H

#include <vector>
#include "BaseWall.h"
#include "InfiniteWall.h"
#include "Math/Vector.h"

/*!
 * \brief A IntersectionOfWalls is convex polygon defined as an intersection of InfiniteWall's.
 * \details It can be defined as the intersection of a set
 * of \ref InfiniteWall's, defined by the normal vector into the wall and a point on
 * the wall.
 * For example, the following gives a cube |x|<1 and |y|<1:
 * \code
 *   IntersectionOfWalls w;
 *   //for each wall, specify a normal and position vector
 *   w.addObject(Vec3D(-1, 0, 0),Vec3D(1, 0, 0));
 *   w.addObject(Vec3D(1, 0, 0),Vec3D(0, 0, 0));
 *   w.addObject(Vec3D(0, -1, 0),Vec3D(0, 1, 0));
 *   w.addObject(Vec3D(0, 1, 0),Vec3D(0, 0, 0));
 *   wallHandler.copyAndAddObject(w);
 * \endcode
 * A particle of radius *r* and position *x* touches an #InfiniteWall with normal *n* and position *p* if
 * \f$p-n\cdot x\leq r\f$
 * (note 'touching particles' also includes particles that are completely enclosed inside the wall).
 * A particle touches an #IntersectionOfWalls if it touches all InfiniteWall objects (shown in the image below).
 * <img src="T8_fig2_finitewall.jpg" height="250px">
 *
 * For a demonstration on how to use this class, see \ref T8 and \ref HourGlass3DDemo (shown in the image below).
 * <img src="HourGlass2DDemo.png" height="250px">
 */
class IntersectionOfWalls : public BaseWall
{
public:
    struct normalAndPosition
    {
        Vec3D normal;
        Vec3D position;
    };
    
    /*!
     * \brief Default constructor.
     */
    IntersectionOfWalls();
    
    /*!
     * \brief Copy constructor.
     */
    IntersectionOfWalls(const IntersectionOfWalls& other);
    
    /*!
     * \brief Constructor setting values.
     */
    IntersectionOfWalls(const std::vector<normalAndPosition>& walls, const ParticleSpecies* species);
    
    /*!
     * \brief Destructor.
     */
    ~IntersectionOfWalls() override;
    
    /*!
     * Copy assignment operator.
     */
    IntersectionOfWalls& operator=(const IntersectionOfWalls& other);
    
    /*!
     * \brief Wall copy method. It calls the copy constructor of this Wall, useful for polymorphism
     */
    IntersectionOfWalls* copy() const override;
    
    /*!
     * \brief Removes all parts of the walls.
     */
    void clear();
    
    /// sets species of subwalls as well
    void setSpecies(const ParticleSpecies* species);
    
    void setHandler(WallHandler* wallHandler) override;
    
    /*! 
     * \brief Returns the number of objects 
     */
    unsigned int getNumberOfObjects();
    
    /*!
     * \brief Adds a wall to the set of infinite walls, given a normal vector pointing into the wall (i.e. out of the simulation domain), going through the point, 
     * so that normal*x=normal*point
     */
    void addObject(Vec3D normal, Vec3D point);
    
    void addObject(Quaternion orientation, Vec3D position);
    
    void add3PointObject(Vec3D PointA, Vec3D PointB, Vec3D PointC);
    
    void setPointsAndLines(unsigned int n);
    
    /*!
     * \brief constructs a tetrahedron for an STL file input
     */
    void addTetraSTL(Vec3D PointA, Vec3D PointB, Vec3D PointC, Vec3D WallNormal, Mdouble Thickness, int wallidentifier);
    
    /*!
     * \brief constructs a tetrahedron from 3 input coordinates
     */
    void addTetra(const Vec3D& PointA, const Vec3D& PointB, const Vec3D& PointC, Mdouble& Thickness);
    
    
    void addPlate(const Vec3D& PointA, const Vec3D& PointB, const Vec3D& PointC, const Vec3D& WallNormal,
                  const Mdouble& Thickness, int wallidentifier);
    
    /*!
     * \brief Adds a wall to the set of finite walls, given an normal vector pointing into the wall (i.e. out of the flow domain), 
     * to give a plane defined by normal*x=position
     * \deprecated Don't use this function, instead use the function addObject(Vec3D, Vec3D).
     */
    MERCURY_DEPRECATED
    void addObject(Vec3D normal, Mdouble position);
    
    /*!
     * \brief Creates an open prism which is a polygon between the points, except the first and last point,
     * and extends infinitely in the PrismAxis direction. Note that if you view from inside of your geometry,
     * the shape formed by points has to be convex, otherwise it will not create the wall correctly.
     */
    void createOpenPrism(std::vector<Vec3D> points, Vec3D prismAxis);
    
    /*!
     * \brief Creates an open prism which is a polygon between the points and extends infinitely in the PrismAxis
     * direction. Note that if you view from inside of your geometry,
     * the shape formed by points has to be convex, otherwise it will not create the wall correctly.
     */
    void createPrism(std::vector<Vec3D> points, Vec3D prismAxis);
    
    /*!
     * \brief Creates an open prism which is a polygon between the points, except the first and last point, and extends
     * infinitely in the direction perpendicular to the first and second wall. Note that if you view from inside of your geometry,
     * the shape formed by points has to be convex, otherwise it will not create the wall correctly.
     */
    void createOpenPrism(std::vector<Vec3D> points);
    
    /*!
     * \brief Creates an open prism which is a polygon between the points and extends infinitely in the direction
     * perpendicular to the first and second wall. Note that if you view from inside of your geometry,
     * the shape formed by points has to be convex, otherwise it will not create the wall correctly.
     */
    void createPrism(std::vector<Vec3D> points);
    
    /*!
     * \brief Compute the distance from the wall for a given BaseParticle and return if there is a collision.
     * If there is a collision, also return the normal vector.
     */
    bool getDistanceAndNormal(const BaseParticle& p, Mdouble& distance, Vec3D& normal_return) const override;
    
    /*!
     * \brief Compute the distance from the wall for a given BaseParticle and return if there is an interaction.
     * If there is an interaction, also return the normal vector.
     */
    bool getDistanceAndNormal(const Vec3D& position, Mdouble wallInteractionRadius, Mdouble& distance,
                              Vec3D& normal_return) const;
    
    /*!
     * \brief Move the IntersectionOfWalls to a new position, which is a Vec3D from the old position.
     * \todo currently the IntersectionOfWall is special for moving and rotating; should we remove that specialty?
     */
//    void move(const Vec3D& move) override;
    
    /*!
     * \brief Reads an IntersectionOfWalls from an input stream, for example a restart file.
     */
    void read(std::istream& is) override;
    
    /*!
     * \brief Writes an IntersectionOfWalls to an output stream, for example a restart file.
     */
    void write(std::ostream& os) const override;
    
    /*!
     * \brief Returns the name of the object, here the string "IntersectionOfWalls".
     */
    std::string getName() const override;
    
    void writeVTK(VTKContainer& vtk) const override;

//    void rotate(const Vec3D& rotate) override;

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

protected:
    /*!
     * \brief A vector that stores the intersection point of three different InfiniteWall.
     * \details C[(n-2)*(n-1)*n/6+(m-1)*m/2+l] is a point intersecting walls 
     * l, m and n, l<m<n
     */
    std::vector<Vec3D> C_;
};

#endif
