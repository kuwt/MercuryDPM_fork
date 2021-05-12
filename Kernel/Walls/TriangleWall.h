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

#ifndef TriangleWall_H
#define TriangleWall_H

#include <vector>
#include "BaseWall.h"
#include "InfiniteWall.h"
#include "Math/Vector.h"

/*!
 * \brief A TriangleWall is convex polygon defined as an intersection of InfiniteWall's.
 * \details It can be defined as the intersection of a set
 * of \ref InfiniteWall's, defined by the normal vector into the wall and a point on the wall:
 * \code
 *   TriangleWall w;
 *   //for each wall, specify a normal and position vector
 *   w.addObject(Vec3D(-1.0,0.0,0.0),Vec3D(1.0,0.0,0.0));
 *   w.addObject(Vec3D(1.0,0.0,0.0),Vec3D(0.0,0.0,0.0));
 *   w.addObject(Vec3D(0.0,-1.0,0.0),Vec3D(0.0,1.0,0.0));
 *   w.addObject(Vec3D(0.0,1.0,0.0),Vec3D(0.0,0.0,0.0));
 *   wallHandler.copyAndAddObject(w);
 * \endcode
 * A particle of radius *r* and position *x* touches an #InfiniteWall with normal *n* and position *p* if
 * \f$p-n\cdot x\leq r\f$
 * (note 'touching particles' also includes particles that are completely enclosed inside the wall).
 * A particle touches an #TriangleWall if it touches all InfiniteWall objects (shown in the image below).
 * <img src="T8_fig2_finitewall.jpg" height="250px">
 *
 * For a demonstration on how to use this class, see \ref T8 and \ref HourGlass3DDemo (shown in the image below).
 * <img src="HourGlass2DDemo.png" height="250px">
 */
class TriangleWall : public BaseWall
{
public:
    
    /*!
     * \brief Default constructor.
     */
    TriangleWall() = default;
    
    /*!
     * \brief Copy constructor.
     */
    TriangleWall(const TriangleWall& other) = default;
    
    /*!
     * \brief Destructor.
     */
    ~TriangleWall() override = default;
    
    /*!
     * \brief Wall copy method. It calls the copy constructor of this Wall, useful for polymorphism
     */
    TriangleWall* copy() const override
    { return new TriangleWall(*this); }
    
    /*!
     * \brief Returns the name of the object, here the string "TriangleWall".
     */
    std::string getName() const override
    { return "TriangleWall"; }
    
    /*!
     * \brief Reads an TriangleWall from an input stream, for example a restart file.
     */
    void read(std::istream& is) override;
    
    /*!
     * \brief Writes an TriangleWall to an output stream, for example a restart file.
     */
    void write(std::ostream& os) const override;
    
    /*!
     * \brief Sets member variables such that the wall represents a triangle with vertices A, B, C
     *  - position_ is set to the center of mass of the wall
     *  - vertexInLabFrame_ is set relative to the position
     *  - updateVertexAndNormal is called to set the remaining variables
     */
    void setVertices(Vec3D A, Vec3D B, Vec3D C);

    std::array<Vec3D,3> getVertices() const {return vertex_;}

    void move(const Vec3D& move) override;
    
    /*!
     * \brief Same as #setVertices(A,B,C), but sets the position explicitly.
     * The position is important when you rotate the wall, as the wall will be rotated around this position.
     */
    void setVertices(Vec3D A, Vec3D B, Vec3D C, Vec3D position);

    void writeVTK(VTKContainer& vtk) const override;

    /*!
     * \brief Returns the position of a vertex
     */
    const Vec3D& getVertex(unsigned i) const {return vertex_[i];}

    /*!
     * \brief Compute the distance from the wall for a given BaseParticle and return if there is a collision.
     * If there is a collision, also return the normal vector.
     */
    bool getDistanceAndNormal(const BaseParticle& p, Mdouble& distance, Vec3D& normal_return) const override;

    void rotate(const Vec3D& angularVelocity) override;
    
    bool isLocal(Vec3D& min, Vec3D& max) const override;

    bool isInsideTriangle(const Vec3D &point) const;

private:
    
    void updateVertexAndNormal();
    
    /*!
     * stores the position of the vertices relative to the position of the wall and rotated into the lab frame;
     */
    std::array<Vec3D, 3> vertexInLabFrame_;
    /*!
     * stores the position of the vertices relative to the position of the wall but not rotated into the lab frame;
     * thus, if the wall rotates, these vertices have to be rotated as well
     */
    std::array<Vec3D, 3> vertex_;
    
    /*!
     * stores the min and max coordinate values of the vertices (needed for hGrid)
     */
    Vec3D vertexMin_;
    Vec3D vertexMax_;
    
    /*!
     * stores the wall normal n in n.x=p
     */
    std::array<Vec3D, 3> edgeNormal_;
    std::array<Vec3D, 3> edge_;
    std::array<double, 3> edgeLength_;

    /*!
     * stores the face normal, not rotated into the lab frame; thus, if the wall rotates, this normal has to be rotated as well
     */
    Vec3D faceNormal_;

//    /*!
//     * number of edges that should be treated as edges (instead of ignored)
//     */
//    unsigned nEdges = 3;
};

#endif
