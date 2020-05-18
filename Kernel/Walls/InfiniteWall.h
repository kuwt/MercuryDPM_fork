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

///This is a class defining walls. It defines the 
///interaction of regular walls and periodic walls
///with particles as defined in Particle
///Modifications:

#ifndef INFINITEWALL_H
#define INFINITEWALL_H

#include "BaseWall.h"
#include "Math/Vector.h"

/*!
 * \brief A infinite wall fills the half-space {point: (position_-point)*normal_<=0}. 
 * \details Thus, the surface of the wall is a plane through position position_
 * with normal_ the outward unit normal vector of the wall 
 * (pointing away from the particles, into the wall).
 * Please note that this wall is infinite and straight.
 * 
 * A particle touches an infinite wall if (position_-point)*normal_<=radius.
 */

class InfiniteWall final : public BaseWall
{
public:
    
    /*!
     * \brief Default constructor, the normal is infinitely long.
     */
    InfiniteWall();
    
    /*!
     * \brief Copy constructor, copy the given wall.
     */
    InfiniteWall(const InfiniteWall& w);
    
    /*!
     * \brief Constructor setting species.
     */
    explicit InfiniteWall(const ParticleSpecies* species);
    
    /*!
     * \brief Constructor setting values.
     */
    InfiniteWall(Vec3D normal, Vec3D point, const ParticleSpecies* species);
    
    /*!
     * \brief Constructor setting values if 3 coordinates are given
     */
    InfiniteWall(Vec3D PointA, Vec3D PointB, Vec3D PointC, const ParticleSpecies* species);
    
    /*!
     * \brief Default destructor.
     */
    ~InfiniteWall() override;
    
    /*!
     * \brief Wall copy method. It calls the copy constructor of this Wall, useful for polymorphism
     */
    InfiniteWall* copy() const override;
    
    /*!
     * \brief Defines a standard wall, given an outward normal vector s.t. normal*x=normal*point for all x of the wall.
     */
    void set(Vec3D normal, Vec3D point);
    
    /*!
     * \brief Changes the normal of the InfiniteWall.
     */
    void setNormal(Vec3D normal);
    
    /*!
     * \brief Defines a standard wall by computing normal*position = point and using the overloaded function set(Vec3D, vec3D).
     * \deprecated In Mercury 2, the user will have to use the new interface, 
     *             namely set(Vec3D, Vec3D).
     */
    MERCURY_DEPRECATED
    void set(Vec3D normal, Mdouble position);
    
    using BaseWall::move;
    
    /*!
     * \brief Returns the distance of the wall to the particle.
     */
    Mdouble getDistance(Vec3D position) const;
    
    /*!
     * \brief Compute the distance from the wall for a given BaseParticle and return if there is a collision. If there is a collision, also return the normal vector.
     */
    bool getDistanceAndNormal(const BaseParticle& p, Mdouble& distance, Vec3D& normal_return) const override;
    
    /*!
    * \brief Compute the distance from the wall for a given BaseParticle and return if there is a collision. If there is a collision, also return the normal vector.
    */
    bool getDistanceNormalOverlapSuperquadric(const SuperQuadricParticle& p, Mdouble& distance, Vec3D& normal_return,
                                              Mdouble& overlap) const override;
    
    /*!
     * \brief Reads InfiniteWall from a restart file.
     */
    void read(std::istream& is) override;
    
    /*!
     * \brief Reads InfiniteWall from an old-style restart file.
     */
    void oldRead(std::istream& is);
    
    /*!
     * \brief Writes the InfiniteWall to an output stream, usually a restart file.
     */
    //void write(std::ostream& os) const override;
    
    /*!
     * \brief Returns the name of the object, in this case the string "InfiniteWall".
     */
    std::string getName() const override;
    
    /*!
     * \brief Access function for normal.
     */
    Vec3D getNormal() const;
    
    /*!
     * Returns all intersection points of the infinite wall with the domain boundary.
     */
    void createVTK(std::vector<Vec3D>& myPoints) const;
    
    /*!
     * Same as createVTK(), but with a self-defined domain size (useful for plotting AxisymmetricWall's).
     */
    void createVTK(std::vector<Vec3D>& myPoints, Vec3D max, Vec3D min) const;
    
    /*!
     * Adds the vtk wall representation to the VTK container
     */
    void writeVTK(VTKContainer& vtk) const override;
    
    Vec3D
    getFurthestPointSuperQuadric(const Vec3D& normalBodyFixed, const Vec3D& axes, Mdouble eps1, Mdouble eps2) const override;
};

#endif
