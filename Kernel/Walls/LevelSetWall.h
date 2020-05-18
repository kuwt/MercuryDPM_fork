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

#ifndef LevelSetWall_H
#define LevelSetWall_H

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

class LevelSetWall final : public BaseWall
{
public:
    
    // a set of predefined shapes for which level set values can be set easily
    enum class Shape
    {
        Sphere,
        Cube,
        Diamond,
        FourSided,
        Cylinder
    };
    
    // Constructor; currently only allows predefined shapes
    LevelSetWall(Shape s, double radius, ParticleSpecies* sp = nullptr);
    
    /*!
     * \brief Default destructor.
     */
    ~LevelSetWall() override;
    
    /*!
     * \brief Wall copy method. It calls the copy constructor of this Wall, useful for polymorphism
     */
    LevelSetWall* copy() const override;
    
    using BaseWall::move;

//    /*!
//     * \brief Returns the distance of the wall to the particle.
//     */
//    Mdouble getDistance(Vec3D position) const;
    
    /*!
     * \brief Compute the distance from the wall for a given BaseParticle and return if there is a collision. If there is a collision, also return the normal vector.
     */
    bool getDistanceAndNormal(const BaseParticle& p, Mdouble& distance, Vec3D& normal_return) const override;
    
    /*!
     * \brief Reads LevelSetWall from a restart file.
     */
    void read(std::istream& is) override;
    
    /*!
     * \brief Returns the name of the object, in this case the string "LevelSetWall".
     */
    std::string getName() const override;
    
    /*!
     * Adds the vtk wall representation to the VTK container
     */
    void writeVTK(VTKContainer& vtk) const override;
    
    /*!
    * \param n how many points should be used to interpolate
    * \param radiusContact how much should be added to the radius
    */
    void writeToFile(int n, double radiusContact) const;

private:
    
    // Evaluates the distance and normal direction to the surface defined by the level set. Needed for contact detection.
    // If distance is bigger than radius_+radiusContact, no normal & distance is returned because no contact is possible.
    // Else, it returns the overlap (level set value) and normal direction (-gradient).
    bool getDistanceAndNormalLabCoordinates(Vec3D position, Mdouble interactionRadius, Mdouble& distance,
                                            Vec3D& normal) const;
    
    void setShapeSphere();
    
    void setShapeCube();
    
    //Better interpolated than a square as the edges of the level set align with the mesh.
    void setShapeDiamond();

    void setShapeCylinder();

    void setShapeFourSided();
    
    void createVTKSphere();
    
    void createVTKCube();
    
    void createVTKDiamond();
    
    void createVTK();
    
    
    // N determines the number of level set values (2N+1)^3;
    // \todo template the level set with N.
    static const int N = 10;
    
    // discrete set of level-set values; levelSet_[i,j,k] is the value of the level set at (x,y,z)=(i,j,k)/N*radius_.
    double levelSet_[2 * N + 1][2 * N + 1][2 * N + 1];
    
    // the radius of a sphere that envelopes the object
    double radius_ = 1;
    
    VTKContainer vtkLabFrame_;
};

#endif
