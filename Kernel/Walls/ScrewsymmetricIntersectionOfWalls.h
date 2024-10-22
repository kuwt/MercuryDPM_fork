//Copyright (c) 2013-2023, The MercuryDPM Developers Team. All rights reserved.
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

#ifndef SCREWSYMMETRICINTERSECTIONOFWALLS_H
#define SCREWSYMMETRICINTERSECTIONOFWALLS_H

#include "IntersectionOfWalls.h"
#include "InteractionHandler.h"
#include "Math/Vector.h"

/*!
 * \brief Use #ScrewsymmetricIntersectionOfWalls to define screwsymmetric walls, such as cylinders, cones, etc
 *
 * \details An #ScrewsymmetricIntersectionOfWalls is equivalent to an #IntersectionOfWalls where the Cartesian coordinate system (x,y,z) is replaced by a cylindrical coordinate system \f$(\hat{r},\theta,\hat{z})\f$. The origin and orientation of the cylindrical coordinate system is defined by the position and orientation of the wall, respectively.
 *
 * In other words, a particle touches an #ScrewsymmetricIntersectionOfWalls, if it touches the #IntersectionOfWalls object in the \f$(r,\theta,z)\f$ coordinate system.
 *
 * Thus, you need to define:
 *  - the position *p* of the wall, which is also the origin of the cylindrical coordinate system
 *  - the orientation *o* of the wall, which is the \f$\hat{z}\f$ direction of the cylindrical coordinate system
 *  - a set of walls in the \f$(\hat{r},\theta,\hat{z})\f$ coordinate system, defined by a normal *n* and position *p*. Only screwsymmetric objects can be defined, thus the \f$\theta \f$ value of the normals has to be zero.
 *
 * ### Example 1
 *
 * Say you want to define a cylindrical wall as in the left image below. If you define the origin and orientation of the cylindrical coordinate system as *(0,0,0)* and *(0,0,1)*, respectively, then the cylinder is a rectangle in the cylindrical coordinate system. Thus, you need to intersect three walls, with normals and position as indicated in the right figure below.
 *
 * \image html ScrewsymmetricWalls.png "A cylindric wall that repels particles"
 *
 * The following code defines such a cylinder:
 * \code
 * ScrewsymmetricIntersectionOfWalls w;
 * w.setSpecies(species);
 * w.setPosition(Vec3D(0,0,0));
 * w.setOrientation(Vec3D(0,0,1));
 * //arguments of addObject are normal and position of the intersected walls
 * w.addObject(Vec3D(-1,0,0), Vec3D(1,0,0));  //Cylindric wall
 * w.addObject(Vec3D(0,0,1), Vec3D(.5,0,-1)); //Bottom wall
 * w.addObject(Vec3D(0,0,-1), Vec3D(.5,0,1)); //Top wall
 * wallHandler.copyAndAddObject(w);
 * \endcode
 *
 * ### Example 2
 *
 * Note, one can also define a cylindric casing that can be filled with particles, see image below.
 *
 * \image html ScrewsymmetricWallsOuter.png "A cylindric casing that can be filled with particles"
 *
 * In this case, you don't have to intersect the walls; instead you need to create three separate walls. A sample code:
 * \code
 * ScrewsymmetricIntersectionOfWalls w;
 * w.setSpecies(species);
 * w.setPosition(Vec3D(0,0,0));
 * w.setOrientation(Vec3D(0,0,1));
 * w.addObject(Vec3D(1,0,0), Vec3D(1,0,0));  //Cylindric wall
 * wallHandler.copyAndAddObject(w);
 *
 * InfiniteWall w1;
 * w1.set(Vec3D(0,0,-1), Vec3D(0,0,-1)); //Bottom wall
 * wallHandler.copyAndAddObject(w1);
 * w1.set(Vec3D(0,0,1), Vec3D(0,0,1)); //Top wall
 * wallHandler.copyAndAddObject(w1);
 * \endcode
 *
 * ### Example 3
 *
 * Say you want a cylindrical casing with an outflow at the base. In this case, you need to define three walls:
 * - The outer cylinder of radius R, height H
 * - The flat top wall
 * - A bottom wall which is a outer cylinder of radius r<R, with a flat top wall at z=0
 *
 * \image html ScrewsymmetricWallsSilo.png "A cylindric casing that can be filled with particles"
 *
 * This can be done as follows:
    \code
    double R = 2;
    double H = 1;
    double r = 1;

    ScrewsymmetricIntersectionOfWalls w;
    w.setSpecies(species);
    w.setPosition(Vec3D(0,0,0));
    w.setOrientation(Vec3D(0,0,1));
    //normal and position of the outer shell in cylindrical coordinates
    w.addObject(Vec3D(1,0,0), Vec3D(R,0,0));
    wallHandler.copyAndAddObject(w);

    InfiniteWall w1;
    w1.set(Vec3D(0,0,1), Vec3D(0,0,H)); //Top wall

    ScrewsymmetricIntersectionOfWalls w0;
    w0.setSpecies(species);
    w0.setPosition(Vec3D(0,0,0));
    w0.setOrientation(Vec3D(0,0,1));
    // bottom wall is an intersection of two walls, an outer cylinder of radius r and a flat top wall at z=0
    w0.addObject(Vec3D(1,0,0), Vec3D(r,0,0));
    w0.addObject(Vec3D(0,0,1), Vec3D(r,0,0));
    wallHandler.copyAndAddObject(w0);
    \endcode
 *

 * For a demonstration on how to use this class, see \ref HourGlass3DDemo.
 */
class ScrewsymmetricIntersectionOfWalls : public IntersectionOfWalls
{
public:
    /*!
     * \brief Default constructor.
     */
    ScrewsymmetricIntersectionOfWalls();
    
    /*!
     * \brief Copy constructor.
     */
    ScrewsymmetricIntersectionOfWalls(const ScrewsymmetricIntersectionOfWalls& p);
    
    /*!
     * \brief Constructor setting values.
     */
    ScrewsymmetricIntersectionOfWalls(Vec3D position, Vec3D normal, std::vector<normalAndPosition> walls,
                                    const ParticleSpecies* species);
    
    /*!
     * \brief Destructor.
     */
    ~ScrewsymmetricIntersectionOfWalls() override;
    
    /*!
     * \brief Copy assignment operator.
     */
    ScrewsymmetricIntersectionOfWalls& operator=(const ScrewsymmetricIntersectionOfWalls& other);
    
    /*!
     * \brief Wall copy method. It calls the copy constructor of this Wall, useful for polymorphism
     */
    ScrewsymmetricIntersectionOfWalls* copy() const final;
    
    /*!
     * \brief Computes the distance from the wall for a given BaseParticle and 
     * returns true if there is a collision. If there is a collision, also 
     * return the normal vector.
     */
    bool getDistanceAndNormal(const BaseParticle& P, Mdouble& distance, Vec3D& normal_return) const final;

    void computeNormalRadialDeltaN(const Vec3D& positionLabFrame, Vec3D& normalVector, Vec3D& radialVector, Mdouble& deltaN) const;

    Mdouble computeDeltaZ(const Vec3D& positionLabFrame, Mdouble h, Mdouble pitch) const;

    /*!
     * \brief reads wall
     */
    void read(std::istream& is) final;
    
    /*!
     * \brief outputs wall
     */
    void write(std::ostream& os) const final;
    
    /*!
     * \brief Returns the name of the object
     */
    std::string getName() const final;
    
    void setAxis(Vec3D a);
    
    /*! converts XYZ limits into RZ limits, to properly limit the VTK plotting area. */
    void convertLimits(Vec3D& min, Vec3D& max) const;
    
    void writeVTK(VTKContainer& vtk) const override;

    void setRightHandedness(bool rightHandedness) {rightHandedness_=rightHandedness;}
    bool getRightHandedness() const {return rightHandedness_;}

    void setPitch(Mdouble pitch) {pitch_=pitch;}
    Mdouble getPitch() const {return pitch_;}

    void setThickness(Mdouble thickness) {thickness_=thickness;}
    Mdouble getThickness() const {return thickness_;}

    void setLength(Mdouble length) {length_=length;}
    Mdouble getLength() const {return length_;}

    ///\brief The right handedness of the screw, i.e. the direction of the screw-blade.
    bool rightHandedness_ = true;

//    void setPointsAndLines(unsigned int n) override {
//        IntersectionOfWalls::setPointsAndLines(n);
//        rMax_ = 0;
//        for (auto p : C_) {
//            rMax_ = std::max(rMax_,p.X);
//        }
//        logger(INFO,"rMax=%",rMax_);
//    }


private:

    /**
     * axial length of one rotation
     */
    Mdouble pitch_ = 1;

    ///\brief The thickness of the screw blade.
    Mdouble thickness_ = 0;

    ///\brief The length of the screw.
    Mdouble length_ = 1;

//    ///\brief The outer screw blade radius.
//    Mdouble rMax_ = 0;
};


#endif
