//Copyright (c) 2013-2018, The MercuryDPM Developers Team. All rights reserved.
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

#ifndef AXISYMMETRICINTERSECTIONOFWALLS_H
#define AXISYMMETRICINTERSECTIONOFWALLS_H

#include "IntersectionOfWalls.h"
#include "InteractionHandler.h"
#include "Math/Vector.h"

/*!
 * \brief A AxisymmetricIntersectionOfWalls is an axisymmetric wall, defined by
 * rotating a twodimensional IntersectionOfWalls around a symmetry axis.
 * \details It is defined by defining an axis position *p* and orientation *o*, around which the wall is axisymmetric, and an IntersectionOfWalls object *w*.
 * To determine if a particle of touches an #AxisymmetricIntersectionOfWalls,
 * the particle's location is determined is a cylindrical coordinate system \f$(r,\theta,z)\f$,
 * where *r* is the radial distance from the axis, *z* is the axial distance along the axis, and \f$\theta\f$ the angular position around the axis (which is ignored as the object is axisymmetric).
 * A particle touches an #AxisymmetricIntersectionOfWalls if it touches the #IntersectionOfWalls *w* in the \f$(r,\theta,z)\f$ coordinate system.
 *
 * The code below defines a cylindrical outer wall of radius *r* with an axis position *p* and orientation *o* and contact properties *s*:
 * \code
 * ParticleSpecies s = setSpecies(speciesHandler.getObject(0));
 * Vec3D p = Vec3D(0.0, 0.0, 0.0);
 * Vec3D o = Vec3D(0.0, 0.0, 1.0);
 * Mdouble r = 1.0;
 *
 * AxisymmetricIntersectionOfWalls w;
 * w.setSpecies(s);
 * w.setPosition(p);
 * w.setOrientation(o);
 * w.addObject(Vec3D(1.0, 0.0, 0.0), Vec3D(r, 0.0, 0.0));
 * wallHandler.copyAndAddObject(w);
 * \endcode
 *
 * For a demonstration on how to use this class, see \ref HourGlass3DDemo (shown in the image below).
 * <img src="HourGlass3DDemo.png" height="250px">
 */
class AxisymmetricIntersectionOfWalls : public IntersectionOfWalls
{
public:
    /*!
     * \brief Default constructor.
     */
    AxisymmetricIntersectionOfWalls();

    /*!
     * \brief Copy constructor.
     */
    AxisymmetricIntersectionOfWalls(const AxisymmetricIntersectionOfWalls& p);

    /*!
     * \brief Constructor setting values.
     */
    AxisymmetricIntersectionOfWalls(Vec3D position, Vec3D normal, std::vector<normalAndPosition> walls, const ParticleSpecies* species);

    /*!
     * \brief Destructor.
     */
    ~AxisymmetricIntersectionOfWalls() override;

    /*!
     * \brief Copy assignment operator.
     */
    AxisymmetricIntersectionOfWalls& operator=(const AxisymmetricIntersectionOfWalls& other);

    /*!
     * \brief Wall copy method. It calls the copy constructor of this Wall, useful for polymorphism
     */
    AxisymmetricIntersectionOfWalls* copy() const final;

    /*!
     * \brief Computes the distance from the wall for a given BaseParticle and 
     * returns true if there is a collision. If there is a collision, also 
     * return the normal vector.
     */
    bool getDistanceAndNormal(const BaseParticle& P, Mdouble& distance, Vec3D& normal_return) const final;

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
    void convertLimits (Vec3D& min, Vec3D& max) const;

    void writeVTK (VTKContainer& vtk) const override;

};


#endif
