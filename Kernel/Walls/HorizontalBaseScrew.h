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

#ifndef HorizontalBaseScrew_H
#define HorizontalBaseScrew_H

#include "Walls/IntersectionOfWalls.h"
#include "InteractionHandler.h"
#include "Math/Vector.h"

/*!
 * \brief A HorizontalBaseScrew is a copy of AxisymmetricIntersectionOfWalls, with an additional, angle-dependent component.
 *
 * \details All changes w.r.t. AxisymmetricIntersectionOfWalls are encapsulated in //BEGIN CHANGE and //BEGIN CHANGE
 */
class HorizontalBaseScrew : public IntersectionOfWalls
{
public:
    /*!
     * \brief Default constructor.
     */
    HorizontalBaseScrew();

    /*!
     * \brief Copy constructor.
     */
    HorizontalBaseScrew(const HorizontalBaseScrew& p);

    /*!
     * \brief Constructor setting values.
     */
    HorizontalBaseScrew(Vec3D position, Vec3D orientation, std::vector<normalAndPosition> walls, const ParticleSpecies* species, Mdouble sinA2Max, Mdouble timeMin);

    /*!
     * \brief Destructor.
     */
    ~HorizontalBaseScrew();

    /*!
     * \brief Copy assignment operator.
     */
    HorizontalBaseScrew& operator=(const HorizontalBaseScrew& other);

    /*!
     * \brief Wall copy method. It calls the copy constructor of this Wall, useful for polymorphism
     */
    HorizontalBaseScrew* copy() const final;

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

    //BEGIN CHANGE
    Mdouble sinA2Max_;
    Mdouble timeMin_;
    //END CHANGE
};


#endif
