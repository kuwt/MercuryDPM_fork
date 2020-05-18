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

#ifndef HorizontalScrew_H
#define HorizontalScrew_H

#include "Walls/BaseWall.h"
#include "Math/Vector.h"

/*!
 * \brief This function defines an Archimedes' screw in the z-direction from a (constant)  starting point, a (constant) length L, a (constant) radius r, a (constant) number or revelations N and a (constant) rotation speed (rev/s)
 *
 * \details q is a new coordinate going from 0 to 1 and t is the time, x=xs+r*cos(2*pi*(offset+N*q)), y=ys+r*sin(2*pi*(offset+N*q)), z=zs+q*L
 * \todo IFCD: Can these details about class HorizontalScrew be made more clear? I don't understand them.
 */
class HorizontalScrew : public BaseWall
{
public:
    
    /*!
     * \brief Default constructor: make a screw with default parameters.
     */
    HorizontalScrew();
    
    /*!
     * \brief Copy constructor, copies another HorizontalScrew.
     */
    HorizontalScrew(const HorizontalScrew& other);

    /*!
     * \brief Constructor in which all parameters of the screw are set.
     */
    HorizontalScrew(Vec3D start, Mdouble l, Mdouble minR, Mdouble lowerR, Mdouble diffR, Mdouble n, Mdouble omega, Mdouble thickness, const ParticleSpecies* s);

    /*!
     * \brief Default destructor.
     */
    ~HorizontalScrew();

    /*!
     * \brief Copy this screw and return a pointer to the copy.
     */
    HorizontalScrew* copy() const final;

    /*!
     * \brief Compute the distance from the HorizontalScrew for a given BaseParticle and return if there is a collision. If there is a collision, also return the normal vector of the interaction point.
     */
    bool getDistanceAndNormal(const BaseParticle& P, Mdouble& distance, Vec3D& normal_return) const final;

    /*!
     * \brief Returns the length of the HorizontalScrew.
     */
    Mdouble getLength() const;

    /*!
     * \brief Returns the starting position of the HorizontalScrew.
     */
    Vec3D getStart() const;

    /*!
     * \brief Rotate the HorizontalScrew for a period dt, so that the offset_ changes with omega_*dt.
     */
    void move_time(Mdouble dt);

    /*!
     * \brief Reads a HorizontalScrew from an input stream, for example a restart file.
     */
    void read(std::istream& is) override;

    /*!
     * \brief Writes this HorizontalScrew to an output stream, for example a restart file.
     */
    void write(std::ostream& os) const override;

    /*!
     * \brief Returns the name of the object, here the string "HorizontalScrew".
     */
    std::string getName() const final;

    void writeVTK (VTKContainer& vtk) const override;

    void setBlades(const Mdouble bladeWidth, const Mdouble bladeLength,const std::vector<Mdouble> bladeMounts);

private:
    /*!
     * \brief The centre of the lower end of the screw.
     */
    Vec3D start_;
    /*!
     * \brief The length of the HorizontalScrew.
     */
    Mdouble l_;
    /*!
     * \brief The outer radius of the HorizontalScrew.
     */
    Mdouble minR_;   // min max radius
    Mdouble lowerR_; // max Radius at the base
    Mdouble diffR_; // max radius at the top
    /*!
     * \brief The number of revelations.
     */
    Mdouble n_;
    /*!
     * \brief Rotation speed in rev/s.
     */
    Mdouble omega_;
    /*!
     * \brief The angle that describes how much the HorizontalScrew has turned, going from 0 to 1 for a rotation.
     */
    Mdouble offset_;
    /*!
     * \brief The thickness of the HorizontalScrew.
     */
    Mdouble thickness_;
    /*!
     * \brief The maximum radial width of a blade (in r).
     */
    Mdouble bladeWidth_;
    /*!
     * \brief The length of a blade (in the q-coordinate, which is a linear mapping from start.Z<z<start.Z+l to 0<q<1).
     */
    Mdouble bladeLength_;
    /*!
     * \brief The starting point of a blade (in the q-coordinate, which is a linear mapping from start.Z<z<start.Z+l to 0<q<1)
     */
    std::vector<Mdouble> bladeMounts_;
};

#endif
