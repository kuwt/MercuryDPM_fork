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

#ifndef NURBSWALL_H
#define NURBSWALL_H

#include "Nurbs/NurbsSurface.h"
#include "BaseWall.h"

/*!
 * \brief This function defines a wall via a NurbsSurface
 */
class NurbsWall : public BaseWall
{
public:

    /*!
     * \brief Default constructor: make a wall with default parameters.
     */
    NurbsWall();

    /*!
     * \brief Copy constructor, copies another wall.
     */
    NurbsWall(const NurbsWall& other);

    /*!
     * \brief Constructor in which all parameters of the wall are set.
     */
    NurbsWall(const NurbsSurface& nurbsSurface);
    
    /*!
     * \brief Default destructor.
     */
    ~NurbsWall();

    /*!
     * \brief Copy this wall and return a pointer to the copy.
     */
    NurbsWall* copy() const final;

    /*!
     * \brief Defines a standard wall, given an outward normal vector s.t. normal*x=normal*point for all x of the wall.
     */
    void set(const NurbsSurface& nurbsSurface);

    /*!
     * \brief Compute the distance from the Screw for a given BaseParticle and return if there is a collision. If there is a collision, also return the normal vector of the interaction point.
     */
    bool getDistanceAndNormal(const BaseParticle& P, Mdouble& distance, Vec3D& normal_return) const final;

    /*!
     * \brief Reads this wall from an input stream, for example a restart file.
     */
    void read(std::istream& is) override;

    /*!
     * \brief Writes this wall to an output stream, for example a restart file.
     */
    void write(std::ostream& os) const override;

    /*!
     * \brief Returns the name of the object, here the string "Screw".
     */
    std::string getName() const final;

    void writeVTK (VTKContainer &vtk) const override;

private:
    /*!
     * \brief The centre of the lower end of the screw.
     */
    NurbsSurface nurbsSurface_;
};

#endif
