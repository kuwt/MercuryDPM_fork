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

#ifndef MERCURY_NURBSSURFACE_H
#define MERCURY_NURBSSURFACE_H

#include <vector>
#include <array>
#include "Math/Vector.h"

class NurbsSurface
{
public:
    /**
      Create a non-rational NurbsSurface
    */
    NurbsSurface();

    /**
      Create a non-rational NurbsSurface
      @param degreeU Degree of the surface in u-direction
      @param degreeV Degree of the surface in v-direction
      @param knotsU Knot vector in u-direction
      @param knotsV Knot vector in v-direction
      @param controlPoints 2D vector of control points
    */
    NurbsSurface(const std::vector<double>& knotsU, const std::vector<double>& knotsV,
                 const std::vector<std::vector<Vec3D>>& controlPoints,
                 const std::vector<std::vector<double>>& weights);
    
    /**
    Evaluate point on a nonrational NURBS surface
    @param[in] u Parameter to evaluate the surface at.
    @param[in] v Parameter to evaluate the surface at.
    @param[in] degreeU Degree of the given surface in u-direction.
    @param[in] degreeV Degree of the given surface in v-direction.
    @param[in] knotsU Knot vector of the surface in u-direction.
    @param[in] knotsV Knot vector of the surface in v-direction.
    @param[in] controlPoints Control points of the surface.
    @param[in, out] point Resulting point on the surface at (u, v).
    */
    Vec3D evaluate(double u, double v) const;

    void evaluateDerivatives(double u, double v, std::array<std::array<Vec3D,3>,3>& S) const;

    /**
     * Find projection onto surface, return distance (and contactPoint)
     */
    bool getDistance(Vec3D P, double radius, double& distance, Vec3D& normal) const;

    /**
      Create a non-rational NurbsSurface
      @param degreeU Degree of the surface in u-direction
      @param degreeV Degree of the surface in v-direction
      @param knotsU Knot vector in u-direction
      @param knotsV Knot vector in v-direction
      @param controlPoints 2D vector of control points
    */
    void set(const std::vector<double>& knotsU, const std::vector<double>& knotsV,
             const std::vector<std::vector<Vec3D>>& controlPoints, const std::vector<std::vector<double>>& weights);

    void setClosedInU(bool closedInU);

    void setClosedInV(bool closedInV);

    void flipOrientation();

    void closestPoint(Vec3D position, double& u, double& v) const;

    void splitSurface (int spanU, int spanV) {

    }

    /*!
     * \brief Adds elements to an output stream
     */
    friend std::ostream& operator<<(std::ostream& os, const NurbsSurface& a);

    /*!
     * \brief Adds elements to an input stream
     */
    friend std::istream& operator>>(std::istream& is, NurbsSurface& a);


private:
    ///mu knots
    std::vector<double> knotsU_;
    ///mv knots
    std::vector<double> knotsV_;
    ///nu x nv control points
    std::vector<std::vector<Vec3D>> controlPoints_;
    ///nu x nv weights
    std::vector<std::vector<double>> weights_;
    ///degree pu = mu-nu-1, pv = mv-nv-1
    unsigned int degreeU_, degreeV_;
    ///make it a periodic system
    bool closedInU_, closedInV_;
};

#endif //MERCURY_NURBSSURFACE_H
