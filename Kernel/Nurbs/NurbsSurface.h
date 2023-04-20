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

#ifndef MERCURY_NURBSSURFACE_H
#define MERCURY_NURBSSURFACE_H

#include <vector>
#include <array>
#include "Math/Vector.h"

class NurbsSurface
{
public:
    /**
     * Create a NurbsSurface
     */
    NurbsSurface();
    
    /**
     * Create a NurbsSurface
     * @param knotsU Knot vector in u-direction
     * @param knotsV Knot vector in v-direction
     * @param controlPoints 2D vector of control points
     * @param weights 2D vector of weights
     */
    NurbsSurface(const std::vector<double>& knotsU, const std::vector<double>& knotsV,
                 const std::vector<std::vector<Vec3D>>& controlPoints,
                 const std::vector<std::vector<double>>& weights);
    
    /*!
     * \brief Create a NurbsSurface. This will create uniform knot vectors in u and v (clamped by default).
     * @param controlPoints 2D vector of control points
     * @param weights 2D vector of weights
     * @param degreeU Degree in u
     * @param degreeV Degree in v
     * @param clampedAtStartU Clamp knots in u at start
     * @param clampedAtEndU Clamp knots in u at end
     * @param clampedAtStartV Clamp knots in v at start
     * @param clampedAtEndV Clamp knots in v at end
     */
    NurbsSurface(const std::vector<std::vector<Vec3D>>& controlPoints,
                 const std::vector<std::vector<Mdouble>>& weights,
                 unsigned int degreeU, unsigned int degreeV,
                 bool clampedAtStartU = true, bool clampedAtEndU = true,
                 bool clampedAtStartV = true, bool clampedAtEndV = true);
    
    /**
     * Evaluate point on a NURBS surface
     * @param u Parameter to evaluate the surface at
     * @param v Parameter to evaluate the surface at
     * @return Resulting point on the surface at (u, v)
     */
    Vec3D evaluate(double u, double v) const;

    void evaluateDerivatives(double u, double v, std::array<std::array<Vec3D,3>,3>& S) const;

    /**
     * Find projection onto surface, return distance (and contactPoint)
     */
    bool getDistance(Vec3D P, double radius, double& distance, Vec3D& normal) const;
    
    /**
     * Create a NurbsSurface
     * @param knotsU Knot vector in u-direction
     * @param knotsV Knot vector in v-direction
     * @param controlPoints 2D vector of control points
     * @param weights 2D vector of weights
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
    
    /*!
     * @return Lowest allowed value for u
     */
    double getLowerBoundU() const{
        return knotsU_[degreeU_];
    }
    /*!
     * @return Highest allowed value for u
     */
    double getUpperBoundU() const{
        return knotsU_[knotsU_.size() - degreeU_ - 1]; // same as knotsU_[controlPoints.size()];
    }
    /*!
     * @return Lowest allowed value for v
     */
    double getLowerBoundV() const{
        return knotsV_[degreeV_];
    }
    /*!
     * @return Highest allowed value for v
     */
    double getUpperBoundV() const{
        return knotsV_[knotsV_.size() - degreeV_ - 1]; // same as knotsV_[controlPoints[0].size()];
    }
    
    const std::vector<std::vector<Vec3D>>& getControlPoints() const
    { return controlPoints_; }
    
    const std::vector<std::vector<Mdouble>>& getWeights() const
    { return weights_; }
    
    const std::vector<Mdouble>& getKnotsU() const
    { return knotsU_; }
    
    const std::vector<Mdouble>& getKnotsV() const
    { return knotsV_; }
    
    /*!
     * \brief This will make the surface repeat itself and ensure continuity over periodic boundaries
     * \details The first and last control point are assumed to be close to or exactly on the periodic boundaries.
     * Degree-1 amount of control points are copied from both sides to the other. On both sides this results in degree
     * amount of control points "overlapping", which ensures continuity.
     * When both ends of the original knot vector weren't uniform, this will change the shape of the surface a bit, however the inner most surface remains intact.
     */
    void makePeriodicContinuousInU();
    
    /*!
     * \brief This will make the surface repeat itself and ensure continuity over periodic boundaries
     * \details The first and last control point are assumed to be close to or exactly on the periodic boundaries.
     * Degree-1 amount of control points are copied from both sides to the other. On both sides this results in degree
     * amount of control points "overlapping", which ensures continuity.
     * When both ends of the original knot vector weren't uniform, this will change the shape of the surface a bit, however the inner most surface remains intact.
     */
    void makePeriodicContinuousInV();
    
    /*!
     * \brief This will make the surface close around on itself and ensure continuity
     * \details The first and last control point are assumed to already be in the same position, then degree-1 control
     * points are copied from the start to the end. Resulting in degree amount of control points are overlapping, which ensures continuity.
     * When both ends of the original knot vector weren't uniform, this will change the shape of the surface a bit, however the inner most surface remains intact.
     */
    void makeClosedInU();
    
    /*!
     * \brief This will make the surface close around on itself and ensure continuity
     * \details The first and last control point are assumed to already be in the same position, then degree-1 control
     * points are copied from the start to the end. Resulting in degree amount of control points are overlapping, which ensures continuity.
     * When both ends of the original knot vector weren't uniform, this will change the shape of the surface a bit, however the inner most surface remains intact.
     */
    void makeClosedInV();
    
    /*!
     * \brief Unclamps the knot vector by changing the control points, weights and knots.
     * @param inU Whether to unclamp the knots in u-direction or the knots in v-direction
     * @param atStart Whether to unclamp the knots at the start or the end
     */
    void unclampKnots(bool inU, bool atStart);
    
    void moveControlPoint(unsigned int indexU, unsigned int indexV, Vec3D dP, bool includingClosedOrPeriodic);

private:
    ///mu knots
    std::vector<Mdouble> knotsU_;
    ///mv knots
    std::vector<Mdouble> knotsV_;
    ///nu x nv control points
    std::vector<std::vector<Vec3D>> controlPoints_;
    ///nu x nv weights
    std::vector<std::vector<Mdouble>> weights_;
    ///degree pu = mu-nu-1, pv = mv-nv-1
    unsigned int degreeU_, degreeV_;
    ///make it a periodic system
    bool closedInU_, closedInV_;
    // For when using periodic boundaries and the surface is extended to make it continuous
    bool periodicInU_, periodicInV_;
    
    // Storing starting points for a quick start of the getDistance method
    std::vector<Vec3D> startingPoints_;
    std::vector<Mdouble> startingKnotsU_;
    std::vector<Mdouble> startingKnotsV_;
    
    /*!
     * \brief Copies control points from the start and adds to the end and vice versa.
     *        The first and last control points are ignored, as they are used the indicate the distance to be shifted by.
     * @param numStartToEnd Amount to copy from start and add to end
     * @param numEndToStart Amount to copy from end and insert before start
     * @param forceBothEndsUniform When only copying to one side, whether or not the other end should remain untouched
     */
    void wrapAroundInU(unsigned int numStartToEnd, unsigned int numEndToStart, bool forceBothEndsUniform = false);
    
    /*!
     * \brief Copies control points from the start and adds to the end and vice versa.
     *        The first and last control points are ignored, as they are used the indicate the distance to be shifted by.
     * @param numStartToEnd Amount to copy from start and add to end
     * @param numEndToStart Amount to copy from end and insert before start
     * @param forceBothEndsUniform When only copying to one side, whether or not the other end should remain untouched
     */
    void wrapAroundInV(unsigned int numStartToEnd, unsigned int numEndToStart, bool forceBothEndsUniform = false);
};

#endif //MERCURY_NURBSSURFACE_H
