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

#include <array>
#include "NurbsSurface.h"
#include "NurbsUtils.h"
#include "Logger.h"
#include "Math/ExtendedMath.h"

using namespace NurbsUtils;

NurbsSurface::NurbsSurface() {
    std::vector<double> knotsU = {0,0,1,1};
    std::vector<double> knotsV = {0,0,1,1};
    std::vector<std::vector<Vec3D>> controlPoints = {{{0,0,0},{0,1,0}},{{1,0,0},{1,1,0}}} ;
    std::vector<std::vector<Mdouble>> weights = {{1,1},{1,1}};
    set(knotsU,knotsV,controlPoints,weights);
    //todo think of good default
}

NurbsSurface::NurbsSurface(const std::vector<double>& knotsU, const std::vector<double>& knotsV,
             const std::vector<std::vector<Vec3D>>& controlPoints, const std::vector<std::vector<double>>& weights)
{
    set(knotsU,knotsV,controlPoints,weights);
    logger(INFO,"Created Nurbs surface.");
    logger(INFO,"  %x% knots",knotsU.size(), knotsV.size());
    logger(INFO,"  %x% control points",controlPoints.size(), controlPoints[0].size());
    logger(INFO,"  %x% degrees",degreeU_,degreeV_);
}

void NurbsSurface::set(const std::vector<double>& knotsU, const std::vector<double>& knotsV,
                           const std::vector<std::vector<Vec3D>>& controlPoints, const std::vector<std::vector<double>>& weights)
{
    controlPoints_ = controlPoints;
    weights_ = weights;
    knotsU_ = knotsU;
    knotsV_ = knotsV;

    if ( controlPoints.size()<2 || controlPoints[0].size()<2 ) {
        logger(ERROR,"At least two control points are necessary");
    }
    
    if (std::any_of(controlPoints.begin(), controlPoints.end(),
                    [controlPoints](auto v) { return v.size() != controlPoints[0].size(); }))
    {
        logger(ERROR, "All rows of the control matrix must have the same size");
    }

    if (controlPoints.size() != weights.size() ||
        std::any_of(weights.begin(), weights.end(),
                    [controlPoints](auto v) { return v.size() != controlPoints[0].size(); }))
    {
        logger(ERROR, "Number of control points and weights must match");
    }

    if (knotsU.size() < 2 || knotsV.size() < 2)
    {
        logger(ERROR, "At least two knots are necessary");
    }

    if (!isKnotVectorMonotonic(knotsU) || !isKnotVectorMonotonic(knotsV))
    {
        logger(ERROR, "Knot vector(s) is not monotonic");
    }

    if (knotsU.size() - controlPoints.size()  < 2 || knotsV.size() - controlPoints[0].size() - 1 < 1)
    {
        logger(ERROR, "Degree has to be at least 1");
    }

    //Reset the u/v interval to [0,1]
    const int minKU = knotsU_.front();
    const int maxKU = knotsU_.back();
    for (auto& k : knotsU_) {
        k = (k - minKU) / (maxKU - minKU);
    }
    const int minKV = knotsV_.front();
    const int maxKV = knotsV_.back();
    for (auto& k : knotsV_) {
        k = (k - minKV) / (maxKV - minKV);
    }

    degreeU_ = knotsU.size() - controlPoints.size() - 1;
    degreeV_ = knotsV.size() - controlPoints[0].size() - 1;

    logger(INFO,"Created Nurbs surface.");
    logger(INFO,"  %x% knots",knotsU.size(), knotsV.size());
    logger(INFO,"  %x% control points",controlPoints.size(), controlPoints[0].size());
    logger(INFO,"  %x% degrees",degreeU_,degreeV_);

    //\todo TW
    closedInU_ = false;
    closedInV_ = false;
}

void NurbsSurface::setClosedInU(bool closedInU) {
    closedInU_ = closedInU;
}

void NurbsSurface::setClosedInV(bool closedInV) {
    closedInV_ = closedInV;
}

void NurbsSurface::flipOrientation()
{
    for (auto& cp : controlPoints_) {
        std::reverse(cp.begin(), cp.end());
    }
    for (auto& w : weights_) {
        std::reverse(w.begin(), w.end());
    }
    std::reverse(knotsV_.begin(), knotsV_.end());
    double maxK = knotsV_[0];
    for (auto& k : knotsV_) {
        k = maxK-k;
    }
    // \todo Fix logger messages
    //logger(INFO,"Created Nurbs surface.");
    //logger(INFO,"  %x% knots",knotsU.size(), knotsV.size());
    //logger(INFO,"  %x% control points",controlPoints.size(), controlPoints[0].size());
    //logger(INFO,"  %x% degrees",degreeU_,degreeV_);
}

Vec3D NurbsSurface::evaluate(double u, double v) const {

    Vec3D point = {0,0,0};
    double temp = 0;

    // Find span and non-zero basis functions
    int spanU = findSpan(degreeU_, knotsU_, u);
    int spanV = findSpan(degreeV_, knotsV_, v);
    std::vector<double> Nu, Nv;
    bsplineBasis(degreeU_, spanU, knotsU_, u, Nu);
    bsplineBasis(degreeV_, spanV, knotsV_, v, Nv);
    // make linear combination
    for (int l = 0; l <= degreeV_; ++l)
    {
        for (int k = 0; k <= degreeU_; ++k)
        {
            double weight = Nv[l] * Nu[k] * weights_[spanU - degreeU_ + k][spanV - degreeV_ + l];
            point += weight * controlPoints_[spanU - degreeU_ + k][spanV - degreeV_ + l];
            temp += weight;
        }
    }
    point /= temp;
    return point;
}

/**
 * Find projection onto surface, return distance (and contactPoint)
 */
bool NurbsSurface::getDistance(Vec3D P, double radius, double& distance, Vec3D& normal) const {
    // find the closest control point
    double u;
    double v;
    double minDist2 = constants::inf;
    for (int i=0; i<controlPoints_.size(); ++i) {
        for (int j=0; j<controlPoints_[i].size(); ++j) {
            const double dist2 = Vec3D::getLengthSquared(controlPoints_[i][j] - P);
            if (dist2<minDist2) {
                u = knotsU_[i];
                v = knotsV_[j];
            }
        }
    }
    ///\todo here we should use the convex hull argument to rule out certain contactse quickly
    //Now do Newton-Raphson: Find closest point C(u) to P, see (6.6) in Nurbs book
    //  define f(u)=C'(u).(C(u)-P), find f(u)=0
    //  iterate u <- u-f(u)/f'(u)
    //in 2D:
    //  define r(u,v) = S(u,v)-P
    //         f = r. dSdu, g = r. dSdv
    //  iterate u <- u + du, v <- v + dv
    //          J.[du;dv] = k
    //          J=[fu fv;gu gv], k=-[f;g]
    std::array<std::array<Vec3D,3>,3> S;
    const double tol = 2*std::numeric_limits<double>::epsilon();
    const double tolSquared = mathsFunc::square<double>(tol);
    Vec3D r;
    double r1;
    for (int i=0;i<15; ++i) {
        //TW this algorithm does not compute contacts with the boundary
        if (u<0) {
            u = closedInU_?(u-static_cast<int>(u)):0;
        } else if (u>1) {
            u = closedInU_?(u-static_cast<int>(u)):1;
        }
        if (v<0) {
            v = closedInV_?(v-static_cast<int>(u)):0;
        } else if (v>1) {
            v = closedInV_?(v-static_cast<int>(v)):1;
        }
        evaluateDerivatives(u,v,S);
        r = S[0][0] - P;
        r1 = r.getLengthSquared();
        if (r1<tolSquared) /*first convergence criterium: contact point on surface*/ {
            distance = 0;
            normal = r;
            return true;
        }
        const double r2 = mathsFunc::square<double>(Vec3D::dot(r, S[1][0]))/(r1*S[1][0].getLengthSquared());
        const double r3 = mathsFunc::square<double>(Vec3D::dot(r, S[0][1]))/(r1*S[0][1].getLengthSquared());
        if (std::fabs(r2)<tolSquared && std::fabs(r3)<tolSquared) /*second convergence criterium: zero cosine*/ {
            //you should exit the function here
            const bool inWall = Vec3D::dot(r,Vec3D::cross(S[1][0],S[0][1]))>=0;
            if (inWall || r1<radius*radius) {
                logger(VERBOSE,"contact found at % after % iterations",S[0][0],i);
                distance = inWall ? -sqrt(r1) : sqrt(r1);
                normal = r / distance;
                return true;
            } else {
                logger(VERBOSE,"contact found, but too distant");
                return false;
            }
        }
        const double f = Vec3D::dot(r, S[1][0]);
        const double g = Vec3D::dot(r, S[0][1]);
        const double a = S[1][0].getLengthSquared()   + Vec3D::dot(r, S[2][0]);
        const double b = Vec3D::dot(S[1][0], S[0][1]) + Vec3D::dot(r, S[1][1]);
        const double d = S[0][1].getLengthSquared()   + Vec3D::dot(r, S[0][2]);
        const double det = a * d - b * b;
        const double du = (b * g - d * f) / det;
        const double dv = (b * f - a * g) / det;
        const double r4 = du*du*S[0][1].getLengthSquared()+dv*dv*S[1][0].getLengthSquared();
        if (r4<tolSquared) /*third convergence criterium: parameters fixed*/ {
            const bool inWall = Vec3D::dot(r,Vec3D::cross(S[1][0],S[0][1]))>=0;
            if (inWall || r1<radius*radius) {
                logger(VERBOSE,"parameters fixed, contact at % after % iterations",S[0][0],i);
                distance = inWall ? -sqrt(r1) : sqrt(r1);
                normal = r / distance;
                return true;
            } else {
                logger(VERBOSE,"parameters fixed, but too distant");
                return false;
            }
        }
        u += du;
        v += dv;
    }
    /* iteration fails */
    logger(WARN,"P=%: Number of allowed iterations exceeded; this should never be reached",P);
    const bool inWall = Vec3D::dot(r,Vec3D::cross(S[1][0],S[0][1]))>=0;
    if (inWall || r1<radius*radius) {
        logger(VERBOSE,"contact found at %",S[0][0]);
        distance = inWall ? -sqrt(r1) : sqrt(r1);
        normal = r / distance;
        return true;
    } else {
        logger(VERBOSE,"contact found, but too distant");
        return false;
    }
}


/**
Evaluate derivatives of a rational NURBS curve
@param[in] u Parameter to evaluate the derivatives at.
@param[in] knots Knot vector of the curve.
@param[in] controlPoints Control points of the curve.
@param[in] weights Weights corresponding to each control point.
@param[in] nDers Number of times to differentiate.
@param[in, out] curveDers Derivatives of the curve at u.
E.g. curveDers[n] is the nth derivative at u, where n is between 0 and nDers-1.
*/
void NurbsSurface::evaluateDerivatives(double u, double v, std::array<std::array<Vec3D,3>,3>& S) const
{
    std::array<std::array<Vec3D,3>,3> A;
    std::array<std::array<double,3>,3> w;
    // Find span and non-zero basis functions
    int spanU = findSpan(degreeU_, knotsU_, u);
    int spanV = findSpan(degreeV_, knotsV_, v);
    std::vector<std::vector<double>> Nu, Nv;
    bsplineDerBasis(degreeU_, spanU, knotsU_, u, 2, Nu);
    bsplineDerBasis(degreeV_, spanV, knotsV_, v, 2, Nv);
    // make linear combination (eq 4.21 in Nurbs book)
    for (int i=0; i<3; ++i) {
        for (int j = 0; j < 3; ++j) {
            A[i][j].setZero();
            w[i][j] = 0;
            for (int k = 0; k <= degreeU_; ++k) {
                for (int l = 0; l <= degreeV_; ++l) {
                    double weight = Nu[i][k] * Nv[j][l] * weights_[spanU - degreeU_ + k][spanV - degreeV_ + l];
                    A[i][j] += weight * controlPoints_[spanU - degreeU_ + k][spanV - degreeV_ + l];
                    w[i][j] += weight;
                }
            }
        }
    }
    S[0][0] = A[0][0]/w[0][0];
    S[1][0] = (A[1][0]-w[1][0]*S[0][0])/w[0][0];
    S[0][1] = (A[0][1]-w[0][1]*S[0][0])/w[0][0];
    S[1][1] = (A[1][1]-w[1][1]*S[0][0]-w[1][0]*S[0][1]-w[0][1]*S[1][0])/w[0][0];
    S[2][0] = (A[2][0]-2*w[1][0]*S[1][0]-w[2][0]*S[0][0])/w[0][0];
    S[0][2] = (A[0][2]-2*w[0][1]*S[0][1]-w[0][2]*S[0][0])/w[0][0];
//    std::cout << S[0][0] << std::endl;
//    std::cout << S[1][0] << std::endl;
//    std::cout << S[0][1] << std::endl;
//    std::cout << S[2][0] << std::endl;
//    std::cout << S[1][1] << std::endl;
//    std::cout << S[0][2] << std::endl;
}


/*!
 * \details Adds all elements of the vector to an output stream.
 * NB: this is a global function and a friend of the Vec3D class!
 * \param[in] os    the output stream,
 * \param[in] a     The vector of interest
 * \return          the output stream with vector elements added
 */
std::ostream& operator<<(std::ostream& os, const NurbsSurface& a)
{
    os << "mu " << a.knotsU_.size() << ' ';
    os << "mv " << a.knotsV_.size() << ' ';
    os << "nu " << a.controlPoints_.size() << ' ';
    os << "nv " << a.controlPoints_[0].size() << ' ';
    os << "knotsU ";
    for (const auto k : a.knotsU_) os << k << ' ';
    os << "knotsV ";
    for (const auto k : a.knotsV_) os << k << ' ';
    os << "controlPoints ";
    for (const auto cp0 : a.controlPoints_) for (const auto cp : cp0) os << cp << ' ';
    os << "weights ";
    for (const auto w0 : a.weights_) for (const auto w : w0) os << w << ' ';
    return os;
}

/*!
 * \details Reads all elements of a given vector from an input stream.
 * NB: this is a global function and a friend of the Vec3D class!
 * \param[in,out] is    the input stream
 * \param[in,out] a     the vector to be read in
 * \return              the input stream from which the vector elements were read
 */
std::istream& operator>>(std::istream& is, NurbsSurface& a)
{
    std::string dummy;
    unsigned mu, mv, nu, nv;
    is >> dummy >> mu;
    is >> dummy >> mv;
    is >> dummy >> nu;
    is >> dummy >> nv;

    is >> dummy;
    std::vector<double> knotsU;
    knotsU.resize(mu);
    for (auto& k : knotsU) is >> k;

    is >> dummy;
    std::vector<double> knotsV;
    knotsV.resize(mv);
    for (auto& k : knotsV) is >> k;

    is >> dummy;
    std::vector<std::vector<Vec3D>> controlPoints;
    controlPoints.resize(nu);
    for (auto& cp0 : controlPoints) {
        cp0.resize(nv);
        for (auto& cp : cp0) is >> cp;
    }

    is >> dummy;
    std::vector<std::vector<double>> weights;
    weights.resize(nu);
    for (auto& w0 : weights) {
        w0.resize(nv);
        for (auto& w : w0) is >> w;
    }

    a.set(knotsU,knotsV,controlPoints,weights);
    return is;
}
