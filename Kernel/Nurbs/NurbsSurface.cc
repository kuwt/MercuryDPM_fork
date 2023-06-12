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
}

NurbsSurface::NurbsSurface(const std::vector<std::vector<Vec3D>>& controlPoints,
                           const std::vector<std::vector<Mdouble>>& weights,
                           unsigned int degreeU, unsigned int degreeV,
                           bool clampedAtStartU, bool clampedAtEndU, bool clampedAtStartV, bool clampedAtEndV)
{
    std::vector<Mdouble> knotsU = createUniformKnotVector(controlPoints.size(), degreeU, clampedAtStartU, clampedAtEndU);
    std::vector<Mdouble> knotsV = createUniformKnotVector(controlPoints[0].size(), degreeV, clampedAtStartV, clampedAtEndV);
    
    set(knotsU, knotsV, controlPoints, weights);
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

    // Reset the u/v interval to [0,1]
    normalizeKnotVector(knotsU_);
    normalizeKnotVector(knotsV_);

    degreeU_ = knotsU.size() - controlPoints.size() - 1;
    degreeV_ = knotsV.size() - controlPoints[0].size() - 1;

    logger(INFO,"Created Nurbs surface.");
    logger(INFO,"  %x% knots",knotsU.size(), knotsV.size());
    logger(INFO,"  %x% control points",controlPoints.size(), controlPoints[0].size());
    logger(INFO,"  %x% degrees",degreeU_,degreeV_);

    //\todo TW
    closedInU_ = false;
    closedInV_ = false;
    periodicInU_ = false;
    periodicInV_ = false;
    
    // This simply evaluates the positions for a bunch of u's and v's and stores them.
    // As it is now, the u's and v's are simply taken uniform.
    // The more points the better (sort of). Too little points can cause a possible convergence to a non-global minimum,
    // or possible no convergence at all.
    // For smooth surfaces this already helps a lot. For non-smooth surfaces it most likely doesn't suffice.
    //\todo JWB The convergence still fails from time to time
    //\todo JWB These aren't updated when the surface changes shape during simulation (moving control points)
    unsigned  nu = knotsU_.size() * 3;
    unsigned  nv = knotsV_.size() * 3;
    startingPoints_.clear();
    startingKnotsU_.clear();
    startingKnotsV_.clear();
    for (double i = 0; i <= nu; i++) {
        double u = getLowerBoundU() + (getUpperBoundU() - getLowerBoundU()) * i / nu;
        for (double j = 0; j <= nv; j++) {
            double v = getLowerBoundV() + (getUpperBoundV() - getLowerBoundV()) * j / nv;
            Vec3D p = evaluate(u, v);
            startingPoints_.push_back(p);
            startingKnotsU_.push_back(u);
            startingKnotsV_.push_back(v);
        }
    }
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
    //JWB Known reasons why convergence might fail (or worse it doesn't fail but finds a wrong point):
    //    1. Multiple same value knots (together with a low degree).
    //    2. Multiple control points at the exact same position.
    //    3. Really high value weights, or big difference in weights.
    //    4. Too little starting points to start off iteration.
    //    5. In general possibly a too low degree.
    // In some cases there is a clear discontinuity, which then is the obvious reason for failure.
    // In other cases this is not really clear, so the reason isn't quite clear.
    // One last reason the iteration might fail is when it actually got quite close to the closest point, but the
    // tolerance is set too tight. This failure does not cause huge problems though, as when the iteration fails the
    // distance is calculated anyway. The tolerance value is pretty much picked from thin air, so improving it is fine.
    
    // Find the closest starting point
    double u, v, minDist2 = constants::inf;
    for (int i = 0; i < startingPoints_.size(); i++) {
        const double dist2 = Vec3D::getLengthSquared(startingPoints_[i] - P);
        if (dist2 < minDist2) {
            u = startingKnotsU_[i];
            v = startingKnotsV_[i];
            minDist2 = dist2;
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
    //JWB The tolerance originally was pretty much grabbed from thin air. Multiplying it by 100 allowed for additional
    // convergence (of the third convergence criterion: parameters fixed) for a bunch of positions.
    const double tol = 2*std::numeric_limits<double>::epsilon() * 100;
    const double tolSquared = mathsFunc::square<double>(tol);
    Vec3D r;
    double r1;
    double previousU, previousV;
    for (int i = 0; i < 15; ++i) {
        evaluateDerivatives(u,v,S);
        r = S[0][0] - P;
        r1 = r.getLengthSquared();
        if (r1 < tolSquared) /*first convergence criterion: contact point on surface*/ {
            distance = 0;
            normal = r;
            return true;
        }
        const double r2 = mathsFunc::square<double>(Vec3D::dot(r, S[1][0]))/(r1*S[1][0].getLengthSquared());
        const double r3 = mathsFunc::square<double>(Vec3D::dot(r, S[0][1]))/(r1*S[0][1].getLengthSquared());
        if (std::fabs(r2) < tolSquared && std::fabs(r3) < tolSquared) /*second convergence criterion: zero cosine*/ {
            //you should exit the function here
            if (r1 < radius * radius) {
                logger(VERBOSE,"contact found at % after % iterations", S[0][0], i);
                distance = sqrt(r1);
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
        
        previousU = u;
        previousV = v;
        u += du;
        v += dv;
        
        // Keeping within bounds
        Mdouble lbU = getLowerBoundU();
        Mdouble ubU = getUpperBoundU();
        Mdouble lbV = getLowerBoundV();
        Mdouble ubV = getUpperBoundV();
        if (u < lbU) {
            u = closedInU_ ? ubU - fmod(lbU - u, ubU - lbU) : lbU;
        }
        else if (u > ubU) {
            u = closedInU_ ? lbU + fmod(u - ubU, ubU - lbU) : ubU;
        }
        if (v < lbV) {
            v = closedInV_ ? ubV - fmod(lbV - v, ubV - lbV) : lbV;
        }
        else if (v > ubV) {
            v = closedInV_ ? lbV + fmod(v - ubV, ubV - lbV) : ubV;
        }
    
        const double r4 = ((u - previousU) * S[1][0] + (v - previousV) * S[0][1]).getLengthSquared();
        if (r4 < tolSquared) /*third convergence criterion: parameters fixed*/ {
            if (r1 < radius * radius) {
                logger(VERBOSE,"parameters fixed, contact at % after % iterations", S[0][0], i);
                distance = sqrt(r1);
                normal = r / distance;
                return true;
            } else {
                logger(VERBOSE,"parameters fixed, but too distant");
                return false;
            }
        }
    }
    /* iteration fails */
    //JWB See comments at the start of this method for possible reasons for failure and how it might be solved.
    // For now a warning is given and whatever was found so far is returned (which is basically rubbish).
    logger(WARN,"P=%: Number of allowed iterations exceeded; this should never be reached", P);
    if (r1 < radius * radius) {
        logger(VERBOSE,"contact found at %", S[0][0]);
        distance = sqrt(r1);
        normal = r / distance;
        return true;
    } else {
        logger(VERBOSE,"contact found, but too distant");
        return false;
    }
}


/**
 * Evaluate derivatives of a NURBS curve
 * @param[in] u Parameter to evaluate derivatives at
 * @param[in] v Parameter to evaluate derivatives at
 * @param[out] S Contains position, first order derivatives and second order derivatives
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
    // See equations around eq 4.21 in Nurbs book
    // [0][0] =      e.g. S[0][0] = S      w[0][0] = w
    // [1][0] = u    e.g. S[1][0] = S_u    w[1][0] = w_u     d/du
    // [0][1] = v    e.g. S[0][1] = S_v    w[0][1] = w_v     d/dv
    // [1][1] = uv   e.g. S[1][1] = S_uv   w[1][1] = w_uv    d^2/dudv
    // [2][0] = uu   e.g. S[2][0] = S_uu   w[2][0] = w_uu    d^2/du^2
    // [0][2] = vv   e.g. S[0][2] = S_vv   w[0][2] = w_vv    d^2/dv^2
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
    for (const auto& cp0 : a.controlPoints_) for (const auto cp : cp0) os << cp << ' ';
    os << "weights ";
    for (const auto& w0 : a.weights_) for (const auto w : w0) os << w << ' ';
    os << "closedInUV " << a.closedInU_ << ' ' << a.closedInV_;
    os << " periodicInUV " << a.periodicInU_ << ' ' << a.periodicInV_;
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
    
    // After setting, since in that method the defaults are set to false.
    // Also, check if dummy variable exist, for backwards compatibility.
    is >> dummy;
    if (dummy == "closedInUV")
    {
        is >> a.closedInU_;
        is >> a.closedInV_;
    }
    is >> dummy;
    if (dummy == "periodicInUV")
    {
        is >> a.periodicInU_;
        is >> a.periodicInV_;
    }
    
    return is;
}

void NurbsSurface::makePeriodicContinuousInU()
{
    // Add degree-1 amount of control points to the front and back.
    wrapAroundInU(degreeU_ - 1, degreeU_ - 1);
    periodicInU_ = true;
}

void NurbsSurface::makePeriodicContinuousInV()
{
    // Add degree-1 amount of control points to the front and back.
    wrapAroundInV(degreeV_ - 1, degreeV_ - 1);
    periodicInV_ = true;
}

void NurbsSurface::makeClosedInU()
{
    // Add degree-1 amount of control points to the back.
    wrapAroundInU(degreeU_ - 1, 0, true);
    setClosedInU(true);
}

void NurbsSurface::makeClosedInV()
{
    // Add degree-1 amount of control points to the back.
    wrapAroundInV(degreeV_ - 1, 0, true);
    setClosedInV(true);
}

void NurbsSurface::wrapAroundInU(unsigned int numStartToEnd, unsigned int numEndToStart, bool forceBothEndsUniform)
{
    // This method copies the given number of control points from the start to the end and vice versa (not counting the
    // first and last control point).
    // The first and last control point are used to know the amount that the control points to be copied have to be
    // shifted by. In case of closing a surface (circle like shape) the shifted distance is simply 0.
    // The knot vector is updated in such a way that only the start and end of the shape might be influenced (both ends
    // should be uniform), however the inner shape remains intact.
    
    // A vector of offsets, which is the difference between the first and last row of control points.
    // These are the values the copied control points need to be shifted by.
    std::vector<Vec3D> offsets;
    offsets.reserve(controlPoints_[0].size());
    for (int j = 0; j < controlPoints_[0].size(); j++)
    {
        offsets.push_back(controlPoints_.back()[j] - controlPoints_.front()[j]);
    }
    
    // Temporarily store "ghost" control points and weights to be added in front
    std::vector<std::vector<Vec3D>> frontControlPoints;
    std::vector<std::vector<Mdouble>> frontWeights;
    // Get the numEndToStart amount, in order, ignoring the last one
    for (unsigned int i = numEndToStart; i > 0; i--)
    {
        std::vector<Vec3D> tmpCP = controlPoints_[controlPoints_.size() - 1 - i];
        for (int j = 0; j < tmpCP.size(); j++)
        {
            tmpCP[j] -= offsets[j];
        }
        frontControlPoints.push_back(tmpCP);
        frontWeights.push_back(weights_[weights_.size() - 1 - i]);
    }
    
    // Add "ghost" control points and weights at the back
    // Get the numStartToEnd amount, in order, ignoring the first one
    for (int i = 1; i <= numStartToEnd; i++)
    {
        std::vector<Vec3D> tmpCP = controlPoints_[i];
        for (int j = 0; j < tmpCP.size(); j++)
        {
            tmpCP[j] += offsets[j];
        }
        controlPoints_.push_back(tmpCP);
        weights_.push_back(weights_[i]);
    }
    
    // Now actually add the "ghost" control points and weights in front
    controlPoints_.insert(controlPoints_.begin(), frontControlPoints.begin(), frontControlPoints.end());
    weights_.insert(weights_.begin(), frontWeights.begin(), frontWeights.end());
    
    extendKnotVector(knotsU_, degreeU_, numEndToStart, numStartToEnd, forceBothEndsUniform);
}

void NurbsSurface::wrapAroundInV(unsigned int numStartToEnd, unsigned int numEndToStart, bool forceBothEndsUniform)
{
    // This method copies the given number of control points from the start to the end and vice versa (not counting the
    // first and last control point).
    // The first and last control point are used to know the amount that the control points to be copied have to be
    // shifted by. In case of closing a surface (circle like shape) the shifted distance is simply 0.
    // The knot vector is updated in such a way that only the start and end of the shape might be influenced (both ends
    // should be uniform), however the inner shape remains intact.
    
    for (int i = 0; i < controlPoints_.size(); i++)
    {
        // Current offset
        Vec3D offset = controlPoints_[i].back() - controlPoints_[i].front();
        
        // Temporarily store "ghost" control points and weights to be added in front
        std::vector<Vec3D> frontControlPoints;
        std::vector<Mdouble> frontWeights;
        // Get the numEndToStart amount, in order, ignoring the last one
        for (unsigned int j = numEndToStart; j > 0; j--)
        {
            frontControlPoints.push_back(controlPoints_[i][controlPoints_[i].size() - 1 - j] - offset);
            frontWeights.push_back(weights_[i][weights_[i].size() - 1 - j]);
        }
        
        // Add "ghost" control points and weights at the back
        // Get the numStartToEnd amount, in order, ignoring the first one
        for (int j = 1; j <= numStartToEnd; j++)
        {
            controlPoints_[i].push_back(controlPoints_[i][j] + offset);
            weights_[i].push_back(weights_[i][j]);
        }
        
        // Now actually add the "ghost" control points and weights in front
        controlPoints_[i].insert(controlPoints_[i].begin(), frontControlPoints.begin(), frontControlPoints.end());
        weights_[i].insert(weights_[i].begin(), frontWeights.begin(), frontWeights.end());
    }
    
    extendKnotVector(knotsV_, degreeV_, numEndToStart, numStartToEnd, forceBothEndsUniform);
}

void NurbsSurface::unclampKnots(bool inU, bool atStart)
{
    // Algorithm from Nurbs book A12.1 extended to surfaces
    
    std::vector<Mdouble>& knots = inU ? knotsU_ : knotsV_;
    unsigned int n = inU ? controlPoints_.size() - 1 : controlPoints_[0].size() - 1;
    unsigned int p = inU ? degreeU_ : degreeV_;
    
    if (atStart)
    {
        // Unclamp at left end
        for (int i = 0; i <= p-2; i++)
        {
            knots[p - i - 1] = knots[p - i] - (knots[n - i + 1] - knots[n - i]);
            int k = p - 1;
            for (int j = i; j >= 0; j--)
            {
                Mdouble alpha = (knots[p] - knots[k]) / (knots[p + j + 1] - knots[k]);
                k--;
                
                // Update changed for each control point and weight in other direction
                // Only difference between inU or not is the indices are swapped: [j][l] -> [l][j]
                if (inU)
                {
                    for (int l = 0; l < controlPoints_[0].size(); l++)
                    {
                        controlPoints_[j][l] = (controlPoints_[j][l] - alpha * controlPoints_[j+1][l]) / (1.0 - alpha);
                        weights_[j][l] = (weights_[j][l] - alpha * weights_[j+1][l]) / (1.0 - alpha);
                    }
                }
                else
                {
                    for (int l = 0; l < controlPoints_.size(); l++)
                    {
                        controlPoints_[l][j] = (controlPoints_[l][j] - alpha * controlPoints_[l][j+1]) / (1.0 - alpha);
                        weights_[l][j] = (weights_[l][j] - alpha * weights_[l][j+1]) / (1.0 - alpha);
                    }
                }
            }
        }
        // Set first knot
        knots[0] = knots[1] - (knots[n - p + 2] - knots[n - p + 1]);
    }
    else
    {
        // Unclamp at right end
        for (int i = 0; i <= p-2; i++)
        {
            knots[n + i + 2] = knots[n + i + 1] + (knots[p + i + 1] - knots[p + i]);
            for (int j = 1; j >= 0; j--)
            {
                Mdouble alpha = (knots[n + 1] - knots[n - j]) / (knots[n - j + i + 2] - knots[n - j]);
                
                // Update changed for each control point and weight in other direction
                // Only difference between inU or not is the indices are swapped: [j][l] -> [l][j]
                if (inU)
                {
                    for (int l = 0; l < controlPoints_[0].size(); l++)
                    {
                        controlPoints_[n - j][l] = (controlPoints_[n - j][l] - (1.0 - alpha) * controlPoints_[n - j -1][l]) / alpha;
                        weights_[n - j][l] = (weights_[n - j][l] - (1.0 - alpha) * weights_[n - j -1][l]) / alpha;
                    }
                }
                else
                {
                    for (int l = 0; l < controlPoints_.size(); l++)
                    {
                        controlPoints_[l][n - j] = (controlPoints_[l][n - j] - (1.0 - alpha) * controlPoints_[l][n - j -1]) / alpha;
                        weights_[l][n - j] = (weights_[l][n - j] - (1.0 - alpha) * weights_[l][n - j -1]) / alpha;
                    }
                }
            }
        }
        // Set last knot
        knots[n + p + 1] = knots[n + p] + (knots[2*p] - knots[2*p - 1]);
    }
    
    normalizeKnotVector(knots);
}

void NurbsSurface::moveControlPoint(unsigned int indexU, unsigned int indexV, Vec3D dP, bool includingClosedOrPeriodic)
{
    // When the surface is closed or periodic, there might be other control points which should also move.
    // To that end, store all the indices in u and v in two vectors, starting with the current one.
    // Then later each combination of indices of u and v is moved.
    std::vector<unsigned int> moveIndexU, moveIndexV;
    moveIndexU.push_back(indexU);
    moveIndexV.push_back(indexV);
    
    if (includingClosedOrPeriodic)
    {
        // Assuming either closed or periodic in u/v, not both.
        unsigned int lowestIndexU, lowestIndexV, shiftU, shiftV;
    
        // For closed, degree-1 control points were added to the end
        if (closedInU_)
        {
            lowestIndexU = degreeU_ - 1;
            shiftU = controlPoints_.size() - degreeU_;
        }
        if (closedInV_)
        {
            lowestIndexV = degreeV_ - 1;
            shiftV = controlPoints_[0].size() - degreeV_;
        }
    
        // For periodic, degree-1 control points were added in front and to the end
        if (periodicInU_)
        {
            lowestIndexU = 2 * (degreeU_ - 1);
            shiftU = controlPoints_.size() - 2 * degreeU_ + 1;
        }
        if (periodicInV_)
        {
            lowestIndexV = 2 * (degreeV_ - 1);
            shiftV = controlPoints_[0].size() - 2 * degreeV_ + 1;
        }
    
        if (closedInU_ || periodicInU_)
        {
            // Note, no else if, because technically it might be possible for e.g. the middle control point to be the degree-1
            // from the starting index as well as the end index. Meaning that on both sides other control points need moving.
            if (indexU <= lowestIndexU)
                moveIndexU.push_back(indexU + shiftU);
            if (indexU >= shiftU)
                moveIndexU.push_back(indexU - shiftU);
        }
        if (closedInV_ || periodicInV_)
        {
            if (indexV <= lowestIndexV)
                moveIndexV.push_back(indexV + shiftV);
            if (indexV >= shiftV)
                moveIndexV.push_back(indexV - shiftV);
        }
    }
    
    for (unsigned int u : moveIndexU)
        for (unsigned int v : moveIndexV)
            controlPoints_[u][v] += dP;
}