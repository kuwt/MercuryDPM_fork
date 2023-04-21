//Copyright (c) 2013-2021, The MercuryDPM Developers Team. All rights reserved.
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

#include "WearableNurbsWall.h"
#include "WallHandler.h"
#include "DPMBase.h"
#include "Nurbs/NurbsUtils.h"

WearableNurbsWall::WearableNurbsWall()
{
    logger(DEBUG, "WearableNurbsWall() constructor finished.");
}

/*!
 * \param[in] other The WearableNurbsWall that has to be copied.
 */
WearableNurbsWall::WearableNurbsWall(const WearableNurbsWall& other)
        : NurbsWall(other)
{
    localDebris_ = other.localDebris_;
    logger(DEBUG, "WearableNurbsWall(const WearableNurbsWall&) copy constructor finished.");
}

/*!
 * @param nurbsSurface The NURBS surface which defines the shape of the wall
 */
WearableNurbsWall::WearableNurbsWall(const NurbsSurface& nurbsSurface)
        : NurbsWall(nurbsSurface)
{
    logger(DEBUG, "WearableNurbsWall(const NurbsSurface&) constructor finished.");
}

WearableNurbsWall::~WearableNurbsWall()
{
    logger(DEBUG, "~WearableNurbsWall() finished, destroyed the wall.");
}

WearableNurbsWall::WearableNurbsWall(Mdouble lengthU, Mdouble lengthV, Mdouble resolutionU, Mdouble resolutionV, bool periodicU, bool periodicV)
{
    set(lengthU, lengthV, resolutionU, resolutionV, periodicU, periodicV);
}

void WearableNurbsWall::set(Mdouble lengthU, Mdouble lengthV, Mdouble resolutionU, Mdouble resolutionV, bool periodicU, bool periodicV)
{
    std::vector<std::vector<Vec3D>> controlPoints;
    std::vector<std::vector<Mdouble>> weights;
    
    int numU = static_cast<int>(lengthU / resolutionU);
    if (lengthU / numU > resolutionU)
    {
        numU++;
        resolutionU = lengthU / numU;
    }
    int numV = static_cast<int>(lengthV / resolutionV);
    if (lengthV / numV > resolutionV)
    {
        numV++;
        resolutionV = lengthV / numV;
    }
    
    for (int i = 0; i <= numU; i++)
    {
        std::vector<Vec3D> tempControlPoints;
        std::vector<Mdouble> tempWeights;
        for (int j = 0; j <= numV; j++)
        {
            tempControlPoints.push_back(Vec3D(i * resolutionU, j * resolutionV, 0.0));
            tempWeights.push_back(1.0);
        }
        controlPoints.push_back(tempControlPoints);
        weights.push_back(tempWeights);
    }
    
    // Make periodic when asked for, only clamp when not periodic.
    // Making periodic will unclamp it when needed, but better to start of unclamped already otherwise
    // warning messages are given to the user.
    NurbsSurface ns(controlPoints, weights, 3, 3, !periodicU, !periodicU, !periodicV, !periodicV);
    if (periodicU)
        ns.makePeriodicContinuousInU();
    if (periodicV)
        ns.makePeriodicContinuousInV();
    nurbsSurface_ = ns;
    
    // Initialize local debris with zeros
    const int sizeU = nurbsSurface_.getControlPoints().size();
    const int sizeV = nurbsSurface_.getControlPoints()[0].size();
    for (int i = 0; i < sizeU; i++)
    {
        std::vector<Mdouble> temp(sizeV, 0.0);
        localDebris_.push_back(temp);
    }
}

/*!
 * \param[in,out] is Input stream from which the wall must be read.
 */
void WearableNurbsWall::read(std::istream& is)
{
    NurbsWall::read(is);
    std::string dummy;
    is >> dummy;
    
    unsigned int nu, nv;
    is >> dummy >> nu;
    is >> dummy >> nv;
    
    std::vector<std::vector<Mdouble>> localDebris;
    localDebris.resize(nu);
    for (auto& d0 : localDebris)
    {
        d0.resize(nv);
        for (auto& d : d0) is >> d;
    }
    
    localDebris_ = localDebris;
}

/*!
 * \param[in,out] os Output stream to which the wall must be written.
 */
void WearableNurbsWall::write(std::ostream& os) const
{
    NurbsWall::write(os);
    os << " Debris ";
    os << "nu " << localDebris_.size() << ' ';
    os << "nv " << localDebris_[0].size() << ' ';
    for (const auto d0 : localDebris_) for (const auto d : d0) os << d << ' ';
}

std::string WearableNurbsWall::getName() const
{
    return "WearableNurbsWall";
}

WearableNurbsWall * WearableNurbsWall::copy() const
{
    return new WearableNurbsWall(*this);
}

void WearableNurbsWall::computeWear()
{
    // From https://en.wikipedia.org/wiki/Archard_equation
    // Q = KWL/H
    // Q is the total volume of wear debris produced
    // K is a dimensionless constant (typically mild wear 1e-8, severe wear 1e-2)
    // W is the total normal load
    // L is the sliding distance
    // H is the hardness of the softest contacting surfaces
    
    const Mdouble dt = getHandler()->getDPMBase()->getTimeStep();
    
    for (BaseInteraction* interaction : getInteractions())
    {
        // Ignore interactions from periodic particles. (Although it seems Q evaluates to 0 anyway (absolute normal force equals 0), but just to be sure)
        if (static_cast<BaseParticle*>(interaction->getP())->getPeriodicFromParticle() != nullptr)
            continue;
            
        Mdouble K = 1.0e-6; // For testing only, should be set by user
        Mdouble W = interaction->getAbsoluteNormalForce();
        // Pythagoras to get tangential magnitude from vector and normal magnitude
        Mdouble tangentialRelativeVelocity = std::sqrt(interaction->getRelativeVelocity().getLengthSquared() -
                interaction->getNormalRelativeVelocity() * interaction->getNormalRelativeVelocity());
        Mdouble L = tangentialRelativeVelocity * dt;
        Mdouble H = 1.0; // For testing only, should be set by user
        Mdouble Q = K * W * L / H;
        
        storeDebris(interaction->getContactPoint(), Q);
    }
    
//    if (!getInteractions().empty() && getHandler()->getDPMBase()->getNumberOfTimeSteps() % 100 == 0)
//    {
//        processDebris();
//        // Reset debris
//        for (auto& d0 : localDebris_)
//            for (auto& d : d0)
//                d = 0.0;
//    }
}

void WearableNurbsWall::storeDebris(Vec3D P, const Mdouble debris)
{
    // Rotate global point back to local coordinate system.
    P -= getPosition();
    getOrientation().rotateBack(P);
    
    // Snap the position to the closest grid point. Assuming the control points only ever move in z-direction, so the
    // x- and y-position don't change, and it was created with a given resolution in u and v.
    const std::vector<std::vector<Vec3D>>& controlPoints = nurbsSurface_.getControlPoints();
    const Mdouble resU = std::abs(controlPoints[1][0].X - controlPoints[0][0].X);
    const Mdouble resV = std::abs(controlPoints[0][1].Y - controlPoints[0][0].Y);
    
    // --- Add to closest grid point -----------------------------------------------------------------------------------
    auto indexU = static_cast<const unsigned>(std::round((P.X - controlPoints[0][0].X) / resU));
    auto indexV = static_cast<const unsigned>(std::round((P.Y - controlPoints[0][0].Y) / resV));
    localDebris_[indexU][indexV] += debris;
    
    // --- Add to 4 closest grid points, weighted by distance ----------------------------------------------------------
//    const Mdouble pU = (P.X - controlPoints[0][0].X) / resU;
//    const Mdouble pV = (P.Y - controlPoints[0][0].Y) / resV;
//
//    // Left and right distance in u and v
//    const Mdouble pUL = pU - std::floor(pU);
//    const Mdouble pVL = pV - std::floor(pV);
//    const Mdouble pUR = 1.0 - pUL;
//    const Mdouble pVR = 1.0 - pVL;
//
//    // Power parameter
//    const Mdouble p = 1.0;
//
//    // Normalizing "length" in u and v
//    const Mdouble lU = std::pow(pUL, p) + std::pow(pUR, p);
//    const Mdouble lV = std::pow(pVL, p) + std::pow(pVR, p);
//
//    // Weights for left and right in u and v
//    const Mdouble wUL = 1.0 - std::pow(pUL, p) / lU;
//    const Mdouble wUR = 1.0 - std::pow(pUR, p) / lU;
//    const Mdouble wVL = 1.0 - std::pow(pVL, p) / lV;
//    const Mdouble wVR = 1.0 - std::pow(pVR, p) / lV;
//
//    const unsigned firstIndexU = std::floor(pU);
//    const unsigned firstIndexV = std::floor(pV);
//    localDebris_[firstIndexU    ][firstIndexV    ] += debris * wUL * wVL;
//    localDebris_[firstIndexU    ][firstIndexV + 1] += debris * wUL * wVR;
//    localDebris_[firstIndexU + 1][firstIndexV    ] += debris * wUR * wVL;
//    localDebris_[firstIndexU + 1][firstIndexV + 1] += debris * wUR * wVR;
}

void WearableNurbsWall::processDebris()
{
    std::vector<std::vector<Vec3D>> controlPoints = nurbsSurface_.getControlPoints();
    Mdouble totalDebris = 0.0;
    
    const int numU = controlPoints.size();
    const int numV = controlPoints[0].size();
    
    for (int i = 0; i < numU; i++)
    {
        for (int j = 0; j < numV; j++)
        {
            controlPoints[i][j].Z = -localDebris_[i][j];
            totalDebris += localDebris_[i][j];
        }
    }

    const Mdouble volumeChange = getVolumeUnderSurface(nurbsSurface_.getKnotsU(), nurbsSurface_.getKnotsV(), controlPoints, nurbsSurface_.getWeights());
    const Mdouble ratio = totalDebris / volumeChange;

    for (int i = 0; i < numU; i++)
    {
        for (int j = 0; j < numV; j++)
        {
            nurbsSurface_.moveControlPoint(i, j, Vec3D(0.0, 0.0, -localDebris_[i][j] * ratio), true);
        }
    }
}

Mdouble WearableNurbsWall::getVolumeUnderSurface(const std::vector<Mdouble>& knotsU, const std::vector<Mdouble>& knotsV,
                                                 const std::vector<std::vector<Vec3D>>& controlPoints, const std::vector<std::vector<Mdouble>>& weights) const
{
    // times 10 is arbitrary number for more accuracy.
    const int numU = knotsU.size() * 10; // Was 10 lately a lot for many test simulations
    const int numV = knotsV.size() * 10;
    
    const int degreeU = knotsU.size() - controlPoints.size() - 1;
    const int degreeV = knotsV.size() - controlPoints[0].size() - 1;
    
    const Mdouble lbU = knotsU[degreeU];
    const Mdouble ubU = knotsU[knotsU.size() - degreeU - 1];
    const Mdouble lbV = knotsV[degreeV];
    const Mdouble ubV = knotsV[knotsV.size() - degreeV - 1];
    
    Mdouble volume = 0.0;
    for (int i = 0; i < numU; i++)
    {
        for (int j = 0; j < numV; j++)
        {
            const Mdouble u0 = lbU + (ubU - lbU) * static_cast<Mdouble>(i)   / static_cast<Mdouble>(numU);
            const Mdouble u1 = lbU + (ubU - lbU) * static_cast<Mdouble>(i+1) / static_cast<Mdouble>(numU);
            const Mdouble v0 = lbV + (ubV - lbV) * static_cast<Mdouble>(j)   / static_cast<Mdouble>(numV);
            const Mdouble v1 = lbV + (ubV - lbV) * static_cast<Mdouble>(j+1) / static_cast<Mdouble>(numV);
            const Vec3D p1 = NurbsUtils::evaluate(u0, v0, knotsU, knotsV, controlPoints, weights); // Bottom left
            const Vec3D p2 = NurbsUtils::evaluate(u0, v1, knotsU, knotsV, controlPoints, weights); // Top left
            const Vec3D p3 = NurbsUtils::evaluate(u1, v0, knotsU, knotsV, controlPoints, weights); // Bottom right
            const Vec3D p4 = NurbsUtils::evaluate(u1, v1, knotsU, knotsV, controlPoints, weights); // Top right

            const Mdouble projectedArea1 = 0.5 * (p1.X * p3.Y - p3.X * p1.Y + p2.X * p1.Y - p1.X * p2.Y + p3.X * p2.Y - p2.X * p3.Y);
            const Mdouble projectedArea2 = 0.5 * (p4.X * p3.Y - p3.X * p4.Y + p2.X * p4.Y - p4.X * p2.Y + p3.X * p2.Y - p2.X * p3.Y);
            const Mdouble averageHeight1 = (p1.Z + p2.Z + p3.Z) / 3.0;
            const Mdouble averageHeight2 = (p4.Z + p2.Z + p3.Z) / 3.0;
            // Area can be negative (depending on order of points considered), and z is always negative assuming only
            // downward movement of control points, therefore take absolute value.
            const Mdouble volume1 = std::abs(averageHeight1 * projectedArea1);
            const Mdouble volume2 = std::abs(averageHeight2 * projectedArea2);
            volume += volume1 + volume2;
        }
    }
    return volume;
}

Mdouble WearableNurbsWall::getVolumeUnderSurfaceX(const std::vector<Mdouble>& knotsU, const std::vector<Mdouble>& knotsV,
                                                  const std::vector<std::vector<Vec3D>>& controlPoints, const std::vector<std::vector<Mdouble>>& weights) const
{
    // times 10 is arbitrary number for more accuracy.
    const int numU = knotsU.size() * 10; // Was 10 lately a lot for many test simulations
    const int numV = knotsV.size() * 10;
    
    const int degreeU = knotsU.size() - controlPoints.size() - 1;
    const int degreeV = knotsV.size() - controlPoints[0].size() - 1;
    
    const Mdouble lbU = knotsU[degreeU];
    const Mdouble ubU = knotsU[knotsU.size() - degreeU - 1];
    const Mdouble lbV = knotsV[degreeV];
    const Mdouble ubV = knotsV[knotsV.size() - degreeV - 1];
    
    Mdouble volume = 0.0;
    for (int i = 0; i < numU; i++)
    {
        for (int j = 0; j < numV; j++)
        {
            const Mdouble u0 = lbU + (ubU - lbU) * static_cast<Mdouble>(i)   / static_cast<Mdouble>(numU);
            const Mdouble u1 = lbU + (ubU - lbU) * static_cast<Mdouble>(i+1) / static_cast<Mdouble>(numU);
            const Mdouble v0 = lbV + (ubV - lbV) * static_cast<Mdouble>(j)   / static_cast<Mdouble>(numV);
            const Mdouble v1 = lbV + (ubV - lbV) * static_cast<Mdouble>(j+1) / static_cast<Mdouble>(numV);
            const Vec3D p1 = NurbsUtils::evaluate(u0, v0, knotsU, knotsV, controlPoints, weights); // Bottom left
            const Vec3D p2 = NurbsUtils::evaluate(u0, v1, knotsU, knotsV, controlPoints, weights); // Top left
            const Vec3D p3 = NurbsUtils::evaluate(u1, v0, knotsU, knotsV, controlPoints, weights); // Bottom right
            const Vec3D p4 = NurbsUtils::evaluate(u1, v1, knotsU, knotsV, controlPoints, weights); // Top right
            
            // Only place with differences, x and z interchanged
            const Mdouble projectedArea1 = 0.5 * (p1.Z * p3.Y - p3.Z * p1.Y + p2.Z * p1.Y - p1.Z * p2.Y + p3.Z * p2.Y - p2.Z * p3.Y);
            const Mdouble projectedArea2 = 0.5 * (p4.Z * p3.Y - p3.Z * p4.Y + p2.Z * p4.Y - p4.Z * p2.Y + p3.Z * p2.Y - p2.Z * p3.Y);
            const Mdouble averageHeight1 = (p1.X + p2.X + p3.X) / 3.0;
            const Mdouble averageHeight2 = (p4.X + p2.X + p3.X) / 3.0;
            // Area can be negative (depending on order of points considered), and z is always negative assuming only
            // downward movement of control points, therefore take absolute value.
            const Mdouble volume1 = std::abs(averageHeight1 * projectedArea1);
            const Mdouble volume2 = std::abs(averageHeight2 * projectedArea2);
            volume += volume1 + volume2;
        }
    }
    return volume;
}

void WearableNurbsWall::writeWallDetailsVTK(VTKData& data) const
{
    // Visualised as a flat surface, which can be coloured according to the debris at each point.
    // The points are simply the control points with the z-position set to 0.
    
    const std::vector<std::vector<Vec3D>>& points = nurbsSurface_.getControlPoints();
    
    // Reserve memory for number of points and cells about to be added
    const unsigned np = points.size() * points[0].size();
    const unsigned nc = (points.size() - 1) * (points[0].size() - 1);
    data.reservePoints(np, { "Debris" });
    data.reserveCells(nc);
    
    // Number of points in v-direction
    size_t nv = points[0].size();
    // Point index offset, the point indices have to be offset by the number of points already present at the start (from previously added data)
    size_t pio = data.getPoints().size();
    
    for (int i = 0; i < points.size(); i++)
    {
        for (int j = 0; j < points[i].size(); j++)
        {
            Vec3D p = points[i][j];
            p.Z = 0;
            getOrientation().rotate(p);
            p += getPosition();
            data.addToPoints(p);
            data.addToPointData("Debris", localDebris_[i][j]);
            
            if (i > 0 && j > 0)
            {
                // Basic 2D/1D mapping for indexing: nv * i + j (+ point index offset)
                // 4 points to form a rectangle: point, point to the left, point down, point to the left and down
                data.addToConnectivity({ nv*i+j+pio, nv*(i-1)+j+pio, nv*i+j-1+pio, nv*(i-1)+j-1+pio });
                data.addToTypes(8);
            }
        }
    }
}