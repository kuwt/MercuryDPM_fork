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

#include "NurbsWall.h"
#include "array"
#include "InteractionHandler.h"
#include "WallHandler.h"
#include "DPMBase.h"
#include "Math/ExtendedMath.h"
#include "Particles/BaseParticle.h"
using mathsFunc::square;
/*!
 * \details Make a default NurbsWall which is centered in the origin, has a length of 1, one
 * revelation, a radius of 1, turns with 1 revelation per second, is infinitely thin
 * and starts at its normal initial point.
 */
NurbsWall::NurbsWall()
{
    logger(DEBUG, "NurbsWall() constructor finished.");
}

/*!
 * \param[in] other The NurbsWall that has to be copied.
 */
NurbsWall::NurbsWall(const NurbsWall& other)
    : BaseWall(other), nurbsSurface_(other.nurbsSurface_)
{
    logger(DEBUG, "NurbsWall(const NurbsWall&) copy constructor finished.");
}

NurbsWall::NurbsWall(const NurbsSurface& nurbsSurface)
    : nurbsSurface_(nurbsSurface)
{
    logger(DEBUG, "NurbsWall(const NurbsSurface&) constructor finished.");
}

NurbsWall::~NurbsWall()
{
    logger(DEBUG, "~NurbsWall() finished, destroyed the wall.");
}

NurbsWall* NurbsWall::copy() const
{
    return new NurbsWall(*this);
}

void NurbsWall::set(const NurbsSurface& nurbsSurface) {
    nurbsSurface_ = nurbsSurface;
}


bool NurbsWall::getDistanceAndNormal(const BaseParticle& p, Mdouble& distance, Vec3D& normal_return) const
{
    //transform coordinates into position-orientation frame
    Vec3D position = p.getPosition() - getPosition();
    getOrientation().rotateBack(position);
    BaseSpecies* s = getHandler()->getDPMBase()->speciesHandler.getMixedObject(p.getSpecies(),getSpecies());
    if (nurbsSurface_.getDistance(position, p.getRadius()+s->getInteractionDistance(), distance, normal_return)) {
        getOrientation().rotate(normal_return);
        return true;
    } else {
        return false;
    }
}

/*!
 * \param[in,out] is Input stream from which the wall must be read.
 */
void NurbsWall::read(std::istream& is)
{
    BaseWall::read(is);
    std::string dummy;
    is >> dummy;
    is >> nurbsSurface_;
}

/*!
 * \param[in,out] os Output stream to which the wall must be written.
 */
void NurbsWall::write(std::ostream& os) const
{
    BaseWall::write(os);
    os << " NurbsSurface " << nurbsSurface_;
}

/*!
 * \return The string "NurbsWall".
 */
std::string NurbsWall::getName() const
{
    return "NurbsWall";
}

void NurbsWall::writeVTK (VTKContainer& vtk) const
{
    unsigned nu = 5 * nurbsSurface_.getControlPoints().size();
    unsigned nv = 5 * nurbsSurface_.getControlPoints()[0].size();
    
    const double lbU = nurbsSurface_.getLowerBoundU();
    const double ubU = nurbsSurface_.getUpperBoundU();
    const double lbV = nurbsSurface_.getLowerBoundV();
    const double ubV = nurbsSurface_.getUpperBoundV();
    
    //create points
    size_t nPoints = vtk.points.size();
    for (double u = 0; u <= nu; u++) {
        for (double v = 0; v <= nv; v++) {
            Vec3D p = nurbsSurface_.evaluate(lbU+(ubU-lbU)*u/nu, lbV+(ubV-lbV)*v/nv);
            getOrientation().rotate(p);
            p += getPosition();
            vtk.points.push_back(p);
        }
    }

    //create connectivity matrix
    //vtk.triangleStrips.reserve(vtk.triangleStrips.size()+nu);
    for (unsigned i = 0; i < nu; ++i) {
        std::vector<double> cell;
        cell.reserve(2*nv);
        for (unsigned j = 0; j <= nv; ++j) {
            cell.push_back(nPoints + j + i * (nv+1));
            cell.push_back(nPoints + j + (i+1) * (nv+1));
        }
        vtk.triangleStrips.push_back(cell);
    }
}

void NurbsWall::writeWallDetailsVTK(VTKData& data) const
{
    const std::vector<std::vector<Vec3D>>& controlPoints = nurbsSurface_.getControlPoints();
    const std::vector<std::vector<double>>& weights = nurbsSurface_.getWeights();
    
    // Reserve memory for number of points and cells about to be added
    const unsigned int np = controlPoints.size() * controlPoints[0].size();
    const unsigned int nc = (controlPoints.size() - 1) * (controlPoints[0].size() - 1);
    data.reservePoints(np, { "Weight", "ID" });
    data.reserveCells(nc);
    
    // Number of points in v-direction
    size_t nv = controlPoints[0].size();
    // Point index offset, the point indices have to be offset by the number of points already present at the start (from previously added data)
    size_t pio = data.getPoints().size();
    // Get last id only when possible and increment it, otherwise start with 0
    double id = data.getPointData()["ID"].empty() ? 0 : data.getPointData()["ID"].back() + 1;

    for (int i = 0; i < controlPoints.size(); i++)
    {
        for (int j = 0; j < controlPoints[i].size(); j++)
        {
            Vec3D p = controlPoints[i][j];
            getOrientation().rotate(p);
            p += getPosition();
            data.addToPoints(p);
            data.addToPointData("Weight", weights[i][j]);
            data.addToPointData("ID", id);
            
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