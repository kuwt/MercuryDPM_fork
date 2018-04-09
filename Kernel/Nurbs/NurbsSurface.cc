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

#include "NurbsSurface.h"
#include "NurbsUtils.h"
#include "Logger.h"
using namespace NurbsUtils;

NurbsSurface::NurbsSurface(const std::vector<double>& knotsU, const std::vector<double>& knotsV,
             const std::vector<std::vector<Vec3D>>& controlPoints, const std::vector<std::vector<double>>& weights)
    : controlPoints_(controlPoints), weights_(weights), knotsU_(knotsU), knotsV_(knotsV)
{
    if ( controlPoints.size()<2 || controlPoints[0].size()<2 ) {
        logger(ERROR,"At least two control points are necessary");
    }

    if ( controlPoints.size()!=weights.size() || controlPoints[0].size()!=weights[0].size() ) {
        logger(ERROR,"Numer of control points and weights must match");
    }

    if ( knotsU.size()<2 || knotsV.size()<2 ) {
        logger(ERROR,"At least two knots are necessary");
    }

    if (!isKnotVectorMonotonic(knotsU) || !isKnotVectorMonotonic(knotsV)) {
        logger(ERROR,"Knot vector(s) is not monotonic");
    }

    if (knotsU.size() - controlPoints.size() -1 < 1 || knotsV.size() - controlPoints[0].size() -1 < 1) {
        logger(ERROR,"Degree has to be at least 1");
    }

    degreeU_ = knotsU.size() - controlPoints.size() - 1;
    degreeV_ = knotsV.size() - controlPoints[0].size() - 1;

    logger(INFO,"Created Nurbs surface.");
    logger(INFO,"  %x% knots",knotsU.size(), knotsV.size());
    logger(INFO,"  %x% control points",controlPoints.size(), controlPoints[0].size());
    logger(INFO,"  %x% degrees",degreeU_,degreeV_);
}

Vec3D NurbsSurface::evaluate(double u, double v) {

    Vec3D point = {0,0,0};
    double temp = 0;

    // Find span and non-zero basis functions
    int spanU = findSpan(degreeU_, knotsU_, u);
    int spanV = findSpan(degreeV_, knotsV_, v);
    std::vector<double> Nu, Nv;
    bsplineBasis(degreeU_, spanU, knotsU_, u, Nu);
    bsplineBasis(degreeV_, spanV, knotsV_, v, Nv);
    // make linear combination
    for (int l = 0; l <= degreeV_; ++l) {
        for (int k = 0; k <= degreeU_; ++k) {
            double weight = Nv[l] * Nu[k] * weights_[spanU - degreeU_ + k][spanV - degreeV_ + l];
            point += weight * controlPoints_[spanU - degreeU_ + k][spanV - degreeV_ + l];
            temp += weight;
        }
    }
    point /= temp;
    return point;
}
