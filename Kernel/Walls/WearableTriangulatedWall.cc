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

#include "WearableTriangulatedWall.h"
#include "DPMBase.h"

WearableTriangulatedWall::WearableTriangulatedWall()
{
    logger(DEBUG, "WearableTriangulatedWall() constructed.");
}

/*!
 * \param[in] other The WearableTriangulatedWall that must be copied.
 */
WearableTriangulatedWall::WearableTriangulatedWall(const WearableTriangulatedWall& other)
        : TriangulatedWall(other)
{
    debris_ = other.debris_;
}

WearableTriangulatedWall::WearableTriangulatedWall(Mdouble lengthU, Mdouble lengthV, Mdouble resolutionU, Mdouble resolutionV)
    : TriangulatedWall()
{
    // Make resolutions fit
    auto numU = static_cast<unsigned>(lengthU / resolutionU);
    if (lengthU / numU > resolutionU)
    {
        numU++;
        resolutionU = lengthU / numU;
    }
    auto numV = static_cast<unsigned>(lengthV / resolutionV);
    if (lengthV / numV > resolutionV)
    {
        numV++;
        resolutionV = lengthV / numV;
    }
    
    // Actual number of points is one more
    numU++;
    numV++;
    
    std::vector<Vec3D> points;
    std::vector<std::vector<unsigned>> cells;
    
    points.reserve(numU * numV);
    cells.reserve(2 * (numU - 1) * (numV - 1));
    
    for (unsigned j = 0; j < numV; j++)
    {
        for (unsigned i = 0; i < numU; i++)
        {
            Vec3D p = Vec3D(i * resolutionU, j * resolutionV, 0.0);
            points.push_back(p);
            
            if (i < numU - 1 && j < numV - 1)
            {
                // Bottom left
                unsigned baseIdx = i + j * numU;
                // Bottom left, top left, top right (clockwise)
                cells.push_back({baseIdx, baseIdx + numU, baseIdx + numU + 1});
                // Bottom left, top right, bottom right (clockwise)
                cells.push_back({baseIdx, baseIdx + numU + 1, baseIdx + 1});
            }
        }
    }
    
    TriangulatedWall::set(points, cells);
    
    // Initialize debris vector
    std::vector<Vec3D> temp0(3, Vec3D(0.0, 0.0, 0.0));
    std::vector<std::vector<Vec3D>> temp1(face_.size(), temp0);
    debris_ = temp1;
}

WearableTriangulatedWall::~WearableTriangulatedWall()
{
    logger(DEBUG, "~WearableTriangulatedWall() has been called.");
}

/*!
 * \param[in] other The WearableTriangulatedWall that must be copied.
 */
WearableTriangulatedWall& WearableTriangulatedWall::operator=(const WearableTriangulatedWall& other)
{
    logger(DEBUG, "WearableTriangulatedWall::operator= called.");
    if (this == &other)
    {
        return *this;
    }
    return *(other.copy());
}

/*!
 * \return pointer to a WearableTriangulatedWall object allocated using new.
 */
WearableTriangulatedWall* WearableTriangulatedWall::copy() const
{
    return new WearableTriangulatedWall(*this);
}

/*!
 * \param[in] is The input stream from which the WearableTriangulatedWall is read, usually a restart file.
 */
void WearableTriangulatedWall::read(std::istream& is)
{
    TriangulatedWall::read(is);
}

/*!
 * \param[in] os The output stream where the WearableTriangulatedWall must be written
 *  to, usually a restart file.
 */
void WearableTriangulatedWall::write(std::ostream& os) const
{
    TriangulatedWall::write(os);
}

/*!
 * \return The string "WearableTriangulatedWall".
 */
std::string WearableTriangulatedWall::getName() const
{
    return "WearableTriangulatedWall";
}

void WearableTriangulatedWall::computeWear()
{
    // Archard wear model, from https://en.wikipedia.org/wiki/Archard_equation
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
    
        storeDebris(interaction->getContactPoint(), interaction->getNormal() * -Q);
//        storeDebris(interaction->getContactPoint(), Vec3D(0.0, 0.0, 1.0) * -Q);
    }
    
    if (!getInteractions().empty())
    {
        processDebris();
        
        // Reset debris to zero
        std::vector<Vec3D> temp0(3, Vec3D(0.0, 0.0, 0.0));
        std::fill(debris_.begin(), debris_.end(), temp0);
    }
}

void WearableTriangulatedWall::storeDebris(Vec3D P, const Vec3D& debris)
{
    // Not sure if needed for triangulated wall?
//    // Rotate global point back to local coordinate system.
//    P -= getPosition();
//    getOrientation().rotateBack(P);
    
    Face* face;
    unsigned it = 0;
    for (auto& f : face_)
    {
        const Vec3D a = *f.vertex[0] - P;
        const Vec3D b = *f.vertex[1] - P;
        const Vec3D c = *f.vertex[2] - P;
        
        const Vec3D u = Vec3D::cross(b, c);
        const Vec3D v = Vec3D::cross(c, a);
        const Vec3D w = Vec3D::cross(a, b);
        
        // When all are pointing in same direction, the point is inside (or on edge of) the triangle.
        if (!(Vec3D::dot(u, v) < 0.0 || Vec3D::dot(u, w) < 0.0))
        {
            face = &f;
            break;
        }
        
        it++;
    }
    
    // temporary because this should not happen
    if (it >= face_.size())
    {
//        logger(INFO, "does this happen?");
        return;
    }
    
    // For now every vertex simply gets 1/3 of debris
    debris_[it][0] += debris / 3.0;
    debris_[it][1] += debris / 3.0;
    debris_[it][2] += debris / 3.0;
}

void WearableTriangulatedWall::processDebris()
{
    for (unsigned i = 0; i < face_.size(); i++)
    {
        *face_[i].vertex[0] += debris_[i][0];
        *face_[i].vertex[1] += debris_[i][1];
        *face_[i].vertex[2] += debris_[i][2];
    }
    
    //set normals and positions
    for (auto& face: face_)
    {
        face.normal = Vec3D::getUnitVector(
                Vec3D::cross(*face.vertex[1] - *face.vertex[0], *face.vertex[2] - *face.vertex[0]));
    }
    
    //set edge normals (inwards facing)
    for (auto& face: face_)
    {
        face.edgeNormal[0] = Vec3D::getUnitVector(Vec3D::cross(face.normal, *face.vertex[1] - *face.vertex[0]));
        face.edgeNormal[1] = Vec3D::getUnitVector(Vec3D::cross(face.normal, *face.vertex[2] - *face.vertex[1]));
        face.edgeNormal[2] = Vec3D::getUnitVector(Vec3D::cross(face.normal, *face.vertex[0] - *face.vertex[2]));
    }
}