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

#include <MercuryBase.h>
#include "TriangleWall.h"
#include "InteractionHandler.h"
#include "WallHandler.h"
#include "DPMBase.h"
//#include "Species/BaseSpecies.h"
#include "Particles/BaseParticle.h"

bool setDistance(Mdouble& distance, const Vec3D& branch, Mdouble distanceMax)
{
    const Mdouble distance2 = branch.getLengthSquared();
    if (distance2 > distanceMax * distanceMax) return false;
    distance = sqrt(distance2);
    return true;
}


/*!
 * \param[in] p BaseParticle we want to calculate the distance and whether it collided of.
 * \param[out] distance The distance of the BaseParticle to this wall.
 * \param[out] normal_return If there was a collision, the normal vector to this wall will be placed here.
 * \return A boolean which says whether or not there was a collision.
 * \details This function computes whether or not there is a collision between
 * a given BaseParticle and this TriangleWall. If there is a collision, this
 * function also computes the distance between the BaseParticle and TriangleWall
 * and the normal of the TriangleWall at the intersection point. It does
 * this by calling TriangleWall::getDistanceAndNormal(const Vec3D& , Mdouble , Mdouble&, Vec3D&) const.
 * Since this function should be called before calculating any
 * Particle-Wall interactions, it can also be used to set the normal vector in
 * case of curved walls.
 */
bool TriangleWall::getDistanceAndNormal(const BaseParticle& p, Mdouble& distance, Vec3D& normal) const
{
    const Vec3D position = p.getPosition();// - getPosition();
    //getOrientation().rotateBack(position);
    const Mdouble signedDistance = Vec3D::dot(position-vertex_[0], faceNormal_);
    distance = fabs(signedDistance);
    if (distance >= p.getWallInteractionRadius(this))
    {
        return false;
    }
    
    const Mdouble distanceMax2 = mathsFunc::square(p.getWallInteractionRadius(this));
    const Vec3D edgeBranch0 = position - vertex_[0];
    const Vec3D edgeBranch1 = position - vertex_[1];
    const Vec3D edgeBranch2 = position - vertex_[2];
    const Mdouble edgeDistance0 = Vec3D::dot(edgeBranch0, edgeNormal_[0]);
    const Mdouble edgeDistance1 = Vec3D::dot(edgeBranch1, edgeNormal_[1]);
    const Mdouble edgeDistance2 = Vec3D::dot(edgeBranch2, edgeNormal_[2]);
    //logger(INFO,"n % d % P % E % e % % %",-faceNormal_,signedDistance, position,edgeNormal_[0],edgeDistance0,edgeDistance1,edgeDistance2);
    if (edgeDistance0 > 0)
    {
        if (edgeDistance1 > 0)
        {
            //logger(INFO,"contact with vertex 1, dist % norm %");
            const Mdouble distance2 = edgeBranch1.getLengthSquared();
            if (distance2 > distanceMax2) return false;
            distance = sqrt(distance2);
            normal = edgeBranch1 / -distance;
        }
        else if (edgeDistance2 > 0)
        {
            //logger(INFO,"contact with vertex 0");
            const Mdouble distance2 = edgeBranch0.getLengthSquared();
            if (distance2 > distanceMax2) return false;
            distance = sqrt(distance2);
            normal = edgeBranch0 / -distance;
        }
        else
        {
            //logger(INFO,"contact with edge 0");
            normal = edge_[0] * Vec3D::dot(edgeBranch0, edge_[0]) - edgeBranch0;
            const Mdouble distance2 = normal.getLengthSquared();
            if (distance2 > distanceMax2) return false;
            distance = sqrt(distance2);
            normal /= distance;
        }
    }
    else if (edgeDistance1 > 0)
    {
        if (edgeDistance2 > 0)
        {
            //logger(INFO,"contact with vertex 2");
            const Mdouble distance2 = edgeBranch2.getLengthSquared();
            if (distance2 > distanceMax2) return false;
            distance = sqrt(distance2);
            normal = edgeBranch2 / -distance;
        }
        else
        {
            //logger(INFO,"contact with edge 1");
            normal = edge_[1] * Vec3D::dot(edgeBranch1, edge_[1]) - edgeBranch1;
            const Mdouble distance2 = normal.getLengthSquared();
            if (distance2 > distanceMax2) return false;
            distance = sqrt(distance2);
            normal /= distance;
        }
    }
    else if (edgeDistance2 > 0)
    {
        //logger(INFO,"contact with edge 2");
        normal = edge_[2] * Vec3D::dot(edgeBranch2, edge_[2]) - edgeBranch2;
        const Mdouble distance2 = normal.getLengthSquared();
        if (distance2 > distanceMax2) return false;
        distance = sqrt(distance2);
        normal /= distance;
    }
    else
    {
        //logger(INFO,"contact with face");
        normal = (signedDistance >= 0) ? -faceNormal_ : faceNormal_;
    }
    return true;
}

void TriangleWall::rotate(const Vec3D& angularVelocityDt)
{
    if (!angularVelocityDt.isZero())
    {
        BaseInteractable::rotate(angularVelocityDt);
        updateVertexAndNormal();
    }
}

/*!
 * \param[in] is The input stream from which the TriangleWall is read, usually a restart file.
 */
void TriangleWall::read(std::istream& is)
{
    BaseWall::read(is);
    std::string dummy;
    is >> dummy;
    for (int i = 0; i < 3; i++)
    {
        is >> vertexInLabFrame_[i];
    }
    updateVertexAndNormal();
}

/*!
 * \param[in] os The output stream where the TriangleWall must be written
 *  to, usually a restart file.
 */
void TriangleWall::write(std::ostream& os) const
{
    BaseWall::write(os);
    os << " vertexInLabFrame ";
    for (int i = 0; i < 3; i++)
    {
        os << ' ' << vertexInLabFrame_[i];
    }
}

void TriangleWall::writeVTK(VTKContainer& vtk) const
{
    const unsigned long s = vtk.points.size();
    vtk.points.reserve(s + 3);
    for (auto v : vertex_)
    {
        vtk.points.push_back(v);
    }
    std::vector<double> cell;
    cell.reserve(3);
    cell.push_back(s);
    cell.push_back(s + 1);
    cell.push_back(s + 2);
    //cell.push_back(s);
    vtk.triangleStrips.push_back(cell);
}

void TriangleWall::setVertices(const Vec3D A, const Vec3D B, const Vec3D C)
{
    setPosition((A + B + C) / 3);
    setOrientation({1, 0, 0, 0});
    vertexInLabFrame_[0] = A - getPosition();
    vertexInLabFrame_[1] = B - getPosition();
    vertexInLabFrame_[2] = C - getPosition();
    //vertexInLabFrame_[0] = A;
    //vertexInLabFrame_[1] = B;
    //vertexInLabFrame_[2] = C;
    updateVertexAndNormal();
}

/*!
 * \details Moves (displaces) the interacable a given distance.
 *          Note, this just updates the position by the move.
 * \param[in] move  Reference to Vec3D which is the distance to move the
 *            interactable.
 */
void TriangleWall::move(const Vec3D& move)
{
    BaseInteractable::move(move);
    updateVertexAndNormal();
}

void TriangleWall::setVertices(const Vec3D A, const Vec3D B, const Vec3D C, const Vec3D position)
{
    setPosition(position);
    setOrientation({1, 0, 0, 0});
    vertexInLabFrame_[0] = A - getPosition();
    vertexInLabFrame_[1] = B - getPosition();
    vertexInLabFrame_[2] = C - getPosition();
    updateVertexAndNormal();
}

/**
 * This function should be called after setting either position_ or vertexInLabFrame_.
 *  - vertex is set to the vertex position in the real coordinate system (rotated and shifted)
 *  - vertexMin_/vertexMax_ is set to define a bounding box around the wall (for contact detection)
 *  - edge_, edgeNormal_ and faceNormal_ is updated (stored for quick computation of contact point)
 */
void TriangleWall::updateVertexAndNormal()
{
    for (int i = 0; i < 3; i++)
    {
        vertex_[i] = vertexInLabFrame_[i];
        getOrientation().rotate(vertex_[i]);
        vertex_[i] += getPosition();
    }
    vertexMin_ = Vec3D::min(Vec3D::min(vertex_[0], vertex_[1]), vertex_[2]);
    vertexMax_ = Vec3D::max(Vec3D::max(vertex_[0], vertex_[1]), vertex_[2]);
    
    edge_ = {vertex_[1] - vertex_[0], vertex_[2] - vertex_[1], vertex_[0] - vertex_[2]};
    faceNormal_ = Vec3D::cross(edge_[0], edge_[1]);
    faceNormal_.normalise();
    
    for (int i = 0; i < 3; i++)
    {
        edgeNormal_[i] = Vec3D::cross(edge_[i], faceNormal_);
        edge_[i].normalise();
        edgeNormal_[i].normalise();
    }
    //logger(INFO,"vertex %,%,% edge %,%,% face %",vertex_[0],vertex_[1],vertex_[2],edgeNormal_[0],edgeNormal_[1],edgeNormal_[2],faceNormal_);
}
