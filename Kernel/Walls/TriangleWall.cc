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
    const Vec3D position = p.getPosition(); // note: no transfer to lab coordinates
    const Mdouble distanceMax = p.getWallInteractionRadius(this);

    // compute distance from face
    //get distance from particle position to the face
    const Mdouble signedDistance = Vec3D::dot(position-vertex_[0], faceNormal_);
    distance = fabs(signedDistance);

    // check if any contact is possible
    if (distance >= distanceMax) return false;

    // compute distance from edges
    const std::array<Vec3D,3> edgeBranch {position - vertex_[0], position - vertex_[1], position - vertex_[2]};
    const std::array<double,3> edgeDistance {Vec3D::dot(edgeBranch[0], edgeNormal_[0]), Vec3D::dot(edgeBranch[1], edgeNormal_[1]), Vec3D::dot(edgeBranch[2], edgeNormal_[2])};

    // find edge with largest distance (this will be the edge if there is a edge contact)
    const size_t id = (edgeDistance[1] > edgeDistance[2]) ?
            (edgeDistance[0] > edgeDistance[1] ? 0 : 1) : (edgeDistance[0] > edgeDistance[2] ? 0 : 2);

    // check if there will be contact
    if (edgeDistance[id] > distanceMax) return false;

    // determine if it is a face contact
    const Vec3D posProjected = position - signedDistance * faceNormal_;
    if (edgeDistance[id] <= 0 && isInsideTriangle(posProjected)){
        normal = (signedDistance >= 0) ? -faceNormal_ : faceNormal_;
        return true;
    }

    // determine if it is an edge or vertex contact
    const double positionAlongEdge = Vec3D::dot(edgeBranch[id], edge_[id]);
    if (positionAlongEdge < 0) {
        //possible contact with left vertex
        distance = edgeBranch[id].getLength();
        if (distance > distanceMax) return false;
        normal = edgeBranch[id] / -distance;
    } else if (positionAlongEdge > edgeLength_[id]) {
        //contact with right vertex
        const size_t idRight = (id + 1) % 3;
        distance = edgeBranch[idRight].getLength();
        if (distance > distanceMax) return false;
        normal = edgeBranch[idRight] / -distance;
    } else {
        // edge contact
        normal = edge_[id] * positionAlongEdge - edgeBranch[id];
        distance = normal.getLength();
        if (distance > distanceMax) return false;
        normal /= distance;
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
        edgeLength_[i] = edge_[i].getLength();
        edge_[i] /= edgeLength_[i];
        edgeNormal_[i].normalise();
    }
    //logger(INFO,"vertex %,%,% edge %,%,% face %",vertex_[0],vertex_[1],vertex_[2],edgeNormal_[0],edgeNormal_[1],edgeNormal_[2],faceNormal_);
}

bool TriangleWall::isLocal(Vec3D& min, Vec3D& max) const
{
    min = vertexMin_;
    max = vertexMax_;
    return true;
}

bool TriangleWall::isInsideTriangle(const Vec3D &point) const
{
    const Vec3D branch0 = point - vertex_[0];
    const Vec3D branch1 = point - vertex_[1];
    const Vec3D branch2 = point - vertex_[2];

    //compute total area
    Mdouble s = sqrt(Vec3D::cross(vertex_[1] - vertex_[0], vertex_[2] - vertex_[1]).getLengthSquared());
    //compute barycentric coordinates
    Mdouble s0 = sqrt(Vec3D::cross(branch0, branch1).getLengthSquared())/s;
    Mdouble s1 = sqrt(Vec3D::cross(branch1, branch2).getLengthSquared())/s;
    Mdouble s2 = sqrt(Vec3D::cross(branch2, branch0).getLengthSquared())/s;
    return (1 > s0 > 0 && 1 > s1 > 0 && 1 > s2 > 0);
}
