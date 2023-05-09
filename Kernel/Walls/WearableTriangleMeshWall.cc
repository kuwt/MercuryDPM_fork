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

#include "WearableTriangleMeshWall.h"
#include "DPMBase.h"
#include <math.h>

WearableTriangleMeshWall::WearableTriangleMeshWall(const std::vector<Vec3D>& points,
                                                   const std::vector<std::array<unsigned , 3>>& cells,
                                                   const ParticleSpecies* species)
                                                   : TriangleMeshWall(points, cells, species)
{
}

WearableTriangleMeshWall::WearableTriangleMeshWall(const Vec3D& P0, const Vec3D& P1, const Vec3D& P2, Mdouble resolutionU,
                                                   Mdouble resolutionV, const ParticleSpecies* species,
                                                   bool periodicInU, bool periodicInV) :
                                                   TriangleMeshWall(P0, P1, P2, resolutionU, resolutionV, species, periodicInU, periodicInV)
{
}

WearableTriangleMeshWall::WearableTriangleMeshWall(const WearableTriangleMeshWall& other) : TriangleMeshWall(other)
{
    wearCoefficient_ = other.wearCoefficient_;
    hardness_ = other.hardness_;
    wearAcceleration_ = other.wearAcceleration_;
}

WearableTriangleMeshWall& WearableTriangleMeshWall::operator=(const WearableTriangleMeshWall& other)
{
    if (this == &other)
    {
        return *this;
    }
    return *(other.copy());
}

WearableTriangleMeshWall* WearableTriangleMeshWall::copy() const
{
    return new WearableTriangleMeshWall(*this);
}

void WearableTriangleMeshWall::read(std::istream& is)
{
    TriangleMeshWall::read(is);
    
    std::string dummy;
    is >> dummy >> wearCoefficient_ >> dummy >> hardness_ >> dummy >> wearAcceleration_;
}

void WearableTriangleMeshWall::write(std::ostream& os) const
{
    TriangleMeshWall::write(os);
    
    os << " wearCoefficient " << wearCoefficient_ << " hardness " << hardness_ << " wearAcceleration " << wearAcceleration_;
}

std::string WearableTriangleMeshWall::getName() const
{
    return "WearableTriangleMeshWall";
}

void WearableTriangleMeshWall::computeWear()
{
    // Archard wear model, from https://en.wikipedia.org/wiki/Archard_equation
    // Q = KWL/H
    // Q is the total volume of wear debris produced
    // K is a dimensionless constant (typically mild wear 1e-8, severe wear 1e-2)
    // W is the total normal load
    // L is the sliding distance
    // H is the hardness of the softest contacting surfaces
    
    // Initialize debris vector with all zero Vec3Ds.
    std::vector<Vec3D> debrisContainer(vertices_.size(), Vec3D(0.0, 0.0, 0.0));
    Mdouble totalDebris = 0.0;
    
    for (auto& t : triangles_)
    {
        for (BaseInteraction* interaction : t.wall.getInteractions())
        {
            // Ignore old interactions and interactions from periodic particles.
            if (interaction->getTimeStamp() <= getHandler()->getDPMBase()->getNumberOfTimeSteps() ||
                static_cast<BaseParticle*>(interaction->getP())->getPeriodicFromParticle() != nullptr)
                continue;
            
            const Mdouble W = interaction->getAbsoluteNormalForce();
            // Pythagoras to get tangential magnitude from the normal magnitude and the relative velocity vector.
            const Mdouble tangentialRelativeVelocity = std::sqrt(interaction->getRelativeVelocity().getLengthSquared() -
                    interaction->getNormalRelativeVelocity() * interaction->getNormalRelativeVelocity());
            const Mdouble L = tangentialRelativeVelocity * getHandler()->getDPMBase()->getTimeStep();
            const Mdouble Q = wearAcceleration_ * wearCoefficient_ * W * L / hardness_;
            
            totalDebris += Q;
            // The contact point in normal direction is halfway the overlap. However, to properly store the debris the
            // position exactly on the wall should be used.
            Vec3D contactPointOnWall = interaction->getContactPoint() + 0.5 * interaction->getOverlap() * interaction->getNormal();
            // Debris is removed in the direction of the interaction normal (rotated back to lab frame).
            storeDebris(t, getOrientation().rotateBack(contactPointOnWall - getPosition()), getOrientation().rotateBack(interaction->getNormal()) * -Q, debrisContainer);
        }
    }
    
    moveVerticesToMatchVolume(debrisContainer, totalDebris);
}

void WearableTriangleMeshWall::storeDebris(const Triangle& triangle, const Vec3D& position, const Vec3D& debris, std::vector<Vec3D>& debrisContainer)
{
    const unsigned index0 = triangle.vertexIndices[0];
    const unsigned index1 = triangle.vertexIndices[1];
    const unsigned index2 = triangle.vertexIndices[2];
    
//    // Assign 1/3 to each vertex
//    debrisContainer[index0] += debris / 3.0;
//    debrisContainer[index1] += debris / 3.0;
//    debrisContainer[index2] += debris / 3.0;
    

//    // Distance from contact point to each vertex
//    const Mdouble distanceSquared0 = Vec3D::getLengthSquared(vertices_[index0].position - position);
//    const Mdouble distanceSquared1 = Vec3D::getLengthSquared(vertices_[index1].position - position);
//    const Mdouble distanceSquared2 = Vec3D::getLengthSquared(vertices_[index2].position - position);
//
//    // Assign to closest vertex
//    const unsigned indexFinal = distanceSquared0 < distanceSquared1 ?
//            (distanceSquared0 < distanceSquared2 ? index0 : index2) :
//            (distanceSquared1 < distanceSquared2 ? index1 : index2);
//    debrisContainer[indexFinal] += debris;

    
    // Barycentric interpolation
    // Get vertices positions and shift by debris position.
    const Vec3D v0 = vertices_[index0].position - position;
    const Vec3D v1 = vertices_[index1].position - position;
    const Vec3D v2 = vertices_[index2].position - position;
    
    // Area of each section, and total area.
    // When the debris position does not lie perfectly in plane, the sum of the area of the sections differs from the
    // area of the triangle itself. It is therefore better to use the sum, even though still a small error in the
    // fractions remain, depending on how far off plane the position is.
    Mdouble A0 = 0.5 * Vec3D::cross(v0, v1).getLength();
    Mdouble A1 = 0.5 * Vec3D::cross(v1, v2).getLength();
    Mdouble A2 = 0.5 * Vec3D::cross(v2, v0).getLength();
    Mdouble A = A0 + A1 + A2;
    
    // Fraction of debris is equal to the  area opposite to vertex divided by total area.
    debrisContainer[index0] += debris * A1 / A;
    debrisContainer[index1] += debris * A2 / A;
    debrisContainer[index2] += debris * A0 / A;
}

void WearableTriangleMeshWall::setWearCoefficient(Mdouble wearCoefficient)
{
    wearCoefficient_ = wearCoefficient;
}

void WearableTriangleMeshWall::setHardness(Mdouble hardness)
{
    hardness_ = hardness;
}

void WearableTriangleMeshWall::setWearAcceleration(Mdouble wearAcceleration)
{
    wearAcceleration_ = wearAcceleration;
}