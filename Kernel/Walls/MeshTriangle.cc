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

#include <MercuryBase.h>
#include <algorithm>

#include "MeshTriangle.h"
#include "InteractionHandler.h"
#include "WallHandler.h"
#include "DPMBase.h"
#include "Particles/BaseParticle.h"

/*!
 * \param[in] p Pointer to the BaseParticle which we want to check the interaction for.
 * \param[in] timeStamp The time at which we want to look at the interaction.
 * \param[in] interactionHandler A pointer to the InteractionHandler in which the interaction can be found.
 * \return A pointer to the BaseInteraction that happened between this BaseWall
 * and the BaseParticle at the timeStamp.
 *
 * \details It is determined, if there is any contact between the particle and
 * the triangular wall. If a contact exists, all neccesary quantities are determined
 * and set. This includes the contact type is determined (corner, edge or 
 * face contact, see laos doi:10.1002/nme.4487) and set with setMultiContactIdentifier.
 */
BaseInteraction*
MeshTriangle::getInteractionWith(BaseParticle* p, unsigned timeStamp, InteractionHandler* interactionHandler)
{
    if (!isActive) return nullptr;
    
    // TODO: Note that this if may lead to contacts beeing ignored in md->checkParticleForInteraction because unadded
    // particles have the Id 0, which might indeed be the Id of an node particle
    if ( std::find(vertexIds_.begin(), vertexIds_.end(), p->getId()) != vertexIds_.end() )
    {
        // Do not detect any contact between particles that correspond to the walls nodes
        return nullptr;
    }
    
    Mdouble distance;
    Vec3D normal;
    Mdouble overlap;
    unsigned int type;

    if (!(p->isSphericalParticle()))
    {
        logger(ERROR, "MeshTriangle::getInteractionWith not implemented for particles of type %", p->getName());
    }

    if (getDistanceNormalOverlapType(*p, distance, normal, overlap, type))
    {
        // look for an existing interaction, or create a new one
        BaseInteraction *c = nullptr;
        if (getGroupId() > 0 && p->getInteractions().size() > 0)
        {
            // Do not care for the result in that case.
            return BaseWall::getInteractionWith(p, timeStamp, interactionHandler);
        }

        // look for an existing interaction, or create a new one
        c = interactionHandler->getInteraction(p, this, timeStamp);

        c->setNormal(-normal);
        c->setDistance(distance);
        c->setOverlap(overlap);
        c->setMultiContactIdentifier(type);
        c->setWallInteraction(1);

        ///\todo{DK: What is the contact point for interactions with walls}
        c->setContactPoint(p->getPosition() - (p->getRadius() - 0.5 * c->getOverlap()) * c->getNormal());

        logger(DEBUG, "Particle contact with wall at %", c->getContactPoint());
        return c;
    }
    return nullptr;
}

/*!
 * \param[in] interactionHandler Pointer to InteractionHandler that contains all the interactions.
 * \param[in] timeStamp The time at which we want to look at the interactions.
 *
 * \details Determine, if a contact os valid, i.e. that the forces due to this 
 * contact should be applied to both the particle in the wall. The evaluation
 * is done by looking at the interaction a specific particle has with the current
 * wall as well as neighboring walls. If a particle has multiple contacts,
 * the selection criteria noted in doi:10.1002/nme.4487 are applied.
 * Note: At this time this leads to issues when the particles are much bigger than
 * the triangles.
 */
void MeshTriangle::checkInteractions(InteractionHandler* interactionHandler, unsigned int timeStamp)
{
    unsigned j, id;
    // Iterate through the interactions and check if there is one with higher priority
    // Note: The id should already be taken into account when creating the interactions
    for (j = 0; j<getInteractions().size(); j++)
    {
        bool interactionValid = true;
        auto i = getInteractions()[j];
        if (i->getTimeStamp()!=timeStamp)
        {
            continue; // Old interactions are not needed
            logger(WARN, "Taking old interaction into account");
        }
        if (i->getMultiContactIdentifier()>3) // Vertex contact
        {
            id = i->getMultiContactIdentifier() - 4;
            
            for (auto triangleId : vertexNeighbors[id])
            {
                // If a neighboring particle has a contact and a lower Id with
                // the same particle, disregard this contact
                if (triangleId<this->getId())
                {
                    interactionValid = false;
                    break;
                }
                // logger(INFO, "id % %", id, getId());
                auto k = getHandler()->getObjectById(triangleId);
                auto interactionNeighbor = interactionHandler->getExistingInteraction(i->getP(), k);
                if (interactionNeighbor && interactionNeighbor->getTimeStamp()==timeStamp)
                {
                    if (interactionNeighbor->getMultiContactIdentifier() <= 3)
                    {
                        // logger(INFO, "Found invalid vertex contact");
                        interactionValid = false;
                        break;
                    }
                }
            }
        }
        else if (i->getMultiContactIdentifier()>0) // Edge contact
        {
            id = i->getMultiContactIdentifier() - 1;
            if (neighbor[id])
            {
                // If a neighboring particle has at least an edge contact with
                // the same particle and a lower Id, disregard this contact
                if (neighbor[id]->getId()<this->getId())
                {
                    // logger(INFO, "Found invalid Edge contact due to Id");
                    interactionValid = false;
                } 
                else
                {
                    auto interactionNeighbor = interactionHandler->getExistingInteraction(i->getP(), neighbor[id]);
                    if (interactionNeighbor && interactionNeighbor->getTimeStamp()==timeStamp)
                    {
                        if (interactionNeighbor->getMultiContactIdentifier() == 0)
                        {
                            // logger(INFO, "Found invalid edge contact");
                            interactionValid = false;
                        }
                    }
                }
            }
        }

        // Undo the interaction
        if (!interactionValid)
        {
            i->getP()->addForce(-i->getForce());
            this->addForce(i->getForce());
            if (getHandler()->getDPMBase()->getRotation())
            {
                i->getP()->addTorque(-i->getTorque() + Vec3D::cross(i->getP()->getPosition() - i->getContactPoint(), i->getForce()));
                this->addTorque(i->getTorque() - Vec3D::cross(this->getPosition() - i->getContactPoint(), i->getForce()));
            }
            i->setForce(Vec3D(0,0,0));
        }
        else
        {
            Vec3D weight = getBaricentricWeight(i->getContactPoint());
            
            // Contact seems valid, distribute its force to the edges
            Vec3D force = -i->getForce();

            for(unsigned k=0; k<3; k++)
            {   
                vertexParticle_[k]->addForce(weight.getComponent(k)*force);
            }
        }
    }
}
/*!
 * \param[in] p BaseParticle we want to calculate the distance and whether it collided of.
 * \param[out] distance The distance of the BaseParticle to this wall.
 * \param[out] normal_return If there was a collision, the normal vector to this wall will be placed here.
 * \return A boolean which says whether or not there was a collision.
 * \details This function computes whether or not there is a collision between
 * a given BaseParticle and this MeshTriangle. If there is a collision, this
 * function also computes the distance between the BaseParticle and MeshTriangle
 * and the normal of the MeshTriangle at the intersection point. It does
 * this by calling MeshTriangle::getDistanceNormalOverlapType.
 * Since this function should be called before calculating any
 * Particle-Wall interactions, it can also be used to set the normal vector in
 * case of curved walls.
 */
bool MeshTriangle::getDistanceAndNormal(const BaseParticle& p, Mdouble& distance, Vec3D& normal) const
{
    if (!isActive) return false;
    Mdouble overlap;
    unsigned int type;
    return getDistanceNormalOverlapType(p, distance, normal, overlap, type);
}


/*!
 * \param[in] p BaseParticle we want to calculate the distance and whether it collided of.
 * \param[out] distance The distance of the BaseParticle to this wall.
 * \param[out] normal_return If there was a collision, the normal vector to this wall will be placed here.
 * \param[out] overlap_return If there was a collision, the overlap will be placed here.
 * \param[out] type_return If there was a collision, the contact type will be placed here.
 * \return A boolean which says whether or not there was a collision.
 * \details This function computes whether or not there is a collision between
 * a given BaseParticle and this MeshTriangle. If there is a collision, this
 * function also computes the distance between the BaseParticle and MeshTriangle
 * and the normal of the MeshTriangle at the intersection point as well as the
 * contact overlap and type.
 * The type is set to the following:
 * 0: face contact
 * 1, 2, 3: contact with type-1th edge
 * 4, 5, 6; contact with the type-4th vertex
 */
bool MeshTriangle::getDistanceNormalOverlapType(const BaseParticle& p, Mdouble& distance, Vec3D& normal, Mdouble& overlap, unsigned int& type) const
{
    if (!isActive) return false;
    
    // TODO: Note that this if may lead to contacts beeing ignored in md->checkParticleForInteraction because unadded
    // particles have the Id 0, which might indeed be the Id of a node particle
    if ( std::find(vertexIds_.begin(), vertexIds_.end(), p.getId()) != vertexIds_.end() )
    {
        // Do not detect any contact between particles that correspond to the walls nodes
        return false;
    }
    
    const Vec3D position = p.getPosition(); // note: no transfer to lab coordinates
    const Mdouble distanceMax = p.getWallInteractionRadius(this);

    // compute distance from face
    // get distance from particle position to the face
    const Mdouble signedDistance = Vec3D::dot(position-vertex_[0], faceNormal_);
    distance = std::abs(signedDistance);

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
        overlap = p.getRadius() - distance;
        type=0; // Face contact
        return true;
    }

    // Then the neighbor will handle this interaction
    // determine if it is an edge or vertex contact
    const double positionAlongEdge = Vec3D::dot(edgeBranch[id], edge_[id]);
    if (positionAlongEdge <= 0) {
        //possible contact with left vertex
        distance = edgeBranch[id].getLength();
        if (distance > distanceMax) return false;
        // check vertex ids
        normal = edgeBranch[id] / -distance;
        type = 4 + id; // Vertex contact
    } else if (positionAlongEdge >= edgeLength_[id]) {
        //contact with right vertex
        const size_t idRight = (id + 1) % 3;
        distance = edgeBranch[idRight].getLength();
        if (distance > distanceMax) return false;
        // check vertex ids
        type = 4 + idRight; // Vertex contact
        normal = edgeBranch[idRight] / -distance;

    } else {
        // edge contact
        normal = edge_[id] * positionAlongEdge - edgeBranch[id];
        distance = normal.getLength();
        if (distance > distanceMax) return false;
        normal /= distance;
        type = 1 + id;
    }

    overlap = p.getRadius() - distance;
    return true;
}

/*!
 * \param[in] contact The point, at which the velocity should be determined.
 * \return A Vec3D giving the calculated velocity.
 * \details Calculates the velocity at the contact position by interpolating
 * the velocity of the triangle nodes using barycentric coordinates. 
 */
const Vec3D MeshTriangle::getVelocityAtContact(const Vec3D& contact) const
{
    Vec3D m = getBaricentricWeight(contact);
    return m.x()*vertexVelocity_[0] + m.y()*vertexVelocity_[1] + m.z()*vertexVelocity_[2];
};

/*!
 * \param[in] contact The point, at which the weights should be calculated.
 * \return A Vec3D giving the weights. (x is the weight of corner 0,...)
 * \details Calculates baricentric weight of a given point in the triangle
 */
const Vec3D MeshTriangle::getBaricentricWeight(const Vec3D& contact) const
{
    Vec3D m;
    // The area of a triangle is
    // Taken from https://gamedev.stackexchange.com/questions/23743/whats-the-most-efficient-way-to-find-barycentric-coordinates
    // Mdouble areaABC = Vec3D::dot( getFaceNormal(), Vec3D::cross(vertex_[1]-vertex_[0], vertex_[2]-vertex_[0] )  ) ;
    // logger.assert_always(std::fabs(areaABC-area_*2)<1e-8, "OOPS, triangle area calculared wrong % %", areaABC, area_*2);
    Mdouble areaPBC = 0.5*Vec3D::dot( getFaceNormal(), Vec3D::cross(vertex_[1]-contact, vertex_[2]-contact )  ) ;
    Mdouble areaPCA = 0.5*Vec3D::dot( getFaceNormal(), Vec3D::cross(vertex_[2]-contact, vertex_[0]-contact )  ) ;

    m.X = areaPBC / area_ ; // alpha
    m.Y = areaPCA / area_ ; // beta
    m.Z = 1.0f - m.X - m.Y ; // gamma
    return m;
}

void MeshTriangle::rotate(const Vec3D& angularVelocityDt)
{
    if (!angularVelocityDt.isZero())
    {
        BaseInteractable::rotate(angularVelocityDt);
        updateVertexAndNormal();
    }
}

/*!
 * \param[in] is The input stream from which the MeshTriangle is read, usually a restart file.
 */
void MeshTriangle::read(std::istream& is)
{
    BaseWall::read(is);
    
    unsigned int i, j, s1, s2, Id;
    std::string dummy;
    
    is >> dummy >> s1;
    vertexNeighbors.reserve(s1);
    for (i=0; i<s1; i++)
    {
        is >> dummy >> s2;
        std::vector<unsigned int> vN;
        vN.reserve(s2);
        
        for (j=0; j<s2; j++)
        {
            is >> Id;
            vN.push_back(Id);
        }
        vertexNeighbors.push_back(vN);
        
    }
    
    
    
    is >> dummy;
    for (int i = 0; i < 3; i++)
    {
        is >> vertex_[i];
    }
    is >> dummy;
    for (int i = 0; i < 3; i++)
    {
        is >> vertexIds_[i];
    }
    
    is >> dummy >> invMass_;
    is >> dummy >> isActive;
    
    // updateVerticesFromParticles();
}

/*!
 * \param[in] os The output stream where the MeshTriangle must be written
 *  to, usually a restart file.
 */
void MeshTriangle::write(std::ostream& os) const
{
    BaseWall::write(os);
    unsigned int i, j, s1, s2;
    s1 = vertexNeighbors.size();
    os << " vertexNeighbors " << s1;
    for (i=0; i<s1; i++)
    {
        s2 = vertexNeighbors[i].size();
        os << " vertexNeighborsElements " << s2;
        for (j=0; j<s2; j++)
        {
            os << " " << vertexNeighbors[i][j];            
        }
    }
    
    os << " vertex ";
    for (int i = 0; i < 3; i++)
    {
        os << ' ' << vertex_[i];
    }
    os << " edgeParticleIds ";
    for (int i = 0; i < 3; i++)
    {
        os << ' ' << vertexIds_[i];
    }
    
    os << " invMass " << invMass_;
    
    os << " isActive " << isActive;
    
    // Note: need not to write faceNormal_ and area_, as these are recalculated
    // on read.
}

void MeshTriangle::writeVTK(VTKContainer& vtk) const
{
    if (!isActive) return;
    
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
    
    vtk.triangleStrips.push_back(cell);
}

void MeshTriangle::setVertices(const Vec3D A, const Vec3D B, const Vec3D C)
{
    setVertices(A, B, C, (A + B + C) / 3);
}

void MeshTriangle::setVertices(const Vec3D A, const Vec3D B, const Vec3D C, const Vec3D position)
{
    setPosition(position);
    setOrientation({1, 0, 0, 0});
    vertex_[0] = A;
    vertex_[1] = B;
    vertex_[2] = C;
    updateVertexAndNormal();
}

void MeshTriangle::setVertexVelocities(const Vec3D A, const Vec3D B, const Vec3D C)
{
    vertexVelocity_[0] = A;
    vertexVelocity_[1] = B;
    vertexVelocity_[2] = C;
}

/*!
 * \param[in] i, j, k The ids of the particles representing the corners.
 */
void MeshTriangle::setVertexIds(unsigned int i, unsigned int j, unsigned int k)
{
    vertexIds_[0] = i;
    vertexIds_[1] = j;
    vertexIds_[2] = k;
    retrieveVertexParticles();
}

/*!
 * \param[in] handler Pointer to the wallHandler.
 * \details Sets the wall handler and calls retrieveVertexParticles
 */
void MeshTriangle::setHandler(WallHandler* handler)
{
    BaseWall::setHandler(handler);
    
    if (handler->getDPMBase()->particleHandler.getNumberOfObjects()==0)
    {
        // Without this, restart does not work
        return;
    }

    retrieveVertexParticles();
}

/*!
 * \details Update the triangle position and velocity based on the vertex particles,
 * if the triangle is active.
 */
void MeshTriangle::updateVerticesFromParticles()
{
    // Need to get references to the particles
    checkActive();
    if (!isActive) return;
    
    setVertices(vertexParticle_[0]->getPosition(), vertexParticle_[1]->getPosition(), vertexParticle_[2]->getPosition());
    setVertexVelocities(vertexParticle_[0]->getVelocity(), vertexParticle_[1]->getVelocity(), vertexParticle_[2]->getVelocity());
}

/*!
 * \param[in] id The id of the removed particle
 * \details If the given id is equal to one of the vertexParticles, the reference
 * to that particle is removed and the triangle is set inactive
 */
void MeshTriangle::handleParticleRemoval(unsigned int id)
{
    unsigned int i;
    for (i=0; i<3; i++)
    {
        if (vertexIds_[i] == id)
        {
            vertexParticle_[i] = nullptr;
            isActive = false;
        }
    }
}

/*!
 * \param[in] id The id of the added particle
 * \param[in] p Pointer to the particle
 * \details If thie given id is equal to one of the vertex Particles, the reference
 * is added. If this added particle makes the triangle active, the vertex positions
 * are uodated using updateVerticesFromParticles().
 */
void MeshTriangle::handleParticleAddition(unsigned int id, BaseParticle* p)
{
    unsigned int i;
    for (i=0; i<3; i++)
    {
        if (vertexIds_[i] == id)
        {
            vertexParticle_[i] = p;
            checkActive();
            updateVerticesFromParticles();
        }
    }
}

/*!
 * \details Tries to get the pointer to the vertex particles from the particleHandler.
 * Afterwards it checks if the triangle is active and updates the vertex positions.
 */
void MeshTriangle::retrieveVertexParticles()
{
    if (getHandler())
    {
        unsigned int i;
        for (i=0; i<3; i++)
        {
            vertexParticle_[i] = getHandler()->getDPMBase()->particleHandler.getObjectById(vertexIds_[i]);
        }
        checkActive();
        updateVerticesFromParticles();
    }
}

/*!
 * \details Checks if the triangle is active. A triangle is considered active,
 * if pointers to all references are known.
 */
void MeshTriangle::checkActive()
{
    isActive = !(!vertexParticle_[0] || !vertexParticle_[1] || !vertexParticle_[2]);
}

/*!
 * \details On restart, try to get all vertex particle pointers.
 */
void MeshTriangle::actionsOnRestart()
{
    // Note this can not be done in the read sequence, as the particles are not yet available there
    retrieveVertexParticles();
}

/*!
 * \details After the particles get new positions, these need to be retrived to 
 * update the triangle.
 */
void MeshTriangle::actionsAfterParticleGhostUpdate()
{
    updateVerticesFromParticles();
}

/*!
 * \details Moves (displaces) the interacable a given distance.
 *          Note, this just updates the position by the move.
 * \param[in] move  Reference to Vec3D which is the distance to move the
 *            interactable.
 */
void MeshTriangle::move(const Vec3D& move)
{
    BaseInteractable::move(move);
    updateVertexAndNormal();
}

/**
 * This function should be called after setting either position_ or vertexInLabFrame_.
 *  - vertex is set to the vertex position in the real coordinate system (rotated and shifted)
 *  - vertexMin_/vertexMax_ is set to define a bounding box around the wall (for contact detection)
 *  - edge_, edgeNormal_ and faceNormal_ is updated (stored for quick computation of contact point)
 */
void MeshTriangle::updateVertexAndNormal()
{
    vertexMin_ = Vec3D::min(Vec3D::min(vertex_[0], vertex_[1]), vertex_[2]);
    vertexMax_ = Vec3D::max(Vec3D::max(vertex_[0], vertex_[1]), vertex_[2]);

    edge_ = {vertex_[1] - vertex_[0], vertex_[2] - vertex_[1], vertex_[0] - vertex_[2]};
    faceNormal_ = Vec3D::cross(edge_[0], edge_[1]);
    area_ = 0.5*faceNormal_.getLength();
    logger.assert_always(0.5*sqrt(Vec3D::cross(vertex_[1] - vertex_[0], vertex_[2] - vertex_[1]).getLengthSquared())==area_, "OOPS, face area wrong");
    
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

/*!
 * \param[out] min the minimum of the triangle in all directions.
 * \param[out] max the maximum of the triangle in all directions
 * \return true The triangle is local
 */
bool MeshTriangle::isLocal(Vec3D& min, Vec3D& max) const
{
    min = vertexMin_;
    max = vertexMax_;
    return true;
}

/*!
 * \param[in] point Position to check
 * \return boolean specifying if point is within the triangle
 */
bool MeshTriangle::isInsideTriangle(const Vec3D &point) const
{
    Vec3D weights = getBaricentricWeight(point);
    
    Mdouble eps = 1e-12;
    return ((1-eps) > weights.X > eps && (1-eps) > weights.Y > eps && (1-eps) > weights.Z > eps) && ((1-eps) > weights.X+weights.Y > eps && (1-eps) > weights.Y+weights.Z > eps && (1-eps) > weights.Z+weights.X > eps);
    // return ((1-eps) > s0 > eps && (1-eps) > s1 > eps && (1-eps) > s2 > eps) && ((1-eps) > s0+s1 > eps && (1-eps) > s1+s2 > eps && (1-eps) > s2+s0 > eps);
    
}


Mdouble MeshTriangle::getInvMass() const
{
    return invMass_;
}

/*!
 * \param[in] mass Value of the mass assigned to the triangle.
 * \details The mass is neccesary for contact force determination. If not given,
 * infinite mass is assumed.
 */
void MeshTriangle::setMass(Mdouble mass)
{
    logger.assert_always(mass > 0.0,
                        "Error in MeshTriangle::setMass, the given mass to be set must be positive.");
    invMass_ = 1.0 / mass;

}

/*!
 * \param[in] pressure The pressure value in units of 1 Pa
 * \details Calculates the force acting on the triangle using the pressure and 
 * the triangles surface area. Applies the force to the vertex particles.
 */
void MeshTriangle::applyPressure(Mdouble pressure)
{
    if (isActive)
    {
        // Area in Normaldirection
        // Calculate F = p*A/3.0
        // The division by 3.0 is to split the force evenly between the points
        Vec3D pressureForce = getArea()*pressure*getFaceNormal()/3.0;
        for (unsigned int j=0; j<3; j++){
            vertexParticle_[j]->addForce(pressureForce);
        }
    }
}

/*!
 * \param[in] force
 * \details Applies a given force to the triangle, by splitting it between the
 * vertex particles.
 */
void MeshTriangle::applyForce(Vec3D force)
{
    if (isActive)
    {
        force /= 3.0;
        for (unsigned int j=0; j<3; j++){
            vertexParticle_[j]->addForce(force);
        }
    }
    
}
