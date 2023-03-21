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

#include "TriangleMeshWall.h"
#include "WallHandler.h"
#include "DPMBase.h"

TriangleMeshWall::TriangleMeshWall(const TriangleMeshWall& other) : BaseWall(other)
{
    triangles_ = other.triangles_;
    vertices_ = other.vertices_;
    periodicCompanions_ = other.periodicCompanions_;
    
    boundingBoxMin_ = other.boundingBoxMin_;
    boundingBoxMax_ = other.boundingBoxMax_;
    boundingBoxMinGlobal_ = other.boundingBoxMinGlobal_;
    boundingBoxMaxGlobal_ = other.boundingBoxMaxGlobal_;
}

TriangleMeshWall::TriangleMeshWall(const std::vector<Vec3D>& points,
                                   const std::vector<std::array<unsigned, 3>>& cells, const ParticleSpecies* species)
{
    // Set the species of the main wall first, as the set function will use it.
    setSpecies(species);
    set(points, cells);
}

TriangleMeshWall::TriangleMeshWall(const Vec3D& P0, const Vec3D& P1, const Vec3D& P2, Mdouble resolutionU,
                                   Mdouble resolutionV, const ParticleSpecies* species,
                                   bool periodicInU, bool periodicInV)
{
    // Set the species first, as it is used in other methods.
    setSpecies(species);
    createParallelogramMesh(P0, P1, P2, resolutionU, resolutionV, periodicInU, periodicInV);
}

TriangleMeshWall& TriangleMeshWall::operator=(const TriangleMeshWall& other)
{
    if (this == &other)
    {
        return *this;
    }
    return *(other.copy());
}

TriangleMeshWall * TriangleMeshWall::copy() const
{
    return new TriangleMeshWall(*this);
}

void TriangleMeshWall::read(std::istream& is)
{
    BaseWall::read(is);
    
    std::string dummy;
    int num;
    is >> dummy >> num;
    
    std::vector<Vec3D> points(num);
    for (int i = 0; i < num; i++)
    {
        Vec3D p;
        is >> p;
        points[i] = p;
    }
    
    is >> dummy >> num;
    std::vector<std::array<unsigned, 3>> cells(num);
    for (int i = 0; i < num; i++)
    {
        std::array<unsigned, 3> c{};
        is >> c[0] >> c[1] >> c[2];
        cells[i] = c;
    }
    
    is >> dummy >> num;
    std::vector<std::pair<unsigned, unsigned>> periodicCompanions(num);
    for (int i = 0; i < num; i ++)
    {
        std::pair<unsigned, unsigned> p;
        is >> p.first >> p.second;
        periodicCompanions[i] = p;
    }
    
    set(points, cells);
    setPeriodicCompanions(periodicCompanions);
    
    updateBoundingBox();
}

void TriangleMeshWall::write(std::ostream& os) const
{
    BaseWall::write(os);
    
    os << " points " << vertices_.size();
    for (const Vertex& v : vertices_)
        os << ' ' << v.position;
    os << " cells " << triangles_.size();
    for (const Triangle& t : triangles_)
        os << ' ' << t.vertexIndices[0] << ' ' << t.vertexIndices[1] << ' ' << t.vertexIndices[2];
    os << " periodicCompanions " << periodicCompanions_.size();
    for (auto& p : periodicCompanions_)
        os << ' ' << p.first << ' ' << p.second;
}

std::string TriangleMeshWall::getName() const
{
    return "TriangleMeshWall";
}

void TriangleMeshWall::set(const std::vector<Vec3D>& points, const std::vector<std::array<unsigned, 3>>& cells)
{
    // The species for the triangles are the same as the main wall.
    // Note: when the species have not been set yet, this will be a nullptr,
    // so a check is done before actually setting the species for the triangles.
    auto species = getSpecies();
    
    // Clear mesh, as set function can be called multiple times to create different meshes. Useful when creating a
    // complicated mesh with the addToMesh() function.
    triangles_.clear();
    vertices_.clear();
    triangles_.reserve(cells.size());
    vertices_.reserve(points.size());
    
    // Add all the points as Vertex to the vector.
    for (auto& p : points)
        vertices_.emplace_back(p);
    
    for (auto& c : cells)
    {
        TriangleWall tw;
        // The position is the origin of the local coordinate system, which is always (0, 0, 0) (not getPosition()!).
        // This is needed in case the mesh has an angular velocity.
        tw.setVertices(points[c[0]], points[c[1]], points[c[2]], Vec3D(0.0, 0.0, 0.0));
        if (species != nullptr)
            tw.setSpecies(species);
        
        // Add as Triangle with TriangleWall and vertexIndices to the vector.
        triangles_.emplace_back(tw, c);
        
        // Add the Triangle index and the local vertex index to each of the Vertices.
        const unsigned tIndex = triangles_.size() - 1;
        vertices_[c[0]].triangleIndices.emplace_back(tIndex, 0);
        vertices_[c[1]].triangleIndices.emplace_back(tIndex, 1);
        vertices_[c[2]].triangleIndices.emplace_back(tIndex, 2);
    }
    
    updateBoundingBox();
}

bool TriangleMeshWall::getDistanceAndNormal(const BaseParticle& p, Mdouble& distance, Vec3D& normal_return) const
{
    // Input BaseParticle is constant, so copy the necessary properties and change to local coordinates.
    SphericalParticle pLocal;
    pLocal.setSpecies(p.getSpecies());
    pLocal.setRadius(p.getRadius());
    pLocal.setPosition(getOrientation().rotateBack(p.getPosition() - getPosition()));
    
    if (!isWithinBoundingBox(pLocal.getPosition(), pLocal.getWallInteractionRadius(this)))
        return false;
    
    for (const Triangle& t : triangles_)
    {
        if (t.wall.getDistanceAndNormal(pLocal, distance, normal_return))
            return true;
    }
    return false;
}

BaseInteraction* TriangleMeshWall::getInteractionWith(BaseParticle* p, unsigned int timeStamp,
                                                       InteractionHandler* interactionHandler)
{
    // In case of multiple contacts, we store the one with the maximum overlap and only calculate the forces for that interaction.
    // Since the BaseWall::getInteractionWith() automatically adds interactions to the handler, any interactions which
    // are not the one with greatest overlap should be removed from the handler.
    Mdouble maxOverlap = 0.0;
    BaseWall* maxWall = nullptr;
    
    // Change particle to local coordinates.
    //\todo How to handle orientation of non-spherical particles
    const Vec3D posOriginal = p->getPosition();
    p->setPosition(getOrientation().rotateBack(posOriginal - getPosition()));
    
    if (!isWithinBoundingBox(p->getPosition(), p->getWallInteractionRadius(this)))
    {
        p->setPosition(posOriginal);
        return nullptr;
    }
    
    // Loop through all triangles and find interactions and store the one with max overlap
    // (the main idea and the force calculations are literally copied from DPMBase::computeForcesDueToWalls)
    for (auto& t : triangles_)
    {
        // Checks if the particle is interacting with this triangle
        BaseInteraction* i = t.wall.getInteractionWith(p, timeStamp, interactionHandler);
        
        // When an interaction exists
        if (i != nullptr)
        {
            // When the overlap is greater than that of any previous interactions
            if (i->getOverlap() > maxOverlap)
            {
                // Remove the previous maximum overlap interaction from the handler
                if (maxWall != nullptr)
                {
                    BaseInteraction* prevMaxI = interactionHandler->getInteraction(p, maxWall, timeStamp);
                    interactionHandler->removeObject(prevMaxI->getIndex());
                }
                
                // Update the maximum overlap and store the interacting wall
                maxOverlap = i->getOverlap();
                maxWall = &t.wall;
            }
            else
            {
                // Remove the last added interaction, as it is not the one with maximum overlap and should be discarded.
                interactionHandler->removeObject(i->getIndex());
            }
        }
    }
    
    // Change particle back to its original coordinates, before doing any force calculations.
    p->setPosition(posOriginal);
    
    // When an interaction exists
    if (maxWall != nullptr)
    {
        // The interaction with maximum overlap
        BaseInteraction* i = interactionHandler->getInteraction(p, maxWall, timeStamp);
        // The TriangleWall which with the interaction occurred
        BaseInteractable* w = maxWall;
        
        // Rotate normal and change contact point to global coordinates
        i->setNormal(getOrientation().rotate(i->getNormal()));
        i->setContactPoint(getOrientation().rotate(i->getContactPoint()) + getPosition());
        
        // Transform TriangleWall position from local to global coordinates. (Needed for torque calculations etc.)
        const Vec3D wPosOriginal = w->getPosition();
        w->setPosition(getOrientation().rotateBack(w->getPosition()) + getPosition());
        // Add the mesh velocities to the TriangleWall velocities. Needed in the force calculations (relative velocity etc.).
        w->addVelocity(getVelocity());
        w->addAngularVelocity(getAngularVelocity());
        // At the moment the TriangleWalls in the mesh always have position, velocity and angular velocities equal to zero.
        // In case of future changes, however, they are treated as if they do have an actual value and therefore at the
        // end are also reset to those original values.
        
        //...calculates the forces between the two objects...
        i->computeForce();
        
        //...and applies them to each of the two objects (wall and particle).
        p->addForce(i->getForce());
        w->addForce(-i->getForce());
        
        //If the rotation flag is on, also applies the relevant torques (getRotation() returns a boolean).
        if (getHandler()->getDPMBase()->getRotation())
        {
            p->addTorque(i->getTorque() - Vec3D::cross(p->getPosition() - i->getContactPoint(), i->getForce()));
            ///\todo TW: I think this torque has the wrong sign
            w->addTorque(-i->getTorque() + Vec3D::cross(w->getPosition() - i->getContactPoint(), i->getForce()));
        }
        
        // Set the TriangleWall position and velocities back to their original (local) values.
        w->setPosition(wPosOriginal);
        w->addVelocity(-getVelocity());
        w->addAngularVelocity(-getAngularVelocity());
    }
    
    // Return no interaction for the main wall itself.
    return nullptr;
}

void TriangleMeshWall::resetForceTorque(int numberOfOMPthreads)
{
    BaseWall::resetForceTorque(numberOfOMPthreads);
    
    for (auto& t : triangles_)
        t.wall.resetForceTorque(numberOfOMPthreads);
}

void TriangleMeshWall::writeVTK(VTKContainer& vtk) const
{
    const unsigned long s = vtk.points.size();

    for (auto& v : vertices_)
    {
        vtk.points.push_back(getOrientation().rotate(v.position) + getPosition());
    }
    
    for (auto& t : triangles_)
    {
        vtk.triangleStrips.push_back({ static_cast<double>(s + t.vertexIndices[0]),
                                       static_cast<double>(s + t.vertexIndices[1]),
                                       static_cast<double>(s + t.vertexIndices[2]) });
    }
}

void TriangleMeshWall::setSpecies(const ParticleSpecies* species)
{
    if (species == nullptr)
        return;
    BaseWall::setSpecies(species);
    for (auto& t : triangles_)
        t.wall.setSpecies(species);
}

void TriangleMeshWall::moveVertex(const unsigned index, const Vec3D& dP)
{
    // Update position of Vertex in vector.
    vertices_[index].position += dP;
    
    // Update each Triangle wall of the triangles who share this Vertex.
    for (auto& ti : vertices_[index].triangleIndices)
        triangles_[ti.first].wall.moveVertex(ti.second, dP);
    
    // Do the same for each of the possible companions of this index.
    for (unsigned idx : getPeriodicCompanions(index))
    {
        // Update position of Vertex in vector.
        vertices_[idx].position += dP;
    
        // Update each Triangle wall of the triangles who share this Vertex.
        for (auto& ti : vertices_[idx].triangleIndices)
            triangles_[ti.first].wall.moveVertex(ti.second, dP);
    }
    
    updateBoundingBox();
}

void TriangleMeshWall::removeTriangle(const unsigned index, bool removeFreeVertex)
{
    // Remove any interactions the TriangleWall might have.
    // Otherwise later when in DPMBase interactionHandler.eraseOldInteractions(getNumberOfTimeSteps()) is called it
    // might cause crashes, since it might be trying to delete pointers which are already used for something else.
    for (BaseInteraction* i : triangles_[index].wall.getInteractions())
        i->removeFromHandler();
    
    // Swap to Triangle with the last Triangle in the vector and then simply pop_back. This is faster, cleaner and
    // simpler then using erase and having to update all other Triangle indices, since every index higher than the
    // removed one would have shifted one down. Now only the swapped Triangle index references have to be updated.
    const unsigned lastIndex = triangles_.size() - 1;
    Triangle triangleToRemove = triangles_[index];
    
    // Don't swap right way, since the swapped Triangle index references need to be updated, but since it's possible
    // that a free Vertex will be removed this might cause wrong index referencing.
    // So first remove the index reference to this Triangle, for each of its Vertices.
    for (unsigned vi : triangleToRemove.vertexIndices)
    {
        // Get the Triangle indices from this Vertex.
        auto& tis = vertices_[vi].triangleIndices;
        // Erase the index of the removed Triangle. This is only a short vector so erase is fine.
        // This will only remove one item, as it should only contain it once.
        tis.erase(std::remove_if(tis.begin(), tis.end(),[index](auto& p) { return p.first == index; }));
        
        // When no other Triangles are connected to this Vertex, it might be erased from the vector.
        if (removeFreeVertex && tis.empty())
        {
            // Removing a Vertex is done in a similar manner, by swapping and popping back.
            const unsigned lastIndexV = vertices_.size() - 1;
            Vertex vertexToRemove = vertices_[vi];
            vertices_[vi] = vertices_[lastIndexV];
            vertices_[lastIndexV] = vertexToRemove;
            vertices_.pop_back();
            
            // No Triangles with index references to the removed Vertex exist.
            // However, for every Triangle which is referencing the swapped Vertex, update the index with the removed index.
            for (auto& ti : vertices_[vi].triangleIndices)
                for (unsigned& tvi : triangles_[ti.first].vertexIndices)
                    if (tvi == lastIndexV)
                    {
                        tvi = vi;
                        break;
                    }
    
            // If this vertex has any periodic companions, those should be removed from the vector.
            periodicCompanions_.erase(std::remove_if(periodicCompanions_.begin(), periodicCompanions_.end(),
                [vi](auto& p) { return p.first == vi || p.second == vi; }), periodicCompanions_.end());
            // Also the swapped index should be updated with the removed index.
            for (auto& p : periodicCompanions_)
            {
                if (p.first == lastIndexV)
                    p.first = vi;
                else if (p.second == lastIndexV)
                    p.second = vi;
            }
        }
    }
    
    // Now swap the Triangles and pop back.
    triangles_[index] = triangles_[lastIndex];
    triangles_[lastIndex] = triangleToRemove;
    triangles_.pop_back();
    
    // For every Vertex which is referencing the swapped Triangle, update the index with the removed index.
    for (unsigned& vi : triangles_[index].vertexIndices)
        for (auto& ti : vertices_[vi].triangleIndices)
            if (ti.first == lastIndex)
            {
                ti.first = index;
                break;
            }
    
    updateBoundingBox();
}

void TriangleMeshWall::refineTriangle(const unsigned index, unsigned numberOfTimes)
{
    // A new vertex will be added in the middle of the triangle and 3 new triangles will be formed. The original triangle
    // will not be deleted, as erasing items from a vector can be slow. Instead its vertices positions and indices will
    // be updated to form one of the smaller triangles.
    // The local vertex indices are 0, 1, 2 and the newly added 3. The original triangle will go from 0, 1, 2 to 0, 1, 3.
    // The first new triangle will have indices 1, 2, 3, and the second 2, 0, 3.
    
    // (1)                                 (2)
    // _____________________________________
    // \\                                //
    //  \  \                          /  /
    //   \   \                     /    /
    //    \     \               /      /
    //     \       \          /       /
    //      \         \     /        /
    //       \          \ /         /
    //        \          |(3)      /
    //         \         |        /
    //          \        |       /
    //           \       |      /
    //            \      |     /
    //             \     |    /
    //              \    |   /
    //               \   |  /
    //                \  | /
    //                 \ |/
    //                  \/
    //                 (0)
    
    // Reserve space for 2 new triangles and 1 new vertex, otherwise references will be broken when adding to vectors.
    triangles_.reserve(triangles_.size() + 2);
    vertices_.reserve(vertices_.size() + 1);
    
    // The triangle to refine.
    Triangle& t = triangles_[index];
    
    // And its vertices.
    Vertex& v0 = vertices_[t.vertexIndices[0]];
    Vertex& v1 = vertices_[t.vertexIndices[1]];
    Vertex& v2 = vertices_[t.vertexIndices[2]];
    
    // Center of triangle is position for new vertex.
    const Vec3D p3 = (v0.position + v1.position + v2.position) / 3.0;
    vertices_.emplace_back(p3);
    const unsigned v3Index = vertices_.size() - 1;
    Vertex& v3 = vertices_[v3Index];
    
    // First new triangle.
    TriangleWall tw1;
    tw1.setVertices(v1.position, v2.position, v3.position);
    tw1.setSpecies(t.wall.getSpecies());
    Triangle t1(tw1, { t.vertexIndices[1], t.vertexIndices[2], v3Index });
    triangles_.push_back(t1);
    
    // Second new triangle.
    TriangleWall tw2;
    tw2.setVertices(v2.position, v0.position, v3.position);
    tw2.setSpecies(t.wall.getSpecies());
    Triangle t2(tw2, { t.vertexIndices[2], t.vertexIndices[0], v3Index });
    triangles_.push_back(t2);
    
    // Their indices for access in the triangles vector.
    const unsigned t2Index = triangles_.size() - 1;
    const unsigned t1Index = t2Index - 1;
    
    // Store the new triangle index and its local vertex index in the vertices.
    v0.triangleIndices.emplace_back(t2Index, 1);
    v1.triangleIndices.emplace_back(t1Index, 0);
    v2.triangleIndices.emplace_back(t1Index, 1);
    // The original triangle does not share vertex 2 anymore, so overwrite its index with that of the second new triangle.
    for (auto& p : v2.triangleIndices)
    {
        if (p.first == index)
        {
            p = std::pair<unsigned, unsigned>(t2Index, 0);
            break;
        }
    }
    
    // Store the triangles indices and local vertex indices in the newly added vertex.
    v3.triangleIndices.emplace_back(index, 2);
    v3.triangleIndices.emplace_back(t1Index, 2);
    v3.triangleIndices.emplace_back(t2Index, 2);
    
    // Overwrite the second vertex index of the original triangle with the newly added vertex.
    t.vertexIndices[2] = v3Index;
    // Update the original triangle vertex positions.
    t.wall.setVertices(v0.position, v1.position, v3.position);
    
    // Recursion, for multiple refinements.
    if (numberOfTimes > 1)
    {
        numberOfTimes--;
        // The newly added vertex shares all three new triangles, so remember their indices.
        const unsigned idx0 = v3.triangleIndices[0].first;
        const unsigned idx1 = v3.triangleIndices[1].first;
        const unsigned idx2 = v3.triangleIndices[2].first;
        // Then refine each of the new triangles.
        refineTriangle(idx0, numberOfTimes);
        refineTriangle(idx1, numberOfTimes);
        refineTriangle(idx2, numberOfTimes);
    }
}

void TriangleMeshWall::setPeriodicCompanions(const std::vector<std::pair<unsigned, unsigned>>& periodicCompanions)
{
    periodicCompanions_ = periodicCompanions;
}

std::vector<unsigned> TriangleMeshWall::getPeriodicCompanions(const unsigned index)
{
    std::vector<unsigned> indices;
    
    // Loop through the companions, compare the index with both pair values and add the other to the return vector in case of a match.
    for (auto& p : periodicCompanions_)
    {
        if (p.first == index)
            indices.push_back(p.second);
        else if (p.second == index)
            indices.push_back(p.first);
    }
    
    return indices;
}

void TriangleMeshWall::moveVerticesToMatchVolume(std::vector<Vec3D> displacements, Mdouble targetVolume, int maxNumRecursiveCalls)
{
    // Calculating the change in volume due to a moving triangle isn't analytically possible (not that I know of).
    // This is because when two vertices are moved, they don't necessarily lie in a plane with the original positions.
    // An approximation would be to split up the plane by picking one of the two diagonals, but which one to pick?
    // Luckily we are in a mesh and therefore the choice doesn't matter, as the neighbour will correct for any errors.
    // Only the edges of the mesh will then be an approximation, which is just something we have to live with.
    
    // The volume of the shape formed by the original and moved triangle can be calculated by splitting it up into 3
    // tetrahedrons. Each of these tetrahedrons is automatically formed when moving and updating the vertices one by one.
    // Therefore we simply loop through all the vertices in the mesh and for each of the triangles connected to it, we
    // calculate the volume of the tetrahedron formed by the 3 original vertices and the current updated vertex. The
    // latter is only actually updated after these volume calculations, so that they're not influenced by it. However
    // when this vertex is referenced in the remainder of the loop, the updated version is used, since that then forms
    // the second or third tetrahedron.
    
    // The total volume is kept track of and together with the target volume a ratio which is used to correct the vertex
    // displacements. A recursive call is then made with the corrected displacements. The recursion is stopped when the
    // ratio equals 1 and only then we actually know the displacements needed to get a change in volume matching the
    // target volume, which is then used to actually update the real vertices, i.e. the privately stored ones.
    
    if (targetVolume == 0.0)
        return;

    if (displacements.size() != vertices_.size())
    {
        logger(WARN, "Displacements vector must be equal in size as the number of vertices in the mesh.");
        return;
    }
    
    // In order to know if we are on the first run of this recursive method,
    // we store a static boolean and reset it once the recursion ends.
    static bool firstRecursiveCall = true;
    if (firstRecursiveCall)
    {
        firstRecursiveCall = false;
        
        // The vertex periodic companions are handled by adding the displacements to each other (only done once, on first run).
        // As one vertex can have multiple periodic companions, we first copy the displacements vector.
        // To prevent unnecessarily copying, however, we first check if there are any periodic companions.
        if (!periodicCompanions_.empty())
        {
            std::vector<Vec3D> disps = displacements;
            for (auto& p : periodicCompanions_)
            {
                displacements[p.first] += disps[p.second];
                displacements[p.second] += disps[p.first];
            }
        }
    }
    
    // Every call a clean copy of the original real vertices is made, so that they can be updated without affecting anything.
    std::vector<Vertex> vertices = vertices_;
    Mdouble totalVolume = 0.0;
    
    // Loop through vertices to displace them and calculate the change in volume for each triangle connected to it.
    for (int i = 0; i < vertices.size(); i++)
    {
        // Loop through each triangle connected to this vertex.
        for (auto& p : vertices[i].triangleIndices)
        {
            Triangle& t = triangles_[p.first];
            // Calculate the volume of the tetrahedron formed by each of the vertices of the triangle
            // (one of which is the current vertex) and the updated position of the current vertex.
            totalVolume += getVolumeTetrahedron(vertices[t.vertexIndices[0]].position,
                                                vertices[t.vertexIndices[1]].position,
                                                vertices[t.vertexIndices[2]].position,
                                                vertices[i].position + displacements[i]);
        }
        
        // Only now actually update the current vertex position.
        vertices[i].position += displacements[i];
    }
    
    // Ratio between the target volume and actual volume
    Mdouble ratio = targetVolume / totalVolume;

    // Stop the recursive call when the ratio is equal to 1 (within certain tolerance).
    // Or when the ratio does not have a finite value (totalVolume turned out to be 0).
    Mdouble tol = 1.0e-6;
    if (mathsFunc::isEqual(ratio, 1.0, tol) || !std::isfinite(ratio) || maxNumRecursiveCalls == 0)
    {
        if (maxNumRecursiveCalls == 0)
            logger(WARN, "[TriangleMeshWall::moveVerticesToMatchVolume()] Maximum number of recursion calls reached for wall with id %.", getId());

        // Update the real vertices.
        int num = vertices_.size();
        for (int i = 0; i < num; i++)
            vertices_[i].position += displacements[i];
        
        // Update the real triangles wall.
        for (Triangle& t : triangles_)
            t.wall.moveVertices({ displacements[t.vertexIndices[0]], displacements[t.vertexIndices[1]], displacements[t.vertexIndices[2]] });
    
        updateBoundingBox();
        
        // Reset static boolean for the next time this function is first called.
        firstRecursiveCall = true;
        return;
    }
    
    // Multiply all displacements with the ratio.
    for (Vec3D& d : displacements)
        d *= ratio;
    // Recursive call.
    moveVerticesToMatchVolume(displacements, targetVolume, --maxNumRecursiveCalls);
}

Mdouble TriangleMeshWall::getVolumeTetrahedron(const Vec3D& a, const Vec3D& b, const Vec3D& c, const Vec3D& d)
{
    return std::fabs(Vec3D::dot((a-d), Vec3D::cross((b-d), (c-d)))) / 6.0;
}

void TriangleMeshWall::addToMesh(TriangleMeshWall mesh)
{
    // All triangles are simply added to the triangles_ vector. When adding vertices, any triangle index of the added
    // mesh simply becomes: updated index = index + number of triangles originally in this mesh.
    // Similar for the vertices, however vertices which share their position with a vertex already in the mesh are
    // ignored. So when adding triangles, the vertex indices can sometimes reference to a vertex already in the mesh.
    // The easiest way to implement is to temporarily store the updated vertex indices, and to later assign them.
    
    // Remember the number of vertices and triangles originally in this mesh.
    const int numVertices = vertices_.size();
    const int numTriangles = triangles_.size();
    // For each of the added vertices, their updated vertex is temporarily stored in a vector.
    // This is the cleanest/easiest way to later update the vertex indices of each added triangle.
    std::vector<unsigned> updatedVertexIndices;
    updatedVertexIndices.reserve(mesh.vertices_.size());
    
    // Loop through the vertices to be added.
    for (Vertex& v : mesh.vertices_)
    {
        // The triangle indices can simply be updated.
        for (auto& p : v.triangleIndices)
            p.first += numTriangles;
        
        // Initialize updated vertex index as -1, to later check if it has been set or not.
        int index = -1;
        // Loop through original vertices and compare vertex positions.
        for (int i = 0; i < numVertices; i++)
        {
            if (v.position.isEqualTo(vertices_[i].position, std::numeric_limits<Mdouble>::epsilon()))
            {
                // When vertices share position, update vertex index.
                index = i;
                break;
            }
        }
    
        // When the updated index is not set, the vertex is not in the mesh yet and should be added.
        if (index == -1)
        {
            index = vertices_.size();
            vertices_.push_back(v);
        }
        else
        {
            // The vertex is already in the mesh, add the triangle indices of the added mesh to it.
            for (auto& p : v.triangleIndices)
                vertices_[index].triangleIndices.push_back(p);
        }
        
        // Store the updated index.
        updatedVertexIndices.push_back(index);
    }

    // Loop through the triangles to be added.
    for (Triangle& t : mesh.triangles_)
    {
        // Update the vertex indices.
        for (unsigned& idx : t.vertexIndices)
            idx = updatedVertexIndices[idx];
        triangles_.push_back(t);
    }
    
    updateBoundingBox();
}

void TriangleMeshWall::createParallelogramMesh(const Vec3D& P0, const Vec3D& P1, const Vec3D& P2, Mdouble resolutionU,
                                                       Mdouble resolutionV, bool periodicInU, bool periodicInV)
{
    // Length of the sides of the parallelogram
    const Mdouble lengthU = (P2 - P0).getLength();
    const Mdouble lengthV = (P1 - P0).getLength();
    // Direction vectors of the parallelogram sides
    const Vec3D dirU = (P2 - P0) / lengthU;
    const Vec3D dirV = (P1 - P0) / lengthV;

    // Number of points is 1 more than number of segments
    const int numU = getNumberOfSegmentsAndResolution(lengthU, resolutionU) + 1;
    const int numV = getNumberOfSegmentsAndResolution(lengthV, resolutionV) + 1;

    std::vector<Vec3D> points;
    std::vector<std::array<unsigned, 3>> cells;

    for (int j = 0; j < numV; j++)
    {
        for (int i = 0; i < numU; i++)
        {
            const Vec3D p = P0 + dirU * i * resolutionU + dirV * j * resolutionV;
            points.push_back(p);
            addZigZagDiagonalCells(cells, numU, numV, i, j);
        }
    }

    set(points, cells);

    // Set possible periodic companions for first/last row/column.
    std::vector<std::pair<unsigned, unsigned>> periodic;
    if (periodicInU)
    {
        // First/last column
        for (int j = 0; j < numV; j++)
            periodic.emplace_back(j * numU, (j+1) * numU - 1);
    }
    if (periodicInV)
    {
        // First/last row
        for (int i = 0; i < numU; i++)
            periodic.emplace_back(i, i + (numV-1) * numU);
    }
    if (periodicInU && periodicInV)
    {
        // Periodic in diagonals u/v
        periodic.emplace_back(0, numU * numV - 1);
        periodic.emplace_back(numU - 1, numU * (numV-1));
    }
    setPeriodicCompanions(periodic);
}

int TriangleMeshWall::getNumberOfSegmentsAndResolution(const Mdouble length, Mdouble& resolution)
{
    // Round up to get a whole number of segments, then update the resolution (might remain unchanged).
    double num = std::ceil(length / resolution);
    resolution = length / num;
    return static_cast<int>(num);
}

void TriangleMeshWall::addZigZagDiagonalCells(std::vector<std::array<unsigned, 3>>& cells, int numU, int numV, int i, int j)
{
    // Not in last row/column
    if (i < numU - 1 && j < numV - 1)
    {
        // Bottom left
        const unsigned baseIdx = i + j * numU;

        // Check moduli to create zigzag diagonals (both being odd or even)
        if ((i % 2 == 0) != (j % 2 == 0))
        {
            // Diagonal from bottom left to top right
            // Bottom left, top left, top right (clockwise)
            cells.push_back({baseIdx, baseIdx + numU, baseIdx + numU + 1});
            // Bottom left, top right, bottom right (clockwise)
            cells.push_back({baseIdx, baseIdx + numU + 1, baseIdx + 1});
        }
        else
        {
            // Diagonal from top left to bottom right
            // Bottom left, top left, bottom right (clockwise)
            cells.push_back({baseIdx, baseIdx + numU, baseIdx + 1});
            // Top left, top right, bottom right (clockwise)
            cells.push_back({baseIdx + numU, baseIdx + numU + 1, baseIdx + 1});
        }
    }
}

void TriangleMeshWall::createFourPointMesh(const Vec3D& P0, const Vec3D& P1, const Vec3D& P2, const Vec3D& P3, int numSegmentsU, int numSegmentsV)
{
    // To get to a point in the mesh, all we need to do is to walk along the bottom edge and then walk up along the
    // linearly interpolated vector between the left and right edge. We step with a certain fraction of the length of
    // the vectors, dependent on the number of segments in that direction.

    // Bottom direction vector (u-direction)
    const Vec3D dirU = P3 - P0;
    // Left and right direction vector (v-direction)
    const Vec3D dirV0 = P1 - P0;
    const Vec3D dirV1 = P2 - P3;
    // Step size for linear interpolation between left and right direction vectors
    Vec3D deltaDirV = (dirV1 - dirV0) / numSegmentsU;

    // Number of points is 1 more than number of segments
    const int numU = numSegmentsU + 1;
    const int numV = numSegmentsV + 1;

    std::vector<Vec3D> points;
    std::vector<std::array<unsigned, 3>> cells;

    for (int j = 0; j < numV; j++)
    {
        for (int i = 0; i < numU; i++)
        {
            // Update direction vector in v-direction
            Vec3D dirV = dirV0 + i * deltaDirV;
            const Vec3D p = P0 + dirU * i / numSegmentsU + dirV * j / numSegmentsV;
            points.push_back(p);
            addZigZagDiagonalCells(cells, numU, numV, i, j);
        }
    }

    set(points, cells);
}

bool TriangleMeshWall::isLocal(Vec3D& min, Vec3D& max) const
{
    min = boundingBoxMinGlobal_;
    max = boundingBoxMaxGlobal_;
    return true;
}

bool TriangleMeshWall::isWithinBoundingBox(const Vec3D& position, Mdouble radius) const
{
    const Vec3D min = boundingBoxMin_ - Vec3D(radius, radius, radius);
    const Vec3D max = boundingBoxMax_ + Vec3D(radius, radius, radius);
    return position >= min && max >= position;
}

void TriangleMeshWall::updateBoundingBox()
{
    // Initialize min max.
    Vec3D min(constants::inf, constants::inf, constants::inf);
    Vec3D max(-constants::inf, -constants::inf, -constants::inf);
    
    // Loop through vertices and store min max positions.
    for (const Vertex& v : vertices_)
    {
        const Vec3D& pos = v.position;
        min = Vec3D::min(min, pos);
        max = Vec3D::max(max, pos);
    }
    
    boundingBoxMin_ = min;
    boundingBoxMax_ = max;
    
    updateBoundingBoxGlobal();
}

void TriangleMeshWall::updateBoundingBoxGlobal()
{
    // This basically is putting a bounding box around the rotated and translated lab frame bounding box.
    
    // All corners of the lab frame bounding box. Using min max for clarity.
    Vec3D min = boundingBoxMin_, max = boundingBoxMax_;
    Vec3D v0 = Vec3D(min.X, min.Y, min.Z);
    Vec3D v1 = Vec3D(min.X, min.Y, max.Z);
    Vec3D v2 = Vec3D(max.X, min.Y, max.Z);
    Vec3D v3 = Vec3D(max.X, min.Y, min.Z);
    Vec3D v4 = Vec3D(min.X, max.Y, min.Z);
    Vec3D v5 = Vec3D(min.X, max.Y, max.Z);
    Vec3D v6 = Vec3D(max.X, max.Y, max.Z);
    Vec3D v7 = Vec3D(max.X, max.Y, min.Z);
    
    // Store in array so we can do a auto loop by reference.
    std::array<Vec3D, 8> corners = { v0, v1, v2, v3, v4, v5, v6, v7 };
    
    // Reinitialize min max to use in loop.
    min = Vec3D(constants::inf, constants::inf, constants::inf);
    max = Vec3D(-constants::inf, -constants::inf, -constants::inf);
    
    // Loop through corners, rotate and translate them to global frame, get global min max position.
    for (Vec3D& pos : corners)
    {
        getOrientation().rotate(pos);
        pos += getPosition();
        min = Vec3D::min(min, pos);
        max = Vec3D::max(max, pos);
    }
    
    boundingBoxMinGlobal_ = min;
    boundingBoxMaxGlobal_ = max;
}

void TriangleMeshWall::setPosition(const Vec3D& position)
{
    BaseWall::setPosition(position);
    updateBoundingBoxGlobal();
}

void TriangleMeshWall::setOrientation(const Quaternion& orientation)
{
    BaseWall::setOrientation(orientation);
    updateBoundingBoxGlobal();
}

void TriangleMeshWall::move(const Vec3D& move)
{
    BaseWall::move(move);
    updateBoundingBoxGlobal();
}

void TriangleMeshWall::rotate(const Vec3D& angularVelocity)
{
    BaseWall::rotate(angularVelocity);
    updateBoundingBoxGlobal();
}