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

#ifndef MeshTriangle_H
#define MeshTriangle_H

#include <vector>
#include <Walls/BaseWall.h>
#include <Walls/InfiniteWall.h>
#include "Math/Vector.h"

/*!
 * \brief MeshTriangle implements a triangle whose vertex positions are defined 
 * by three particles.
 * \details It tracks the movement of the specified particles and updates its own 
 * potition every timestep. If neighboring objects along the edges and vertices
 * are known, their contacts will be taken into account to consider if particle
 * contacts should produce a force. Calculated contact forces will be transferred
 * to the certex particles.
 *
 * A triangle may e.g. be constructed with the following code.
 * \code
 * p0 = Vec3D(0, 0.1, 1);
 * p1 = Vec3D(0, 0.1, 0.9);
 * p2 = Vec3D(0, 0  , 1);
 *
 * SphericalParticle p;
 * p.setSpecies(someSpecies);
 * p.setRadius(2e-2);
 * p.setVelocity(Vec3D(0.0, 0.0, 0.0));
 * p.setHandler(&particleHandler);
 *
 * p.setPosition(p0);
 * unsigned int Id1 = particleHandler.copyAndAddObject(p)->getId();
 * p.setPosition(p1);
 * unsigned int Id2 = particleHandler.copyAndAddObject(p)->getId();
 * p.setPosition(p2);
 * unsigned int Id3 = particleHandler.copyAndAddObject(p)->getId();
 * 
 * MeshTriangle f;
 * f.setVertices(p0, p1, p2);
 * 
 * // Create all faces with their initial positions
 * f.setVertexIds(0, 1, 2);
 * f.setSpecies(someSpecies);
 * f.setMass(mass);
 * // Retrieve the correct positions from the particles
 * f.actionsAfterParticleGhostUpdate();
 * \endcode
 * For a longer example, please have a look at the class Membrane.
 */
class MeshTriangle : public BaseWall
{
public:

    // TODO write proper constructors (also copy constructors)
    /*!
     * \brief Default constructor.
     */
    MeshTriangle() = default;

    /*!
     * \brief Copy constructor.
     */
    MeshTriangle(const MeshTriangle& other) = default;

    /*!
     * \brief Destructor.
     */
    ~MeshTriangle() override = default;

    /*!
     * \brief Wall copy method. It calls the copy constructor of this Wall, useful for polymorphism
     */
    MeshTriangle* copy() const override
    { return new MeshTriangle(*this); }

    /*!
     * \brief Returns the name of the object, here the string "MeshTriangle".
     */
    std::string getName() const override
    { return "MeshTriangle"; }

    /*!
     * \brief Reads an MeshTriangle from an input stream, for example a restart file.
     */
    void read(std::istream& is) override;

    /*!
     * \brief Writes an MeshTriangle to an output stream, for example a restart file.
     */
    void write(std::ostream& os) const override;

    /*!
     * \brief Sets member variables such that the wall represents a triangle with vertices A, B, C
     *  - position_ is set to the center of mass of the wall
     *  - updateVertexAndNormal is called to set the remaining variables
     */
    void setVertices(Vec3D A, Vec3D B, Vec3D C);
    
    /*!
    * \brief Same as #setVertices(A,B,C), but sets the position explicitly.
    * The position is important when you rotate the wall, as the wall will be rotated around this position.
    */
    void setVertices(Vec3D A, Vec3D B, Vec3D C, Vec3D position);
    
    /*!
     * \brief Sets the velocity of the vertex points
     */
    void setVertexVelocities(Vec3D A, Vec3D B, Vec3D C);
    
    /*!
     * \brief Returns an array of the vertex coordinates
     */
    std::array<Vec3D,3> getVertices() const {return vertex_;}

    /*!
     * \brief Returns the face normal
     */
    Vec3D getFaceNormal() const {return faceNormal_;}
    
    /*!
     * \brief Returns the area of the triangle
     */
    Mdouble getArea() const {return area_;}
    
    
    void move(const Vec3D& move) override;

    /*!
     * \brief sets the ids of the vertex particles. Calls retrieveVertexParticles.
     */
    void setVertexIds(unsigned int i, unsigned int j, unsigned int k);
    
    /*!
     * \brief Returns an array containing the ids of the vertex particles
     */
    std::array<unsigned int, 3> getVertexIds() const { return vertexIds_; }
    
    // 
    void writeVTK(VTKContainer& vtk) const override;

    
    
    BaseInteraction* getInteractionWith(BaseParticle* p, unsigned timeStamp, InteractionHandler* interactionHandler) override;
    
    /*!
     * \brief Checks, if the forces of all interctions should be applied
     */
    void checkInteractions(InteractionHandler* interactionHandler, unsigned int timeStamp) override;
    
    /*!
     * \brief Compute the distance from the wall for a given BaseParticle and return if there is a collision.
     * If there is a collision, also return the normal vector.
     */
    bool getDistanceAndNormal(const BaseParticle& p, Mdouble& distance, Vec3D& normal_return) const override;
    bool getDistanceNormalOverlapType(const BaseParticle& p, Mdouble& distance, Vec3D& normal, Mdouble& overlap, unsigned int& type) const;
    
    /*!
     * \brief Calculates the local velocity at a specified point
     */
    const Vec3D getVelocityAtContact(const Vec3D& contact) const override;
    
    /*!
     * \brief Calculates the barycentric weight of a specified point
     */
    const Vec3D getBaricentricWeight(const Vec3D& contact) const;

    void rotate(const Vec3D& angularVelocity) override;

    /*!
     * \brief Determines if the triangle is considered local
     */
    bool isLocal(Vec3D& min, Vec3D& max) const override;

    /*!
     * \brief Determines if a given point is within the triangle
     */
    bool isInsideTriangle(const Vec3D &point) const;

    Mdouble getInvMass() const override;
    void setMass(Mdouble mass);

    /*!
     * \brief Actions executed on restart.
     */
    void actionsOnRestart() override;
    
    /*!
     * \brief actionsPerformed after the position update of (ghost-) particles.
     */
    void actionsAfterParticleGhostUpdate() override;
    
    /*!
     * \brief Handles the removal of particles to the particle Handler
     */
    void handleParticleRemoval(unsigned int id) override;
    
    /*!
     * \brief Handles the addition of particles to the particle Handler
     */
    void handleParticleAddition(unsigned int id, BaseParticle* p) override;
    
    /*!
     * \brief Tries to get pointers to all vertex particles from the handler
     */
    void retrieveVertexParticles();
    
    /*!
     * \brief Check if the triangle is considered active
     */
    void checkActive();
    
    bool getActive() { return isActive; }
    
    /*!
     * \brief Set the handler 
     */
    void setHandler(WallHandler* handler) override;
    
    /*!
     * \brief Apply a force pointing in normal direction corresponding to the specified pressure.
     */
    void applyPressure(Mdouble presure);
    
    /*!
     * \brief Apply the given force to the triangle
     */
    void applyForce(Vec3D force);
    
    
    /*!
     * stores references to the neighbors along the edges.
     */
    std::array<MeshTriangle*, 3> neighbor = {{nullptr}};
    
    /*!
     * stores references to the neighbors on corners.
     */
    std::vector<std::vector<unsigned int>> vertexNeighbors;

private:
    
    /*!
     * \brief Update vertexMin_, vertexMax_ and faceNormal_ for an updated position.
     */
    void updateVertexAndNormal();
    
    /*!
     * \brief Retrieve new positions from updated vertex particles
     */
    void updateVerticesFromParticles();
    
    /*!
     * stores the position of the vertices relative to the position of the wall but not rotated into the lab frame;
     * thus, if the wall rotates, these vertices have to be rotated as well
     */
    std::array<Vec3D, 3> vertex_;
    std::array<Vec3D, 3> vertexVelocity_;
    std::array<BaseParticle*, 3> vertexParticle_ = {{nullptr}};

    /*!
     * stores the min and max coordinate values of the vertices (needed for hGrid)
     */
    Vec3D vertexMin_;
    Vec3D vertexMax_;

    /*!
     * stores the wall normal n in n.x=p
     */
    std::array<Vec3D, 3> edgeNormal_;
    std::array<Vec3D, 3> edge_;
    // Used to ignore certain contacts
    std::array<unsigned int, 3> vertexIds_;
    std::array<double, 3> edgeLength_;


    /*!
     * stores the face normal, not rotated into the lab frame; thus, if the wall rotates, this normal has to be rotated as well
     */
    Vec3D faceNormal_;
    Mdouble area_;

    // TODO: replace defualt value by default value in proper constructor
    Mdouble invMass_= 0.0;
    
    bool isActive = 0;
};

#endif
