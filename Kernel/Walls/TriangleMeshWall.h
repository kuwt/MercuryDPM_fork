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

#ifndef MERCURYDPM_TRIANGLEMESHWALL_H
#define MERCURYDPM_TRIANGLEMESHWALL_H

#include "BaseWall.h"
#include "TriangleWall.h"

class TriangleMeshWall : public BaseWall
{
public:
    
    /*!
     * \brief Default constructor.
     */
    TriangleMeshWall() = default;
    
     /*!
      * \brief Copy constructor.
      * @param other The TriangleMeshWall that must be copied.
      */
    TriangleMeshWall(const TriangleMeshWall& other);
    
    /*!
     * \brief Constructor creating triangles
     * @param points Vector of vertices
     * @param cells Vector of cells, each cell holds three indices to the points vector
     * @param species
     */
    TriangleMeshWall(const std::vector<Vec3D>& points, const std::vector<std::array<unsigned, 3>>& cells, const ParticleSpecies* species = nullptr);

    /*!
     * \brief Constructor creating parallelogram shaped mesh with a given resolution in u- and v-direction.
     * @param P0, P1, P2 The three points, in clockwise order, forming the parallelogram (fourth point is implied).
     * @param resolutionU Maximum length of a segment in u-direction.
     * @param resolutionV Maximum length of a segment in v-direction.
     * @param species
     * @param periodicInU Whether or not the first and last column of vertices lie in a periodic boundary and should be set as periodic companions.
     * @param periodicInV Whether or not the first and last row of vertices lie in a periodic boundary and should be set as periodic companions.
     */
    TriangleMeshWall(const Vec3D& P0, const Vec3D& P1, const Vec3D& P, Mdouble resolutionU, Mdouble resolutionV,
                     const ParticleSpecies* species = nullptr, bool periodicInU = false, bool periodicInV = false);
    
    /*!
     * \brief Destructor
     */
    ~TriangleMeshWall() override = default;
    
    /*!
     * \brief Copy assignment operator.
     * @param other The TriangleMeshWall that must be copied.
     * @return
     */
    TriangleMeshWall& operator=(const TriangleMeshWall& other);
    
    /*!
     * \brief Wall copy method.
     * @return Pointer to the copy.
     */
    TriangleMeshWall* copy() const override;
    
    /*!
     * \brief Reads a TriangleMeshWall from an input stream, for example a restart file.
     */
    void read(std::istream& is) override;
    
    /*!
     * \brief Writes a TriangleMeshWall to an output stream, for example a restart file.
     */
    void write(std::ostream& os) const override;
    
    /*!
     * \brief Returns the name of the object.
     * @return "TriangleMeshWall"
     */
    std::string getName() const override;
    
    /*!
     * \brief Set function creating triangles
     * @param points Vector of vertices
     * @param cells Vector of cells, each cell holds three indices to the points vector
     */
    void set(const std::vector<Vec3D>& points, const std::vector<std::array<unsigned, 3>>& cells);
    
    /*!
     * \brief Checks if a collision is happening with any of the TriangleWalls. NOTE: this does not handle actual
     * interactions, but can only be used to know IF one (of possible many) interaction occurs.
     * @param p
     * @param distance
     * @param normal_return
     * @return Whether or not there is a collision.
     */
    bool getDistanceAndNormal(const BaseParticle& p, Mdouble& distance, Vec3D& normal_return) const override;
    
     /*!
      * \brief Gets the interaction between a given BaseParticle and all TriangleWalls in this mesh at a given time.
      * In case of multiple contacts, chooses the one with greatest overlap.
      * @param p
      * @param timeStamp
      * @param interactionHandler
      * @return
      * \note Particle positions are rotated to lab frame, however the particle orientation itself is not considered,
      * so most likely doesn't work properly for non-spherical particles.
      */
    BaseInteraction* getInteractionWith(BaseParticle* p, unsigned timeStamp, InteractionHandler* interactionHandler) override;
    
    void resetForceTorque(int numberOfOMPthreads) override;
    
    void writeVTK(VTKContainer& vtk) const override;
    
    void setSpecies(const ParticleSpecies* species);
    
    /*!
     * \brief Updates a vertex by a given change in position.
     * @param index Index of vertex to move.
     * @param dP Change in position.
     */
    void moveVertex(unsigned index, const Vec3D& dP);
    
    /*!
     * \brief Removes a single triangle from the mesh.
     * @param index Index of triangle to remove.
     * @param removeFreeVertex Removes a vertex when it has no more triangles attached to it.
     */
    void removeTriangle(unsigned index, bool removeFreeVertex = true);
    
    /*!
     * \warning For now can only be used in initial setup, i.e. when there are no interactions yet. Interactions have to
     * be projected onto newly formed triangles, which is not yet implemented.
     * \brief Makes a really basic refinement by adding a vertex in the middle of the triangle to create 3 new ones.
     * Note: Multiple refinements will definitely NOT create a nice mesh.
     * @param index Index of triangle to refine.
     * @param numberOfTimes How many times to refine.
     */
    void refineTriangle(unsigned index, unsigned numberOfTimes = 1);
    
    /*!
     * \brief Sets the vector with pairs of indices of vertices which lie on a periodic boundary and are each others
     * 'companion'. E.g. when moving one vertex, the other is also moved the same amount. A vertex can have multiple
     * companions (when lying on multiple periodic boundaries) and those should be in the vector as multiple pairs.
     * @param periodicCompanions Vector of pairs containing the companions.
     */
    void setPeriodicCompanions(const std::vector<std::pair<unsigned, unsigned>>& periodicCompanions);
    
    /*!
     * \brief Moves all vertices by given displacements and corrects to match the change in volume with the target volume.
     * The final displacements will be proportional to the initial ones, but their actual values will have change.
     * \warning This calculates absolute volume and might cause unexpected results when volume is both 'added' and 'removed'.
     * @param displacements Initial change in position for each vertex.
     * @param targetVolume Volume that the change in volume must match.
     * @param maxNumRecursiveCalls Maximum number of recursive calls allowed, in case certain tolerances are set too
     * tight and infinite recursion occurs.
     */
    void moveVerticesToMatchVolume(std::vector<Vec3D> displacements, Mdouble targetVolume, int maxNumRecursiveCalls = 15);
    
    /*!
     * \brief Calculates the volume of a tetrahedron formed by 4 vertices.
     * @return Volume of tetrahedron
     */
    Mdouble getVolumeTetrahedron(const Vec3D& a, const Vec3D& b, const Vec3D& c, const Vec3D& d);
    
    /*!
     * \brief Adds a mesh to this mesh. When vertices share a position, only the one in the current mesh is kept and
     * both meshes are then "connected" at that point.
     * @param mesh Mesh to be added.
     */
    void addToMesh(TriangleMeshWall mesh);

    /*!
     * \brief Creates a parallelogram shaped mesh with a given resolution in u- and v-direction.
     * @param P0, P1, P2 The three points, in clockwise order, forming the parallelogram (fourth point is implied).
     * @param resolutionU Maximum length of a segment in u-direction.
     * @param resolutionV Maximum length of a segment in v-direction.
     * @param periodicInU Whether or not the first and last column of vertices lie in a periodic boundary and should be set as periodic companions.
     * @param periodicInV Whether or not the first and last row of vertices lie in a periodic boundary and should be set as periodic companions.
     */
    void createParallelogramMesh(const Vec3D& P0, const Vec3D& P1, const Vec3D& P2, Mdouble resolutionU, Mdouble resolutionV, bool periodicInU = false, bool periodicInV = false);

    /*!
     * \brief Creates mesh consisting of four points, with linearly interpolated segments. The points therefore don't
     * necessarily have to lie in plane.
     * @param P0, P1, P2, P3 The four points, in clock-wise order.
     * @param numSegmentsU Number of segments in u-direction.
     * @param numSegmentsV Number of segments in v-direction.
     */
    void createFourPointMesh(const Vec3D& P0, const Vec3D& P1, const Vec3D& P2, const Vec3D& P3, int numSegmentsU, int numSegmentsV);
    
    bool isLocal(Vec3D& min, Vec3D& max) const override;
    
    void setPosition(const Vec3D& position) override;
    void setOrientation(const Quaternion& orientation) override;
    void move(const Vec3D& move) override;
    void rotate(const Vec3D& angularVelocity) override;
   
protected:
    struct Triangle
    {
        /*!
         * \brief Default constructor.
         */
        Triangle() = default;
        
        /*!
         * \brief Constructor setting values.
         * @param w The TriangleWall.
         * @param vi Array of indices of the three vertices the triangle is made of.
         */
        Triangle(const TriangleWall& w, std::array<unsigned, 3> vi) : wall(w), vertexIndices(vi) {};
        
        /// The TriangleWall.
        TriangleWall wall;
        
        /// Array of indices to the vertices vector of the 3 vertices the triangle is made of.
        std::array<unsigned, 3> vertexIndices{};
    };
    
    struct Vertex
    {
        /*!
         * \brief Default constructor.
         */
        Vertex() = default;
        
        explicit Vertex(Vec3D P) : position(P) {};
        
        /// The vertex position
        Vec3D position;
        
        /// Vector of pairs, with for each pair first an index to the triangle in the triangles vector and
        /// second a local vertex index to know which of the three triangle vertices corresponds to this one.
        std::vector<std::pair<unsigned, unsigned>> triangleIndices;
    };

protected:
    std::vector<Triangle> triangles_;
    std::vector<Vertex> vertices_;
    
private:
    std::vector<std::pair<unsigned, unsigned>> periodicCompanions_;
    
     /*!
      * \brief Looks up the companions belonging to a vertex index, if it has any.
      * @param index Index of the vertex to lookup the companions for.
      * @return Vector of vertex indices of the companions.
      */
    std::vector<unsigned> getPeriodicCompanions(unsigned index);
    
    /// Bounding box corners in lab frame and global frame. Used to speed up contact detection. HGrid uses the global frame min max.
    Vec3D boundingBoxMin_, boundingBoxMax_;
    Vec3D boundingBoxMinGlobal_, boundingBoxMaxGlobal_;
    
    /*!
     * \brief Checks whether or not a part of the particle is within the bounding box
     * @param position Position of the particle
     * @param radius Interaction radius of the particle
     * @return True when within the box
     */
    bool isWithinBoundingBox(const Vec3D& position, Mdouble radius) const;
    
    /*!
     * \brief Sets the local and global bounding boxes. This method should be called any time a vertex is updated in lab frame.
     */
    void updateBoundingBox();
    
    /*!
     * \brief Sets the global bounding box, using the bounding box in lab frame. This method should be called when there
     * are changes in global frame, but not in lab frame. E.g. when calling setPosition().
     */
    void updateBoundingBoxGlobal();

    /*!
     * \brief Finds the number of segments and the resolution needed to perfectly fit a length.
     * @param length Total length
     * @param resolution Maximum length of a segment, will be made smaller when needed.
     * @return Resulting number of segments.
     */
    int getNumberOfSegmentsAndResolution(Mdouble length, Mdouble& resolution);

    /*!
     * \brief Helper function that adds two cells (triangles) to the provided cells vector. Used by createParallelogramMesh and createFourPointMesh.
     * @param cells The vector to which the cells are added.
     * @param numU The number of segments in u-direction.
     * @param numV The number of segments in v-direction.
     * @param i The cell index in u-direction.
     * @param j The cell index in v-direction.
     */
    void addZigZagDiagonalCells(std::vector<std::array<unsigned, 3>>& cells, int numU, int numV, int i, int j);
};


#endif //MERCURYDPM_TRIANGLEMESHWALL_H
