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

#ifndef MEMBRANE_H
#define MEMBRANE_H

#include <Mercury3D.h>
#include <vector>
#include <Math/Vector.h>
#include <array>

#include <Particles/BaseParticle.h>
#include <Walls/MeshTriangle.h>


/*!
 * \class Membrane
 * \brief A Membrane consists of masses connected by springs.
 * \details This class assumes, that the masses are connected in a triangular mesh.
 * Additionally, for correct calculation of the spring constant for the distance
 * spring, a hexagonal unit cell is assumed.
 *
 * This class creates the particles required for the mass spring system. 
 * For contact detections, further triangular wall elements are created
 *
 * Please make sure to provide appropriate a species for the particles as well as
 * the wall elements (e.g giving the particles a species so that they do not interact
 * with other particles in the simulation).
 *
 * An instance of the membrane may be initialized with the following lines, e.g.
 * by using a STL file:
 * \code
 * // Initialize the membrane and set neccesary quantities
 * Membrane membrane = Membrane(membraneSpecies_, membraneParticleSpecies_);
 * membrane.setDPMBase(this);
 * membrane.setKnAndCrittDampCoeff(Kn_, critDampCoeff_);
 * membrane.setKeAndKd(Ke_, Kd_);
 *
 * // Set the (approximate) radius of the vertex particles
 * membrane.setParticleRadius(membraneParticleRadius_);
 *
 * // Create an initial particle for the vertex particles
 * SphericalParticle p0;
 * p0.setSpecies(membraneParticleSpecies_);
 * p0.setRadius(membrane.getParticleRadius());
 * p0.setVelocity(Vec3D(0.0, 0.0, 0.0));
 *
 * // Create the system using the specified particle and stl file
 * membrane_.readFromSTL(p0, "geometry.stl");
 * \endcode
 *
 * Note that the method computeAdditionalForces needs to be called from the overridden
 * DPMBase::computeAdditionalForces every timestep for the Membrane to work.
 *
 * As an example please also have a look at the Driver code MembraneDemo
 */
class Membrane : public BaseObject
{
public:
    
    /*!
     * A struct containing the information needed for calculations along one edge
     */
    struct Edge
    {
        // 3 and 4
        /*!
         * Stores the ids of the particles of the edge
         */
        std::array<unsigned int, 2> vertexId;
        
        /*!
         * Stores a pointer to the particles of the edge
         */
        std::array<BaseParticle*, 2> vertex = {{nullptr}};

        // 1 and 2
        /*!
         * Stores the ids of the remaining particles of faces connected by the edge
         */
        std::array<unsigned int, 2> faceVertexId;
        
        /*!
         * Stores a pointer to the remaining particles of faces connected by the edge
         */
        std::array<BaseParticle*, 2> faceVertex = {{nullptr}};
        
        /*!
         * Stores a pointer to the faces connected by the edge
         */
        std::array<MeshTriangle*, 2> face = {{nullptr}};

        /*!
         * Stores the initial area of the connected triangles
         */
        std::array<Mdouble, 2> faceInitialArea;
        
        /*!
         * Stores values calculated for bending
         */
        std::array<Mdouble, 8> uPre;

        /*!
         * Stores the initial difference between the edge particles
         */
        Vec3D initialState;
        
        /*!
         * Stores initial length of the edge
         */
        Mdouble initialLength;
        
        /*!
         * Stores the initial sine of half the bending angle
         */
        Mdouble initialSineHalfTheta;
        
        /*!
         * Stores the effective mass of the edge
         */
        Mdouble effectiveMass;
        
        /*!
         * Stores if the edge should calculate the disntance spring forces
         */
        bool isStretchActive;
        
        /*!
         * Stores if the edge should calculate the bending penalty forces
         */
        bool isBendActive;

        /*!
         * \brief Calculates and applies all neccesary forces
         */
        void applyForce(Mdouble Kn, Mdouble critDampCoeff, Mdouble Ke, Mdouble Kd, bool bendingAreaConstant);
        
        /*!
         * \brief Apply the force due to stretching only.
         */
        void applyStretchForce(Mdouble Kn, Mdouble critDampCoeff);
        
        /*!
         * \brief Apply a force due to bending
         */
        void applyBendForce(Mdouble Kn, Mdouble Kd, bool bendingAreaConstant);
        
        /*!
         * \brief Calculate some prefactors for the bending penalty
         */
        void calculateUPre(Vec3D& state, Mdouble& length, std::array<Mdouble, 2>& faceArea);
        
        /*!
         * \brief Calculate the sine of half the bending angle
         */
        Mdouble getSineHalfTheta();
        
        /*!
         * \brief handle the partical removal
         */
        void handleParticleRemoval(unsigned int id);
        
        /*!
         * \brief Handle the particle addition
         */
        void handleParticleAddition(unsigned int id, BaseParticle* p);
        
        /*!
         * \brief check if the edge should calculate bending or stretch forces
         */
        void checkActive();
    };

    /*!
     * \brief Default constructor.
     */
    Membrane();

    /*
     * \brief Constructor to create a membrane with two species
    */
    Membrane(ParticleSpecies* membraneSpecies, ParticleSpecies* membraneParticleSpecies);

    /*!
     * \brief Copy constructor.
     */
    // Membrane(const Membrane& other);


    /*!
     * \brief Destructor.
     */
    ~Membrane() override;

    /*!
     * Copy assignment operator.
     */
    // Membrane& operator=(const Membrane& other);

    /*!
     * \brief Copy method. It calls the copy constructor of this Object
     */
    Membrane* copy() const;

    /*!
     * \brief Reads a Membrane from an input stream, for example a restart file.
     */
    void read(std::istream& is) override;

    /*!
     * \brief Writes a Membrane to an output stream, for example a restart file.
     */
    void write(std::ostream& os) const override;

    /*!
     * \brief Returns the name of the object, here the string "Membrane".
     */
    std::string getName() const override;

    /*!
     * \brief Set the parameters needed for the stretching forces
     */
    void setKnAndCrittDampCoeff(Mdouble Kn, Mdouble critDampCoeff);
    
    /*!
     * \brief Set the elastic modulus and thickness of the membrane
     * \deltails The supplied values are used to calculate the spring constant
     * of the membrane.
     */
    void setElasticModulusAndThickness(Mdouble E, Mdouble thickness);
    
    /*!
     * \brief Set the spring constant of the membrane
     */
    void setSpringConstant(Mdouble k);
    
    /*!
     * \brief Set damping coefficient for the distance springs
     */
    void setCriticalDampingCoefficient(Mdouble coeff);
    
    /*!
    * \brief Set the parameters needed for the bending forces
    */
    void setKeAndKd(Mdouble Ke, Mdouble Kd);
    
    /*!
    * \details If set to true, the bending penalty is calulated using the initial
    * triangle area instead of the current one.
    */
    void setBendingAreaConstant(bool areaConstant);
    
    /*!
    * \brief Set the thickness of the membrane
    */
    void setThickness(Mdouble thickness);
    
    
    Mdouble getKn() { return Kn_; }
    Mdouble getCriticalDampingCoefficient() { return critDampCoeff_; }
    Mdouble getKe() { return Ke_; }
    Mdouble getKd() { return Kd_; }
    Mdouble getThickness() { return thickness_; }
    Mdouble getBendingAreaConstant() { return bendingAreaConstant_; }
    
    /*!
     * \brief Set the radius of the vertex particles
     */
    void setParticleRadius(Mdouble radius);
    
    /*!
     * \brief Returns the radius of the vertex particles
     */
    Mdouble getParticleRadius() { return particleRadius_; }
    
    /*!
     * \brief Set a pointer to DPMBase
     */
    void setDPMBase(DPMBase* dpm) { DPMBase_ = dpm; }
    
    /*!
     * \brief Get the stored pointer to DPMBase
     */
    DPMBase* getDPMBase() { return DPMBase_; }
    
    /*!
     * \brief Returns a vector with pointers to the vertex particles
     */
    std::vector<BaseParticle*> getVertexParticles(){ return vertexParticle_; };
    
    /*!
     * \brief Returns a vecter with pointers to the mesh triangles
     */
    std::vector<MeshTriangle*> getFaces(){ return face_; };    
    
    /*!
     * \brief save the current positions of the vertex particles to a stream
     */
    void saveVertexPositions(std::ostream& os);
    
    /*!
     * \brief Load the positions of the vertex particles from a stream and apply 
     * them to existing particles.
     */
    void loadVertexPositions(std::istream& is);
    
    /*!
     * \brief Save the Membrane as a .off file.
     */
    void saveAsOFF(unsigned int d);
    
    /*!
     * \brief Save the Membrane as a .stl file
     */
    void saveAsSTL(std::string fileName);
    
    /*!
     * \brief Load the Membrane geometry from a .stl file
     */
    void loadFromSTL(BaseParticle& p0, std::string fileName);
    
    /*!
     * \brief Load the Membrane geometry from a .stl file
     */
    void loadFromSTL(BaseParticle& p0, std::string fileName, Mdouble eps);
    
    /*!
     * \brief Build the geometry from specified positions and their connectivity
     */
    void buildMesh(BaseParticle& p0, std::vector<Vec3D> vertexPositions, std::vector<unsigned int> edgeVertices, std::vector<unsigned int> faceVertices);
    
    /*!
     * \brief Calculates an updated radius for every vertex particle in order to 
     * account for the correct mass.
     */
    void adjustVertexParticleSize();
    
    /*!
     * \brief Set the correct edge mass by taking the mass from the conencted vertices.
     */
    void updateEdgeMass();
    
    /*!
     * \brief Compute the forces due to the mass spring system
     */
    void computeAdditionalForces();
    
    /*!
     * \brief Apply a surface pressure to the membrane
     */
    void applyPressure(Mdouble pressure);
    
    /*!
     * \brief Handles the removal of vertex particles from the particle handler
     */
    void handleParticleRemoval(unsigned int id);
    
    /*!
     * \brief Handles the addition of vertex particles to the particle handler
     */
    void handleParticleAddition(unsigned int id, BaseParticle* p);

    /*!
     * \brief Calculate the volume of the membrane
     */
    Mdouble getVolume();

private:
    
    /*!
     * \brief Handles the actual creation of vertex particles
     */
    void createVertexParticles(BaseParticle& p0, std::vector<Vec3D> vertexPositions);
    
    /*!
     * \brief Compute the forces due to the mass spring system
     */
    void initializeEdgeBendingQuantities();
    
    /*!
     * \brief Update the faces to have the correct neighbors
     */
    void updateFaceNeighbors();
    
    /*!
     * \brief Helper function to check if a given vertex already exists
     */
    unsigned int addVertex(std::vector<Vec3D>& vertices, Vec3D pos, Mdouble eps);
    
    /*!
     * \brief Helper function to check if a given edge already exists
     */
    bool addEdge(std::vector<unsigned int>& edges, std::map<std::pair<unsigned int, unsigned int>, int>& map, unsigned int v0, unsigned int v1);
    
    /*!
     * Stores the pointers to the vertex particles
     */
    std::vector<BaseParticle*> vertexParticle_;
    
    /*!
     * Stores the ids of the vertex particles 
     */
    std::vector<unsigned int> vertexParticleId_;
    
    /*!
     * Sores the id of the first particle created for the membrane.
     * This value helps for the eventual MPI parallelization.
     */
    unsigned int vertexInitId_;

    /*!
     * Stores pointers to the MeshTriangle objects
     */
    std::vector<MeshTriangle*> face_;
    
    /*!
     * Stores the 3 indices per face which can be used to get the ids from vertexParticleId_
     * // Todo: remove this
     */
    std::vector<unsigned int> faceVertices_;

    /*!
     * Stores the edges
     */
    std::vector<Edge> edge_;

    /*!
     * Stores the radius of membrane particles
     */
    Mdouble particleRadius_;

    /*!
     * Stores the spring constant for stretching
     */
    Mdouble Kn_;
    
    /*!
     * Stores damping coefficient for the normal springs
     */
    Mdouble critDampCoeff_;
    
    /*!
     * Stores spring constant for the bending
     */
    Mdouble Ke_;
    
    /*!
     * Stores dissipation constant for the bending
     */
    Mdouble Kd_;
    
    /*!
     * Stores thickness of the Membrane
     */
    Mdouble thickness_;
    
    /*!
     * Stores a pointer to DPMBase
     */
    bool bendingAreaConstant_;
    
    /*!
     * Stores a pointer to the membrane species
     */
    ParticleSpecies* membraneSpecies_ = nullptr;
    
    /*!
     * Stores a pointer to the membrane particle species
     */
    ParticleSpecies* membraneParticleSpecies_ = nullptr;
    
    /*!
     * Stores a pointer to DPMBase
     */
    DPMBase* DPMBase_ = nullptr;
};


#endif
