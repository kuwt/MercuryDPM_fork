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

#include "Membrane.h"
#include <Mercury3D.h>
#include "InteractionHandler.h"
#include "Particles/BaseParticle.h"

#include "Math/ExtendedMath.h"
#include "Math/Vector.h"

#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <climits>

Membrane::Membrane()
{
    logger(DEBUG, "Membrane() constructed.");
    vertexInitId_ = 0;
    particleRadius_ = 0;
    Kn_ = 0;
    critDampCoeff_ = 0;
    Ke_ = 0;
    Kd_ = 0;
    thickness_ = 0;
    bendingAreaConstant_ = 0;
}

Membrane::Membrane(ParticleSpecies* membraneSpecies, ParticleSpecies* membraneParticleSpecies) : Membrane()
{
    membraneSpecies_ = membraneSpecies;
    membraneParticleSpecies_ = membraneParticleSpecies;
}

// /*!
//  * \param[in] other The Membrane that must be copied.
//  */
// Membrane::Membrane(const Membrane& other)
//         : BaseObject(other)
// {
//     face_ = other.face_;
//     membraneSpecies_ = other.membraneSpecies_;
//     DPMBase_ = other.DPMBase_;
// 
//     logger(DEBUG, "Membrane(Membrane&) constructed.");
// }


Membrane::~Membrane()
{
    logger(DEBUG, "~Membrane() has been called.");
}

/*!
 * \return pointer to a Membrane object allocated using new.
 */
Membrane* Membrane::copy() const
{
    return new Membrane(*this);
}


/*!
 * \param[in] is The input stream from which the Membrane is read, usually a restart file.
 */
void Membrane::read(std::istream& is)
{
    
    logger.assert_debug(DPMBase_, "Error in Membrane::read: DPMBase is not knwon.");
    
    std::string dummy;
    unsigned int i, n, Id;

    // Reading in the Ids of all the vertex particles
    is >> dummy >> n;
    vertexParticle_.reserve(n);
    for (i=0; i<n; i++)
    {
        is >> Id;
        vertexParticleId_.push_back(Id);
        // this might retrieve a nullptr
        vertexParticle_.push_back(getDPMBase()->particleHandler.getObjectById(Id));
    }
    is >> dummy >> vertexInitId_;

    // Reading in the Ids of the faces
    is >> dummy >> n;
    face_.reserve(n);
    MeshTriangle* f;
    for (i=0; i<n; i++)
    {
        is >> Id;
        f = dynamic_cast<MeshTriangle*>(getDPMBase()->wallHandler.getObjectById(Id));
        face_.push_back(f);
    }

    // Reading in the edge information
    is >> dummy >> n;
    edge_.reserve(n);
    for (i=0; i<n; i++)
    {
        Edge e;
        // the vertex ids
        is >> dummy >> e.vertexId[0];
        e.vertex[0] = getDPMBase()->particleHandler.getObjectById(e.vertexId[0]);
        is >> dummy >> e.vertexId[1];
        e.vertex[1] = getDPMBase()->particleHandler.getObjectById(e.vertexId[1]);

        is >> dummy >> e.faceVertexId[0];
        if (e.faceVertexId[0] != UINT_MAX) {
            e.faceVertex[0] = getDPMBase()->particleHandler.getObjectById(e.faceVertexId[0]);
        } else {
            e.faceVertex[0] = nullptr;
        }
        is >> dummy >> e.faceVertexId[1];
        if (e.faceVertexId[1] != UINT_MAX) {

            e.faceVertex[1] = getDPMBase()->particleHandler.getObjectById(e.faceVertexId[1]);
        } else {
            e.faceVertex[1] = nullptr;
        }
        is >> dummy >> Id;
        if (Id != UINT_MAX) {
            e.face[0] = dynamic_cast<MeshTriangle*>(getDPMBase()->wallHandler.getObjectById(Id));
        } else {
            e.face[0] = nullptr;
        }
        is >> dummy >> Id;
        if (Id != UINT_MAX) {
            e.face[1] = dynamic_cast<MeshTriangle*>(getDPMBase()->wallHandler.getObjectById(Id));
        } else {
            e.face[1] = nullptr;
        }

        is >> dummy >> e.faceInitialArea[0];
        is >> dummy >> e.faceInitialArea[1];
        for (unsigned int i=0; i<8; i++){
            is >> dummy >> e.uPre[i];
        }
        
        // initial state
        is >> dummy >> e.initialState;
        is >> dummy >> e.initialLength;
        
        if (!e.face[0] || !e.face[1])
        {
            is >> dummy >> dummy;
        }
        else
        {
            is >> dummy >> e.initialSineHalfTheta;
        }
        is >> dummy >> e.effectiveMass;
        
        e.checkActive();
        edge_.push_back(e);
    }
    
    // Reading in some general attributes
    is >> dummy >> particleRadius_; // Mdouble

    is >> dummy >> Kn_; // Mdouble
    is >> dummy >> critDampCoeff_; //Mdouble

    is >> dummy >> Ke_; // Mdouble
    is >> dummy >> Kd_; // Mdouble
    
    is >> dummy >> thickness_;
    
    is >> dummy >> bendingAreaConstant_;
    
    is >> dummy >> Id;
    membraneSpecies_ = getDPMBase()->speciesHandler.getObjectById(Id);
    is >> dummy >> Id;
    membraneParticleSpecies_ = getDPMBase()->speciesHandler.getObjectById(Id);

    updateFaceNeighbors();
    
    for (auto& e: edge_)
    {
        e.checkActive();
    }
}

/*!
 * \param[in] os The output stream where the Membrane must be written
 *  to, usually a restart file.
 */
void Membrane::write(std::ostream& os) const
{

    // Write out the vertex-Particle Ids
    os << " vertexParticle " << vertexParticleId_.size();
    for (unsigned int i: vertexParticleId_)
    {
        os << " " << i;
    }    
    os << " vertexInitId_ " << vertexInitId_;


    // Write out the face-wall Ids
    os << " face " << face_.size();
    for (auto f: face_)
    {
        os << " " << f->getId();
    }

    // Write out edge information
    os << " edge " << edge_.size();
    for (auto e: edge_)
    {
        // the vertex ids
        os << " v3 " << e.vertexId[0];
        os << " v4 " << e.vertexId[1];

        
        os << " v1 " << e.faceVertexId[0];    
        os << " v2 " << e.faceVertexId[1];

        if (e.face[1] != nullptr){
            os << " f1 " << e.face[0]->getId();
        } else {
            os << " f1 " << UINT_MAX;
        }
        if (e.face[1] != nullptr){
            os << " f2 " << e.face[1]->getId();
        } else {
            os << " f2 " << UINT_MAX;
        }

        os << " area1 " << e.faceInitialArea[0];
        os << " area2 " << e.faceInitialArea[1];
        for (unsigned int i=0; i<8; i++){
            os << " uPre " << e.uPre[i];
        }
        // initial state
        os << " initState " << e.initialState;
        os << " initDist " << e.initialLength;
        os << " initialSineHalfTheta " << e.initialSineHalfTheta;
        os << " effectiveMass " << e.effectiveMass;
    }

    // Write out general properties
    os << " particleRadius " << particleRadius_; // Mdouble

    os << " Kn " << Kn_; // Mdouble
    os << " critDampCoeff " << critDampCoeff_; //Mdouble
    
    os << " Ke " << Ke_; // Mdouble
    os << " Kd " << Kd_; // Mdouble
    
    os << " thickness " << thickness_;
    
    os << " bendingAreaConstant " << bendingAreaConstant_;
    
    os << " membraneSpecies " << membraneSpecies_->getId();
    os << " membraneParticleSpecies " << membraneParticleSpecies_->getId();
}

/*!
 * \return The string "Membrane".
 */
std::string Membrane::getName() const
{
    return "Membrane";
}

void Membrane::setKnAndCrittDampCoeff(Mdouble Kn, Mdouble critDampCoeff) { 
    Kn_ = Kn; 
    critDampCoeff_ = critDampCoeff;
}

void Membrane::setElasticModulusAndThickness(Mdouble E, Mdouble thickness)
{
    setThickness(thickness);
    Kn_ = E * thickness_/(sqrt(3)*(1-1/3.0));
}

void Membrane::setSpringConstant(Mdouble k)
{
    logger.assert_debug(k >= 0, "Error in Membrane::setSpringConstant: The constant has to be greater than or equal to 0.");
    Kn_ = k;
}

void Membrane::setCriticalDampingCoefficient(Mdouble coeff)
{ 
    critDampCoeff_ = coeff; 
}

/*!
 * \brief Set the parameters needed for the bending forces
 */
void Membrane::setKeAndKd(Mdouble Ke, Mdouble Kd)
{ 
    Ke_ = Ke; 
    Kd_ = Kd;
}

void Membrane::setThickness(Mdouble thickness)
{ 
    logger.assert_debug(thickness > 0, "Error in Membrane::setThickness: The thickness has to be greater than 0.");
    thickness_ = thickness; 
}

void Membrane::setBendingAreaConstant(bool areaConstant)
{ 
    bendingAreaConstant_ = areaConstant; 
}

void Membrane::setParticleRadius(Mdouble radius)
{
    logger.assert_debug(radius > 0, "Error in Membrane::setParticleRadius: Radius has to be greater than 0.");
    particleRadius_ = radius;
}

void Membrane::saveVertexPositions(std::ostream& os)
{
    #ifdef MERCURY_USE_MPI
    if (PROCESSOR_ID==0)
    {
        os << " nMembraneParticles " << vertexParticleId_.size();
        os << " membraneParticlePositions ";
    }
    
    Vec3D pos, pos1;
    Mdouble active;
    BaseParticle* p0;
    for (unsigned int i=0; i < vertexParticleId_.size(); i++)
    {
        pos = Vec3D(0,0,0);
        p0 = getDPMBase()->particleHandler.getObjectById(vertexParticleId_[i]);
        if (p0 && !p0->isMPIParticle())
        {
            pos = p0->getPosition();
            active = getMPISum(1.0);
        }
        else
        {
            active = getMPISum(0.0);
        }
        
        pos1 = getMPISum(pos);
        if (PROCESSOR_ID==0 && active > 0.1)
        {
            os << pos1 << " ";
        }
    }    
    if (PROCESSOR_ID==0)
    {
        os << "\n";
    }
    #else
    os << " nMembraneParticles " << vertexParticleId_.size();
    os << " membraneParticlePositions ";
    BaseParticle* p0;
    for (unsigned int i=0; i < vertexParticleId_.size(); i++)
    {
        p0 = getDPMBase()->particleHandler.getObjectById(vertexParticleId_[i]);
        os << p0->getPosition() << " ";
    }
    #endif
}


void Membrane::loadVertexPositions(std::istream& is)
{
    unsigned int nVertices, i;
    std::string dummy;
    Vec3D pos;
    is >> dummy >> nVertices;
    is >> dummy;

    #ifdef MERCURY_USE_MPI
    if (PROCESSOR_ID==0)
    #endif
    logger.assert_always(nVertices==vertexParticleId_.size(),
                  "Error in Membrane::loadVertexPositions: The number of saved Particles differs from the number of actual particles");
    for (i=0; i<nVertices; i++)
    {
        is >> pos;
        if (vertexParticle_[i]) 
        {
            vertexParticle_[i]->setPosition(pos);
        }
    }

    // #ifdef MERCURY_USE_MPI
    // if (PROCESSOR_ID==0)
    // #endif
    // logger(INFO, "Loaded % membrane particles", nVertices);
}

/*!
 * \param[in] unsigned int d
 * \details This function saves the current geometry of the membrane to a file 
 * named membrane.$d.off.
 */
void Membrane::saveAsOFF(unsigned int d)
{
    std::stringstream lastName("");
    lastName << "membrane." << d << ".off";
    
    std::ofstream membraneOFF;
    membraneOFF.open(lastName.str(), std::ios_base::out);
    
    membraneOFF << "OFF\n";
    membraneOFF << "# membrane.off\n";
    
    // Number of  Vertices Faces and Edges (edges are optional)
    membraneOFF << vertexParticle_.size() << " " << face_.size() << " 0" << std::endl;
    
    // Print the coordinates of all vertices
    for (auto p: vertexParticle_)
    {
        membraneOFF << p->getPosition() << std::endl;
    }
    
    // Print all faces
    for (auto f: face_)
    {
        membraneOFF << "3 ";
        for (auto id: f->getVertexIds())
        {
            membraneOFF << id - vertexInitId_ << " ";
        } 
        membraneOFF <<  std::endl;
    }
}

/*!
 * \param[in] fileName The name of the file where the data is written to.
 * \details This function writes the current geometry of the Membrane to a
 * STL file.
 */
void Membrane::saveAsSTL(std::string fileName)
{
    std::ofstream os;
    os.open(fileName, std::ios_base::out);
    
    if (!os.is_open())
        logger(ERROR, "Membrane::saveAsSTL Cannot open file %.", fileName);
        
    os << "solid Membrane" << std::endl;
    
    std::array<Vec3D,3> vertexPositions;
    for (auto face: face_)
    {
        os << "  facet normal " << face->getFaceNormal() << std::endl;
        os << "    outer loop" << std::endl;
        vertexPositions = face->getVertices();
        for (auto pos: vertexPositions)
        {
            os << "      vertex " << pos << std::endl;
        }
        os << "    endloop" << std::endl;
        os << "  endfacet" << std::endl;
    }
    
    os << "endsolid Membrane" << std::endl;
}

/*!
 * \param[in] p0 BaseParticle that will be used for the vertex particles
 * \param[in] fileName The name of the file to load.
 * \details This function loads the geometry from the STL file and calls the function
 * buildMesh.
 */
void Membrane::loadFromSTL(BaseParticle &p0, std::string fileName)
{
    loadFromSTL(p0, fileName, 1e-8);
}

/*!
 * \param[in] p0 BaseParticle that will be used for the vertex particles
 * \param[in] fileName The name of the file to load.
 * \param[in] eps The distance between two points below which they are assumed
 * to be the same.
 * \details This function loads the geometry from the STL file and calls the function
 * buildMesh.
 */
void Membrane::loadFromSTL(BaseParticle& p0, std::string fileName, Mdouble eps)
{
    std::vector<Vec3D> vertexPositions;
    std::vector<unsigned int> edgeVertices;
    std::vector<unsigned int> faceVertices;
    
    
    std::ifstream is;
    is.open(fileName, std::ios_base::in);
    
    if (!is.is_open())
        logger(ERROR, "Membrane::loadFromSTL Cannot open file %.", fileName);
        
    std::string dummy;
    unsigned int i;
    Vec3D pos;
    
    // solid name
    is >> dummy >> dummy;
    // endsolid or facet
    is >> dummy;
    while ( dummy.compare("endsolid") )
    {
        // normal Vec3D
        is >> dummy >> pos;
        // outer loop
        is >> dummy >> dummy;
        for (i=0; i<3; i++)
        {
            // vertex pos
            is >> dummy >> pos;
            faceVertices.push_back(addVertex(vertexPositions, pos, eps));
        }
        // endloop
        is >> dummy;
        // endfacet
        is >> dummy;
        // endsolid or facet
        is >> dummy;
    }
    is.close();
    
    // Create the neccesary edges;
    std::map<std::pair<unsigned int, unsigned int>, int> existingEdges;
    for (i=0; i<faceVertices.size()/3; i++)
    {
        addEdge(edgeVertices, existingEdges, faceVertices[3*i+0], faceVertices[3*i+1]);
        addEdge(edgeVertices, existingEdges, faceVertices[3*i+1], faceVertices[3*i+2]);
        addEdge(edgeVertices, existingEdges, faceVertices[3*i+2], faceVertices[3*i+0]);
    }
    
    // logger(INFO, "Got % vertices, % edges % faces", vertexPositions.size(), edgeVertices.size()/2, faceVertices.size()/3);
    buildMesh(p0, vertexPositions, edgeVertices, faceVertices);
}

/*!
 * \param[in] vertices The vector, where the new position should be added
 * \param[in] pos The new position to add to the vector
 * \param[in] eps The distance between two points below which they are assumed
 * to be the same.
 * \returns the index of the position in the vector vertices.
 * \details This function adds a given position pos to the vector vertices, if 
 * it is not already contained.
 */
unsigned int Membrane::addVertex(std::vector<Vec3D>& vertices, Vec3D pos, Mdouble eps)
{
    unsigned int i;
    for (i=0; i<vertices.size(); i++)
    {
        if ((vertices[i]-pos).getLength() < eps)
        {
            return i;
        }
    }
    
    vertices.push_back(pos);
    
    return vertices.size() - 1;
}

/*!
 * \param[in] edges The vector, where the new edge should be added
 * \param[in] map A map used to dtermin, if the edge already exists
 * \param[in] v0 The index of the first point of the edge
 * \param[in] v0 The index of the second point of the edge
 * \returns True if a new edge was added, false if not.
 * \details This function adds a given edge to the vector edges, if 
 * it is not already contained.
 */
bool Membrane::addEdge(std::vector<unsigned int>& edges, std::map<std::pair<unsigned int, unsigned int>, int>& map, unsigned int v0, unsigned int v1)
{
    if ( map.find(std::make_pair(v0, v1)) == map.end() ){
        map[std::make_pair(v0,v1)] = 1;
        map[std::make_pair(v1,v0)] = 1;

        edges.push_back(v0);
        edges.push_back(v1);
        return true;
    }
    return false;
}

/*!
 * \param[in] p0 BaseParticle that will be used for the vertex particles
 * \param[in] vertexPositions The position of the vertices
 * \param[in] edgeVertices A vector of indices, where two consecutive values define
 * the vertex particles connected by the edge
 * \param[in] faceVertices A vector of indices, where three consecutive values define
 * the vertex particles making up the triangle face
 * \details This function creates the actual mesh including the vertex particles
 * and wall elements. The mass of the triangles is set to match the actual weight
 * of the cooresponding triangle. The density of the particle species is set, so
 * that the mass of a vertex particle is equal to 5/3 times the mass of a triangle.
 */
void Membrane::buildMesh(BaseParticle& p0, std::vector<Vec3D> vertexPositions, std::vector<unsigned int> edgeVertices, std::vector<unsigned int> faceVertices)
{    
    unsigned int triangleCount = faceVertices.size() / 3;
    unsigned int edgeCount = edgeVertices.size() / 2;
    // #ifdef MERCURY_USE_MPI
    // if (PROCESSOR_ID==0)
    // #endif
    // logger(INFO, "Creating membrane with % vertices and % faces.", vertexPositions.size(), triangleCount);

    // Some helper variables
    unsigned int i;
    Vec3D v;

    // Calculate the density of the Membrane particles to approximately get the membrane weight
    // Note this should be changed. But for the beginning it is a way of setting an approx. okej density for the particles.
    i = 0;
    Mdouble membraneDensity = membraneSpecies_->getDensity();
    

    // Create particles corresponding to the Vertex Positions
    createVertexParticles(p0, vertexPositions);
    
    Mdouble averageTriangleMass = 0;
    Mdouble triangleMass = 0;
    // Create all faces with their initial positions
    face_.reserve(triangleCount);
    for (i = 0; i < triangleCount; i++)
    {
        // It seems copyAndAddObject does not actually copy in MPI mode
        MeshTriangle f;
        f.setVertices(vertexPositions[faceVertices[3*i]],
                      vertexPositions[faceVertices[3*i+1]],
                      vertexPositions[faceVertices[3*i+2]]);
        f.setVertexIds(vertexParticleId_[faceVertices[3*i]],
                       vertexParticleId_[faceVertices[3*i+1]],
                       vertexParticleId_[faceVertices[3*i+2]]);
        f.setSpecies(membraneSpecies_);
        
        triangleMass = f.getArea()*thickness_*membraneDensity;
        averageTriangleMass += triangleMass;
        f.setMass(triangleMass);
        
        f.actionsAfterParticleGhostUpdate();
        
        face_.push_back(getDPMBase()->wallHandler.copyAndAddObject(f));
        
    }
    
    averageTriangleMass /=  triangleCount;
    Mdouble particleDensity = 5.0/3.0 * averageTriangleMass/(particleRadius_*particleRadius_*particleRadius_*4/3.0*constants::pi);
    membraneParticleSpecies_->setDensity(particleDensity);

    // Create the edges
    Mdouble invMass1, invMass2;
    
    Mdouble vertexInvMass = 1.0 / ( membraneParticleSpecies_->getDensity() * 4.0 / 3.0 * constants::pi * particleRadius_ * particleRadius_ * particleRadius_ );
    for (i=0; i<edgeCount; i++){
        Edge e;
        e.vertexId[0] = vertexParticleId_[edgeVertices[2*i+0]];
        e.vertexId[1] = vertexParticleId_[edgeVertices[2*i+1]];
        
        invMass1 = vertexInvMass;
        invMass2 = vertexInvMass;
        
        e.vertex[0] = getDPMBase()->particleHandler.getObjectById(e.vertexId[0]);
        e.vertex[1] = getDPMBase()->particleHandler.getObjectById(e.vertexId[1]);
        
        // e.vertex[1] = vertexParticle_[edgeVertices[2*i+1]];
        e.initialState = vertexPositions[edgeVertices[2*i+1]]-vertexPositions[edgeVertices[2*i+0]];
        e.initialLength = e.initialState.getLength();
        e.effectiveMass = 1/(invMass1 + invMass2);
        
        // e.checkActive();
        edge_.push_back(e);
    }
    
    initializeEdgeBendingQuantities();
    
    for (i=0; i<edge_.size(); i++){
        edge_[i].checkActive();
    }
    updateFaceNeighbors();
    for (i=0; i<edge_.size(); i++)
    {
        edge_[i].calculateUPre(edge_[i].initialState, edge_[i].initialLength, edge_[i].faceInitialArea);
    }
}

/*!
 * \param[in] p0 The template particle used to create the vertex particles
 * \details This function creates a particle for every position in the vector vertexPositions
 */
void Membrane::createVertexParticles(BaseParticle& p0, std::vector<Vec3D> vertexPositions)
{
    
    vertexInitId_ = getDPMBase()->particleHandler.getNumberOfRealObjects();
    
    Vec3D pos;
    for (unsigned int i = 0; i < vertexPositions.size(); i++){
        pos = vertexPositions[i];
        // Set the position
        p0.setPosition(pos);

        // Add the particle and save a reference
        getDPMBase()->particleHandler.copyAndAddObject(p0);
        vertexParticleId_.push_back(vertexInitId_+i);
        vertexParticle_.push_back(getDPMBase()->particleHandler.getObjectById(vertexInitId_+i));
    }

}
/*!
 * \details This function initializes the quantities and pointer needed on the edges
 * for calculating the forces due to the mass spring system.
 */
void Membrane::initializeEdgeBendingQuantities()
{
    //set neighbours
    unsigned int i, j;
    
    for (auto& e: edge_)
    {
        e.face[0] = nullptr;
        e.face[1] = nullptr;
        e.faceVertexId[0] = UINT_MAX;
        e.faceVertexId[1] = UINT_MAX;
        e.faceVertex[0] = nullptr;
        e.faceVertex[1] = nullptr;
        e.faceInitialArea[0] = 0;
        e.faceInitialArea[1] = 0;
        e.initialSineHalfTheta = 0;
        i = 0;

        unsigned int edgeId0 = e.vertexId[0];
        unsigned int edgeId1 = e.vertexId[1];

        for (unsigned int f=0; f<face_.size(); f++)
        {
            MeshTriangle* face0 = face_[f];
            std::array<unsigned int, 3> particleIds = face0->getVertexIds();
            j = UINT_MAX;

            if ( particleIds[0] == edgeId0)
            {
                if (particleIds[1] == edgeId1) j = particleIds[2];
                if (particleIds[2] == edgeId1) j = particleIds[1];
            }
            else if ( particleIds[0] == edgeId1)
            {
                if (particleIds[1] == edgeId0) j = particleIds[2];
                if (particleIds[2] == edgeId0) j = particleIds[1];
            }
            else if (   (particleIds[1] == edgeId1 && particleIds[2] == edgeId0)
                      ||(particleIds[1] == edgeId0 && particleIds[2] == edgeId1))
            {
                j = particleIds[0];
            }

            if (j!=UINT_MAX)
            {
                e.faceVertexId[i] = j;
                e.faceVertex[i] = getDPMBase()->particleHandler.getObjectById(j);
                e.face[i] = face0;
                e.faceInitialArea[i] = face0->getArea();
                if (i==1)
                {
                    e.initialSineHalfTheta = e.getSineHalfTheta();
                }
                i += 1;
            }
        }
    }
}

/*!
 * \details This function assigns the coorect neighbors to the wall elements.
 */
void Membrane::updateFaceNeighbors()
{
    //set neighbours
    unsigned int i, j, k, l;

    for (i =0 ; i<face_.size(); i++)
    {
        for (j=0; j<3; j++)
        {
            face_[i]->neighbor[j] = nullptr;
            std::vector<unsigned int> a;
            face_[i]->vertexNeighbors.push_back(a);
        }
    }
    for (i =0 ; i<face_.size(); i++)
    {
        MeshTriangle* face0 = face_[i];
        std::array<unsigned int, 3> particleIds0 = face0->getVertexIds();
        for (j=i+1 ; j<face_.size(); j++)
        {
            MeshTriangle* face1 = face_[j];
            std::array<unsigned int, 3> particleIds1 = face1->getVertexIds();
            if (particleIds0[0] == particleIds1[0])
            {
                if (particleIds0[1] == particleIds1[2])
                { //edge 0=2
                    face0->neighbor[0] = &*face1;
                    face1->neighbor[2] = &*face0;
                }
                else if (particleIds0[2] == particleIds1[1])
                { //edge 2=0
                    face0->neighbor[2] = &*face1;
                    face1->neighbor[0] = &*face0;
                }
            }
            else if (particleIds0[0] == particleIds1[1])
            {
                if (particleIds0[1] == particleIds1[0])
                { //edge 0=0
                    face0->neighbor[0] = &*face1;
                    face1->neighbor[0] = &*face0;
                }
                else if (particleIds0[2] == particleIds1[2])
                { //edge 2=1
                    face0->neighbor[2] = &*face1;
                    face1->neighbor[1] = &*face0;
                }
            }
            else if (particleIds0[0] == particleIds1[2])
            {
                if (particleIds0[1] == particleIds1[1])
                { //edge 0=1
                    face0->neighbor[0] = &*face1;
                    face1->neighbor[1] = &*face0;
                }
                else if (particleIds0[2] == particleIds1[0])
                { //edge 2=2
                    face0->neighbor[2] = &*face1;
                    face1->neighbor[2] = &*face0;
                }
            }
            else if (particleIds0[1] == particleIds1[0])
            {
                if (particleIds0[2] == particleIds1[2])
                { //edge 1=2
                    face0->neighbor[1] = &*face1;
                    face1->neighbor[2] = &*face0;
                }
            }
            else if (particleIds0[1] == particleIds1[1])
            {
                if (particleIds0[2] == particleIds1[0])
                { //edge 1=0
                    face0->neighbor[1] = &*face1;
                    face1->neighbor[0] = &*face0;
                }
            }
            else if (particleIds0[1] == particleIds1[2])
            {
                if (particleIds0[2] == particleIds1[1])
                { //edge 1=1
                    face0->neighbor[1] = &*face1;
                    face1->neighbor[1] = &*face0;
                }
            }

        }
    }

    for (i =0 ; i<face_.size(); i++)
    {
        MeshTriangle* face0 = face_[i];
        std::array<unsigned int, 3> particleIds0 = face0->getVertexIds();
        for (j=i+1 ; j<face_.size(); j++)
        {
            MeshTriangle* face1 = face_[j];
            std::array<unsigned int, 3> particleIds1 = face1->getVertexIds();
            for (k=0; k<3; k++)
            {
                for (l=0; l<3; l++)
                {
                    if( particleIds0[k] == particleIds1[l])
                    {
                        face0->vertexNeighbors[k].push_back(face1->getId());
                        face1->vertexNeighbors[l].push_back(face0->getId());
                    }
                }
            }
        }
    }
}

/*!
 * \details Calculates an updated radius for every vertex particle in order to 
 * account for the correct mass.
 */
void Membrane::adjustVertexParticleSize()
{
    std::vector<unsigned int> numberTriangles(vertexParticle_.size(), 0);
    std::vector<Mdouble> massParticle(vertexParticle_.size(), 0.0);
    
    Mdouble mass;
    for (auto f: face_)
    {
        // 3 Nodes per triangle
        mass = 1/f->getInvMass()/3;
        for (auto id: f->getVertexIds())
        {
            numberTriangles[id-vertexInitId_] += 1;
            massParticle[id-vertexInitId_] += mass;
        }
    }
    
    unsigned int i;
    Mdouble volume;
    Mdouble radius;
    Mdouble density = membraneParticleSpecies_->getDensity();
    for (i=0; i<vertexParticle_.size(); i++)
    {
        mass = massParticle[i]/numberTriangles[i];
        volume = mass / density;
        radius = std::cbrt(3.0/4.0 * volume/constants::pi);
        
        vertexParticle_[i]->setRadius(radius);
    }
}

/*!
 * \details Set the correct edge mass by taking the mass from the conencted vertices.
 */
void Membrane::updateEdgeMass()
{
    Mdouble invMass;
    for (auto e: edge_)
    {
        invMass = vertexParticle_[e.vertexId[0]-vertexInitId_]->getInvMass()
                + vertexParticle_[e.vertexId[1]-vertexInitId_]->getInvMass();
        
        if (invMass > 0)
        {
            e.effectiveMass = 1/invMass;
        }
        else
        {
            e.effectiveMass = 0;
        }
        
        e.checkActive();
    }
}

/*!
 * \details This function iterates through all edges to calculate the forces due
 * to the mass spring system. Ideally this function is called from the computeAllForces
 * function of DPMBase.
 */
void Membrane::computeAdditionalForces(){
    // Note: I can use omp on add force, because this function does check for omp, or not?
    // #pragma omp parallel num_threads(getNumberOfOMPThreads())
    {
        // Calculate the forces due to the bonds
        // #pragma omp for
        for (auto e: edge_){
            e.applyForce(Kn_, critDampCoeff_, Ke_, Kd_, bendingAreaConstant_);
        }
    }
}

/*!
 * \param[in] pressure Value of the surface pressure acing on the membrane [Pa]
 * \details This function applies a given pressure to the surface by calling the
 * function applyPressure of a MeshTriangle.
 */
void Membrane::applyPressure(Mdouble pressure){
    // Set the pressureforce to 0
    logger(DEBUG,"Calculating pressure force.");

    // #pragma omp parallel num_threads(getNumberOfOMPThreads())
    {
    // This should also be OMP parralelizable, because addForce pays attention
    // to OMP
        // #pragma omp for
        for (unsigned int i=0; i < face_.size(); i++ ){
            face_[i]->applyPressure(pressure);
        }
    }
}

/*!
 * \param[in] id The id of the removed particle
 * \details This function handles the removal of a particle with the given id from
 * the particle handler. If the particle was contained in an edge, the edge is disabled.
 * This function is potentially usable for MPI parallelization
 */
void Membrane::handleParticleRemoval(unsigned int id)
{
    for (auto& e: edge_)
    {
        e.handleParticleRemoval(id);
    }
    
    // Now find the corresponding vertexParticle_
    unsigned i = id - vertexInitId_;
    if (vertexParticleId_[i] == id)
    {
        vertexParticle_[i] = nullptr;
    } else
    {
        // Need to search for it
        for (i = 0; i < vertexParticleId_.size(); i++)
        {
            if (vertexParticleId_[i]==id)
            {
                // now i is the coorect index
                vertexParticle_[i] = nullptr;
                break;
            }
        }
    }
}

/*!
 * \param[in] id The id of the added particle
 * \param[in] p Pointer to the added particle
 * \details This function handles the addition of a particle with the given id to
 * the particle handler. If the particle is contained in an edge, the edge is potentially.
 * reactivated. This function is potentially usable for MPI parallelization
 */
void Membrane::handleParticleAddition(unsigned int id, BaseParticle* p)
{
    for (auto& e: edge_)
    {
        e.handleParticleAddition(id, p);
    }
    
    // Now find the corresponding vertexParticle_
    unsigned i = id - vertexInitId_;
    if (vertexParticleId_.size() > i && vertexParticleId_[i] == id)
    {
        vertexParticle_[i] = p;
    } else
    {
        // Need to search for it
        for (i = 0; i < vertexParticleId_.size(); i++)
        {
            if (vertexParticleId_[i]==id)
            {
                // now i is the correct index
                vertexParticle_[i] = p;
                break;
            }
        }
    }
}

//////////////////////////////////////////////////////////////////////////////
// from here on there are sub member functions only

/*!
 * \param[in] Kn The spring constant used for stretching.
 * \param[in] critDampCoeff The critical damping coefficient used for stretching.
 * \param[in] Ke The spring constant used for bending.
 * \param[in] Kn The dissipation constant used for bending.
 * \details This function applies the forces due to the mass spring system.
 */
void Membrane::Edge::applyForce(Mdouble Kn, Mdouble critDampCoeff, Mdouble Ke, Mdouble Kd, bool bendingAreaConstant)
{
    applyStretchForce(Kn, critDampCoeff);
    applyBendForce(Ke, Kd, bendingAreaConstant);
}

/*!
 * \param[in] Kn The spring constant used for stretching.
 * \param[in] critDampCoeff The critical damping coefficient used for stretching.
 * \details This function applies the forces due to stretching of the edge
 */
void Membrane::Edge::applyStretchForce(Mdouble Kn, Mdouble critDampCoeff)
{
    if ( !isStretchActive ) return;
    // First stretching
    Vec3D distanceVec = vertex[1]->getPosition() - vertex[0]->getPosition();
    Mdouble distance = distanceVec.getLength();
    Vec3D normal = distanceVec/distance;
    
    if (distance<0.1*initialLength)
    {
        logger(WARN, "Edge length critical: % of original length", distance/initialLength);
    }
    
    Vec3D normalDisplacement = normal*(initialLength-distance);

    Vec3D force = Vec3D(0,0,0);

    // Normal force
    force += normalDisplacement * Kn;

    // Damping
    Vec3D normalRelativeVelocity = normal* Vec3D::dot(normal, vertex[1]->getVelocity() - vertex[0]->getVelocity());
    Mdouble dissipationCoefficientN = 2*critDampCoeff * sqrt(effectiveMass * Kn);
    force += -dissipationCoefficientN*normalRelativeVelocity;

    // Normal force and Tangential Force
    vertex[0]->addForce(-force);
    vertex[1]->addForce(force);
}

/*!
 * \param[in] Ke The spring constant used for bending.
 * \param[in] Kn The dissipation constant used for bending.
 * \details This function applies the forces due to bending of the triangles
 * connected by the edge. The calculation is based on the model described in
 * doi:10.1109/SIBGRAPI.2014.20.
 */
void Membrane::Edge::applyBendForce(Mdouble Ke, Mdouble Kd, bool bendingAreaConstant)
{
    // Calculation according to doi:10.1109/SIBGRAPI.2014.20
    // Now bending
    // normal / face area
    if ( !isBendActive )
    {
        return;
    }

    unsigned int i;

    std::array<Vec3D, 4> u;
    std::array<Vec3D, 4> force;

    Vec3D initialState1 = vertex[1]->getPosition()-vertex[0]->getPosition();
    Mdouble initialLength1 = initialState1.getLength();

    // Do adapt when triangle changes
    std::array<Mdouble, 2> faceArea;
    if (!bendingAreaConstant)
    {
        faceArea[0] = face[0]->getArea();
        faceArea[1] = face[1]->getArea();
    }
    else
    {
        faceArea[0] = faceInitialArea[0];
        faceArea[1] = faceInitialArea[1];
    }
    calculateUPre(initialState1, initialLength1, faceArea);
    
    
    for (i=0; i<4; i++)
    {
        u[i] = face[0]->getFaceNormal() * uPre[i] +  face[1]->getFaceNormal() * uPre[4+i];
    }
    Mdouble sineHalfTheta = getSineHalfTheta();
    Mdouble prefactor = Ke*initialLength1*initialLength1 / (faceArea[0]+faceArea[1]) * ( sineHalfTheta - initialSineHalfTheta );

    for (i=0; i<4; i++)
    {
        force[i] = prefactor*u[i];
    }
    
    Mdouble dThetadT = Vec3D::dot(u[0], faceVertex[0]->getVelocity())
                     + Vec3D::dot(u[1], faceVertex[1]->getVelocity())
                     + Vec3D::dot(u[2], vertex[0]->getVelocity())
                     + Vec3D::dot(u[3], vertex[1]->getVelocity());

    prefactor = -Kd*initialLength1*dThetadT;
    for (i=0; i<4; i++)
    {
        force[i] += prefactor*u[i];
    }

    faceVertex[0]->addForce(force[0]);
    faceVertex[1]->addForce(force[1]);
    vertex[0]->addForce(force[2]);
    vertex[1]->addForce(force[3]);

}


void Membrane::Edge::calculateUPre(Vec3D& state, Mdouble& length, std::array<Mdouble, 2>& faceArea)
{
    if ( face[0]==nullptr || face[1]==nullptr)
    {
        for (unsigned int i=0; i<8; i++)
        {
            uPre[i] = 0;
        }
        return;
    }

    uPre[0] = length/faceArea[0];
    uPre[1] = 0;
    uPre[2] =   Vec3D::dot(state, faceVertex[0]->getPosition()-vertex[1]->getPosition()) / length / faceArea[0];
    uPre[3] = - Vec3D::dot(state, faceVertex[0]->getPosition()-vertex[0]->getPosition()) / length / faceArea[0];

    uPre[4] = 0;
    uPre[5] = length/faceArea[1];

    uPre[6] = + Vec3D::dot(state, faceVertex[1]->getPosition()-vertex[1]->getPosition()) / length / faceArea[1];
    uPre[7] = - Vec3D::dot(state, faceVertex[1]->getPosition()-vertex[0]->getPosition()) / length / faceArea[1];
}

/*!
 * \details This function calculates the sinus of half the angle between the 
 * two connected triangles.
 */
Mdouble Membrane::Edge::getSineHalfTheta()
{
    Mdouble dotProduct = Vec3D::dot(face[0]->getFaceNormal(), face[1]->getFaceNormal());
    Mdouble sineHalfTheta;
    if (dotProduct>=1)
    {
        sineHalfTheta = 0;
    }
    else
    {
        sineHalfTheta = sqrt((1-dotProduct)/2);
    }
    
    Vec3D initialState1 = vertex[1]->getPosition()-vertex[0]->getPosition();
    if (Vec3D::dot(Vec3D::cross(face[0]->getFaceNormal(), face[1]->getFaceNormal()), initialState1) < 0)
    {
        return -sineHalfTheta;
    }
    return sineHalfTheta;
}

/*!
 * \param[in] id The id of the removed particle
 * \details This function handles the removal of a particle with the given id from
 * the particle handler. If the particle was contained in the edge, the edge is disabled.
 * This function is potentially usable for MPI parallelization
 */
void Membrane::Edge::handleParticleRemoval(unsigned int id)
{
    unsigned int i;
    for (i=0; i<2; i++)
    {
        if (vertexId[i] == id)
        {
            vertex[i] = nullptr;
            isStretchActive = false;
            isBendActive = false;
        }
        if (faceVertexId[i] == id)
        {
            faceVertex[i] = nullptr;
            isBendActive = false;
        }
    }
}

/*!
 * \param[in] id The id of the added particle
 * \param[in] p Pointer to the added particle
 * \details This function handles the addition of a particle with the given id to
 * the particle handler. If the particle is contained in the edge, the edge is potentially.
 * reactivated. This function is potentially usable for MPI parallelization
 */
void Membrane::Edge::handleParticleAddition(unsigned int id, BaseParticle* p)
{
    unsigned int i;
    for (i=0; i<2; i++)
    {
        if (vertexId[i] == id)
        {
            vertex[i] = p;
            checkActive();
        }
        if (faceVertexId[i] == id)
        {
            faceVertex[i] = p;
            checkActive();
        }
    }
}

/*!
 * \details This function checks if the edge knows al required pointers to calculate
 * and apply stretch and/or bending forces.
 */
void Membrane::Edge::checkActive()
{
    isStretchActive = !(!vertex[0] || !vertex[1]) && effectiveMass > 0;
    isBendActive = isStretchActive && !(!faceVertex[0] || !faceVertex[1] || !face[0] || !face[1]);
}
