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

#ifndef TriangulatedWall_H
#define TriangulatedWall_H

#include <vector>
#include "BaseWall.h"
#include "Math/Vector.h"
#include <array>

/*!
 * \class TriangulatedWall
 * \brief A TriangulatedWall is a triangulation created from a set of vertices and a n-by-3 connectivity matrix defining n faces.
 * \details It is initialised by a unstructured grid vtk file:
 * \code
 * TriangulatedWall w("TriangulatedWallSimple.vtk",speciesHandler.getLastObject());
 * \endcode
 *
 * The file consists of a set of vertices and a n-by-3 connectivity matrix defining n faces.
 * Three vertices form a face; the face normal is oriented such that the vertices are ordered in anticlockwise direction around the normal.
 * <img src="triangulatedWall.png" height="250px">
 *
 * Particles interact with a TriangulatedWall when they contact a face (from either side), a edge, or a vertex.
 *
 * For a demonstration on how to use this class, see \ref TriangulatedWallSelfTest (shown in the image below).
 * <img src="triangulatedWallsSelfTest.png" height="250px">
 */
class TriangulatedWall : public BaseWall
{
public:
    
    /**
     * Struct used to store the properties of a face needed for contact detection.
     * <img src="triangulatedWall.png" height="250px">
     */
    struct Face
    {
        ///defines the three vertices (anticlockwise direction around the normal)
        std::array<Vec3D*, 3> vertex;
        ///For each edge, stores the neighboring face (nullptr if none)
        std::array<Face*, 3> neighbor = {{nullptr}};
        ///For each edge, stores the vector normal to the face normal and the edge direction (vector between the vertices).
        std::array<Vec3D, 3> edgeNormal;
        ///face normal (vertices are ordered anticlockwise direction around the normal)
        Vec3D normal;
        
        ///computes the signed distance to the face in normal direction
        Mdouble getDistance(const Vec3D& otherPosition) const;
        
        ///Returns true if contact with the face exists, false if not. If contact exists, then the distance and normal is returned as well
        bool getDistanceAndNormal(const BaseParticle& p, Mdouble& distance, Vec3D& normal_return, Mdouble interactionRadius) const;
    };
    
    /*!
     * \brief Default constructor.
     */
    TriangulatedWall();
    
    /*!
     * \brief Copy constructor.
     */
    TriangulatedWall(const TriangulatedWall& other);
    
    /*!
     * \brief Constructor setting values.
     */
    TriangulatedWall(std::string filename, const ParticleSpecies* species);
    
    void readVTK(std::string filename);
    
    void writeVTK(VTKContainer& vtk) const override;
    
    /*!
     * \brief Destructor.
     */
    ~TriangulatedWall() override;
    
    /*!
     * Copy assignment operator.
     */
    TriangulatedWall& operator=(const TriangulatedWall& other);
    
    /*!
     * \brief Wall copy method. It calls the copy constructor of this Wall, useful for polymorphism
     */
    TriangulatedWall* copy() const override;
    
    /*!
     * \brief Compute the distance from the wall for a given BaseParticle and return if there is a collision. If there is a collision, also return the normal vector.
     */
    bool getDistanceAndNormal(const BaseParticle& p, Mdouble& distance, Vec3D& normal_return) const override;
    
    /*!
     * \brief Move the TriangulatedWall to a new position, which is a Vec3D from the old position.
     */
    void move(const Vec3D& move) override;
    
    /*!
     * \brief Reads an TriangulatedWall from an input stream, for example a restart file.
     */
    void read(std::istream& is) override;
    
    /*!
     * \brief Writes an TriangulatedWall to an output stream, for example a restart file.
     */
    void write(std::ostream& os) const override;
    
    /*!
     * \brief Returns the name of the object, here the string "TriangulatedWall".
     */
    std::string getName() const override;
    
    /*!
     * \brief Get the interaction between this TriangulatedWall and given BaseParticle at a given time.
     */
    BaseInteraction* getInteractionWith(BaseParticle* p, unsigned timeStamp,
                                                     InteractionHandler* interactionHandler) override;

private:
    /*!
     * stores the vertex coordinates
     */
    std::vector<Vec3D> vertex_;
    
    /*!
     * stores the face properties
     */
    std::vector<Face> face_;
    
};


#endif
