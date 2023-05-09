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

#ifndef MERCURYDPM_WEARABLETRIANGLEMESHWALL_H
#define MERCURYDPM_WEARABLETRIANGLEMESHWALL_H

#include "TriangleMeshWall.h"

class WearableTriangleMeshWall : public TriangleMeshWall
{
public:
    
    /*!
     * \brief Default constructor.
     */
    WearableTriangleMeshWall() = default;
    
    /*!
     * \brief Constructor creating triangles
     * @param points Vector of vertices
     * @param cells Vector of cells, each cell holds three indices to the points vector
     * @param species
     */
    WearableTriangleMeshWall(const std::vector<Vec3D>& points, const std::vector<std::array<unsigned, 3>>& cells, const ParticleSpecies* species = nullptr);

    /*!
     * \brief Constructor creating parallelogram shaped mesh with a given resolution in u- and v-direction.
     * @param P0, P1, P2 The three points, in clockwise order, forming the parallelogram (fourth point is implied).
     * @param resolutionU Maximum length of a segment in u-direction.
     * @param resolutionV Maximum length of a segment in v-direction.
     * @param species
     * @param periodicInU Whether or not the first and last column of vertices lie in a periodic boundary and should be set as periodic companions.
     * @param periodicInV Whether or not the first and last row of vertices lie in a periodic boundary and should be set as periodic companions.
     */
    WearableTriangleMeshWall(const Vec3D& P0, const Vec3D& P1, const Vec3D& P, Mdouble resolutionU, Mdouble resolutionV,
                             const ParticleSpecies* species = nullptr, bool periodicInU = false, bool periodicInV = false);
    
     /*!
      * \brief Copy constructor.
      * @param other The WearableTriangleMeshWall that must be copied.
      */
    WearableTriangleMeshWall(const WearableTriangleMeshWall& other);
    
    /*!
     * \brief Destructor.
     */
    ~WearableTriangleMeshWall() = default;
    
     /*!
      * \brief Copy assignment operator.
      * @param other The WearableTriangleMeshWall that must be copied.
      * @return
      */
    WearableTriangleMeshWall& operator=(const WearableTriangleMeshWall& other);
    
     /*!
      * \brief Wall copy method.
      * @return Pointer to the copy.
      */
    WearableTriangleMeshWall* copy() const override;
    
    /*!
     * \brief Reads a WearableTriangleMeshWall from an input stream, for example a restart file.
     */
    void read(std::istream& is) override;
    
    /*!
     * \brief Writes a WearableTriangleMeshWall to an output stream, for example a restart file.
     */
    void write(std::ostream& os) const override;
    
     /*!
      * \brief Returns the name of the object.
      * @return "WearableTriangleMeshWall"
      */
    std::string getName() const override;
    
    void computeWear() override;
    
    /*!
     * @param wearCoefficient Non-dimensionless wear constant.
     */
    void setWearCoefficient(Mdouble wearCoefficient);
    
    /*!
     * @param hardness Hardness of the wall.
     */
    void setHardness(Mdouble hardness);
    
    /*!
     * @param wearAcceleration Value to accelerate the wear process with, to reduce simulation time needed to get results.
     */
    void setWearAcceleration(Mdouble wearAcceleration);
    
private:
    /*!
     * \brief Proportionally assigns the debris located at a certain position on a triangle to the triangle vertices.
     * @param triangle Triangle on which the point lies.
     * @param position Point of debris.
     * @param debris The debris, which has a direction, i.e. the vector length is be the actual debris volume.
     * @param debrisContainer Vector storing the debris for all vertices (index wise related).
     */
    void storeDebris(const Triangle& triangle, const Vec3D& position, const Vec3D& debris, std::vector<Vec3D>& debrisContainer);
    
    /// Dimensionless wear coefficient.
    Mdouble wearCoefficient_;
    /// Hardness.
    Mdouble hardness_;
    /// Accelerates the wear process, to reduce simulation time needed to get results.
    Mdouble wearAcceleration_;
};


#endif //MERCURYDPM_WEARABLETRIANGLEMESHWALL_H
