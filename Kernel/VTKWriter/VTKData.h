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

#ifndef MERCURY_VTKDATA_H
#define MERCURY_VTKDATA_H

#include <string>
#include <vector>
#include <fstream>
#include "Math/Vector.h"
#include <unordered_map>

class VTKData
{
public:
    VTKData() = default;
    ~VTKData() = default;
    
    /*!
     * \brief Adds a point to the points vector
     * @param point Point to be added
     */
    void addToPoints(Vec3D point);
    
    /*!
     * \brief Adds a value to the pointData values vector corresponding to the given key
     * @param key To which pointData values vector to add
     * @param value Value to be added
     */
    void addToPointData(const std::string& key, Mdouble value);
    
    /*!
     * \brief Adds a vector of indices to the connectivity vector
     * @param indices Vector of indices to be added
     */
    void addToConnectivity(const std::vector<size_t>& indices);
    
    /*!
     * \brief Adds a type to the types vector
     * @param type Type to be added
     */
    void addToTypes(int type);
    
    /*!
     * @return Vector of points
     */
    std::vector<Vec3D> getPoints() const;
    
    /*!
     * @return Map of pointData, each having a key and a vector of values
     */
    std::unordered_map<std::string, std::vector<Mdouble>> getPointData() const;
    
    /*!
     * @return Vector of vectors with indices
     */
    std::vector<std::vector<size_t>> getConnectivity() const;
    
    /*!
     * @return Vector of types
     */
    std::vector<int> getTypes() const;
    
    /*!
     * \brief Reserves additional memory for the points vector and optionally for the values vectors of the pointData
     * @param n Number of new points to reserve memory for
     * @param keys Keys for pointData values vectors to reserve memory for
     */
    void reservePoints(unsigned int n, const std::vector<std::string>& keys = {});
    
    /*!
     * \brief Reserves additional memory for connectivity and types vectors
     * @param n Number of new cells to reserve memory for
     */
    void reserveCells(unsigned int n);
    
    /*!
     * \brief Writes the data to a file with the given file name
     * @param fileName Full filename, e.g. Example.vtu
     */
    void writeVTKData(std::string fileName) const;

    void writeVTKDataFromVtkContainer(std::string fileName, const std::vector<Vec3D>& points, const std::vector<std::vector<double>>& triangleStrips);
    
private:
    std::fstream makeVTKFileWithHeader(std::string& fileName) const;
    void writePoints(std::fstream& file) const;
    void writePointData(std::fstream& file) const;
    void writeConnectivity(std::fstream& file) const;
    void writeOffsets(std::fstream& file) const;
    void writeTypes(std::fstream& file) const;
    
    /*!
     * The 3D positions
     */
    std::vector<Vec3D> points_;
    
    /*!
     * The point data, corresponding to each point
     */
    std::unordered_map<std::string, std::vector<Mdouble>> pointData_;
    
    /*!
     * The connectivity 2D vector, connecting points
     */
    std::vector<std::vector<size_t>> connectivity_;
    
    /*!
     * The types, corresponding to each row of the connectivity 2D vector
     */
    std::vector<int> types_;
};


#endif //MERCURY_VTKDATA_H
