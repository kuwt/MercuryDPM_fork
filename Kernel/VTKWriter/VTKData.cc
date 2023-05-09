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

#include "VTKData.h"

void VTKData::addToPoints(Vec3D point)
{
    points_.push_back(point);
}

void VTKData::addToPointData(const std::string& key, Mdouble value)
{
    pointData_[key].push_back(value);
}

void VTKData::addToConnectivity(const std::vector<size_t>& indices)
{
    connectivity_.push_back(indices);
}

void VTKData::addToTypes(int type)
{
    types_.push_back(type);
}

std::vector<Vec3D> VTKData::getPoints() const
{
    return points_;
}

std::unordered_map<std::string, std::vector<Mdouble>> VTKData::getPointData() const
{
    return pointData_;
}

std::vector<std::vector<size_t>> VTKData::getConnectivity() const
{
    return connectivity_;
}

std::vector<int> VTKData::getTypes() const
{
    return types_;
}

void VTKData::reservePoints(const unsigned int n, const std::vector<std::string>& keys)
{
    unsigned int nt = points_.size() + n;
    points_.reserve(nt);
    for (const std::string& key : keys)
        pointData_[key].reserve(nt);
}

void VTKData::reserveCells(unsigned int n)
{
    unsigned int nt = connectivity_.size() + n;
    connectivity_.reserve(nt);
    types_.reserve(nt);
}

void VTKData::writeVTKData(std::string fileName) const
{
    std::fstream file = makeVTKFileWithHeader(fileName);
    
    file << "<Piece NumberOfPoints=\"" + std::to_string(points_.size()) + "\" NumberOfCells=\"" + std::to_string(connectivity_.size()) + "\">\n";
    writePoints(file);
    
    writePointData(file);
    
    file << "<Cells>\n";
    writeConnectivity(file);
    writeOffsets(file);
    writeTypes(file);
    file << "</Cells>\n";
    
    file << "</Piece>\n";
    file << "</UnstructuredGrid>\n";
    file << "</VTKFile>\n";
    // Close output file
    file.close();
}

void VTKData::writeVTKDataFromVtkContainer(std::string fileName, const std::vector<Vec3D>& points, const std::vector<std::vector<double>>& triangleStrips)
{
    for (auto& p : points)
        addToPoints(p);
    for (auto& c : triangleStrips)
    {
        std::vector<size_t> cell;
        for (auto& cc : c)
            cell.push_back(static_cast<size_t>(cc));
        addToConnectivity(cell);
        addToTypes(6);
    }
    writeVTKData(fileName);
}

std::fstream VTKData::makeVTKFileWithHeader(std::string& fileName) const
{
    // Open output file
    std::fstream file;
    file.open(fileName.c_str(), std::ios_base::out);
    if (file.fail())
    {
        logger(WARN, "File % could not be opened", fileName);
    }
    
    // Write output file header
    file << "<?xml version=\"1.0\"?>\n";
    //file << "<!-- time " << handler_.getDPMBase()->getTime() << "-->\n"; //\todo JWB do we care about time?
    file << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
    file << "<UnstructuredGrid>\n";
    return file;
}

void VTKData::writePoints(std::fstream &file) const
{
    file << "<Points>\n";
    file << "  <DataArray type=\"Float32\" Name=\"Position\" NumberOfComponents=\"3\" format=\"ascii\">\n";
    
    for (const Vec3D p : points_)
    {
        file << "\t" << p << "\n";
    }
    
    file << "  </DataArray>\n";
    file << "</Points>\n";
}

void VTKData::writePointData(std::fstream &file) const
{
    file << "<PointData Vectors=\"vector\">\n";
    
    for (const auto& da : pointData_)
    {
        const std::string& name = da.first;
        const std::vector<Mdouble>& values = da.second;
        
        unsigned int numberOfComponents = points_.size() / values.size();
        file << "  <DataArray type=\"Float32\" Name=\"" << name << "\" NumberOfComponents=\"" << numberOfComponents << "\" format=\"ascii\">\n";
        
        int i = 0;
        while (i < values.size())
        {
            file << "\t";
            for (int j = 0; j < numberOfComponents; j++)
            {
                file << values[i++] << " ";
            }
            file << "\n";
        }
        
        file << "  </DataArray>\n";
    }
    
    file << "</PointData>\n";
}

void VTKData::writeConnectivity(std::fstream &file) const
{
    file << "  <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n";
    
    for (const std::vector<size_t>& vec : connectivity_)
    {
        file << "\t";
        for (const int i : vec)
        {
            file << i << " ";
        }
        file << "\n";
    }
    
    file << "  </DataArray>\n";
}

void VTKData::writeOffsets(std::fstream &file) const
{
    file << "  <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n";

    unsigned int offset = 0;
    for (const std::vector<size_t>& vec : connectivity_)
    {
        offset += vec.size();
        file << "\t" << offset << "\n";
    }
    
    file << "  </DataArray>\n";
}

void VTKData::writeTypes(std::fstream &file) const
{
    file << "  <DataArray type=\"UInt8\"  Name=\"types\" format=\"ascii\">\n";
    
    for (const int i : types_)
    {
        file << "\t" << i << "\n";
    }
    
    file << "  </DataArray>\n";
}