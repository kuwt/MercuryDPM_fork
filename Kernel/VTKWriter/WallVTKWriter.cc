//Copyright (c) 2013-2018, The MercuryDPM Developers Team. All rights reserved.
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

#include "VTKWriter/WallVTKWriter.h"
#include "DPMBase.h"


void WallVTKWriter::getVTKData(VTKContainer& vtk) const
{
    //set capacity of points and cells based on the previous time step
    static unsigned int capacityPoints = 0;
    static unsigned int capacityTriangleStrips = 0;
    
    vtk.triangleStrips.reserve(capacityPoints);
    vtk.points.reserve(capacityTriangleStrips);
    
    //add all wall data to the point and cell arrays
    for (const auto& w: handler_)
    {
        w->renderWall(vtk);
        logger(DEBUG, "points: %, cells: %", vtk.points.size(), vtk.triangleStrips.size());
    }
}

void WallVTKWriter::writeVTK() const
{
    if (PROCESSOR_ID!=0) return;
    std::fstream file = makeVTKFileWithHeader();
    
    VTKContainer vtk;
    getVTKData(vtk);
    file << "<Piece NumberOfPoints=\"" << vtk.points.size()
         << "\" NumberOfCells=\"" << vtk.triangleStrips.size()
         << "\">\n"
         << "<Points>\n";
    writeVTKPoints(file, vtk);
    file << "</Points>\n"
         << "<Cells>\n";
    writeVTKCells(file, vtk);
    file << "</Cells>\n"
         << "</Piece>\n"
         << "</UnstructuredGrid>\n"
         << "</VTKFile>\n";
    file.close();
}

void WallVTKWriter::writeVTKPoints(std::fstream& file, VTKContainer& vtk) const
{
    file << "  <DataArray type=\"Float32\" Name=\"Position\" NumberOfComponents=\"3\" format=\"ascii\">\n";
    for (const Vec3D& p : vtk.points)
    {
        file << '\t' << p << '\n';
    }
    file << "  </DataArray>\n";
}

void WallVTKWriter::writeVTKCells(std::fstream& file, VTKContainer& vtk) const
{
    file << "  <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n";
    for (const std::vector<double>& c : vtk.triangleStrips)
    {
        file << '\t';
        for (const double& i : c)
        {
            file << i << ' ';
        }
        file << '\n';
    }
    file << "  </DataArray>\n";
    file << "  <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n";
    unsigned count = 0;
    for (const std::vector<double>& c : vtk.triangleStrips)
    {
        count += c.size();
        file << '\t' << count << '\n';
    }
    file << "  </DataArray>\n";
    file << "  <DataArray type=\"UInt8\"  Name=\"types\" format=\"ascii\">\n";
    for (const std::vector<double>& c : vtk.triangleStrips)
    {
        if (c.front() == c.back())
        {
            //polygon
            file << "\t7\n";
        }
        else
        {
            //triangle strips
            file << "\t6\n";
        }
    }
    file << "  </DataArray>\n";
}
