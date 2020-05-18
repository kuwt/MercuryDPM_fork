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

#include "VTKWriter/InteractionVTKWriter.h"
#include "DPMBase.h"

void InteractionVTKWriter::writeVTK() const
{
    std::fstream file = makeVTKFileWithHeader();
    file << "<Piece NumberOfPoints=\"" << handler_.getNumberOfObjects() << "\" NumberOfCells=\"" << 0 << "\">\n";
    file << "<Points>\n";
    writeVTKPoints(file);
    file << "</Points>\n";
    file << "<PointData  Vectors=\"vector\">\n";
    writeVTKPointData(file);
    file << "</PointData>\n";
    writeVTKFooterAndClose(file);
}

void InteractionVTKWriter::writeVTKPoints(std::fstream& file) const
{
    file << "  <DataArray type=\"Float32\" Name=\"Position\" NumberOfComponents=\"3\" format=\"ascii\">\n";
    for (const auto& p: handler_)
    {
        file << '\t' << p->getContactPoint() << '\n';
    }
    file << "  </DataArray>\n";
}

void InteractionVTKWriter::writeVTKPointData(std::fstream& file) const
{
    file << "  <DataArray type=\"Float32\" Name=\"Normal\" NumberOfComponents=\"3\" format=\"ascii\">\n";
    // Add velocity
    for (const auto& p: handler_)
    {
        file << '\t' << p->getNormal() << '\n';
    }
    file << "  </DataArray>\n";
    file << "  <DataArray type=\"Float32\" Name=\"Overlap\" format=\"ascii\">\n";
    
    // Add overlap
    for (const auto& p: handler_)
    {
        file << '\t' << p->getOverlap() << '\n';
    }
    file << "  </DataArray>\n";
    file << "  <DataArray type=\"Float32\" Name=\"ContactRadius\" format=\"ascii\">\n";
    
    // Add radius
    for (const auto& p: handler_)
    {
        file << '\t' << p->getContactRadius() << '\n';
    }
    file << "  </DataArray>\n";
    file << "  <DataArray type=\"Float32\" Name=\"Force\" format=\"ascii\">\n";
    
    // Add species type
    for (const auto& p: handler_)
    {
        file << '\t' << p->getForce() << '\n';
    }
    file << "  </DataArray>\n";
    file << "  <DataArray type=\"Float32\" Name=\"TangentialOverlap\" format=\"ascii\">\n";
    
    // Add species type
    for (const auto& p: handler_)
    {
        file << '\t' << p->getTangentialOverlap() << '\n';
    }
    file << "  </DataArray>\n";
    file << "  <DataArray type=\"Float32\" Name=\"Torque\" format=\"ascii\">\n";
    
    // Add species type
    for (const auto& p: handler_)
    {
        file << '\t' << p->getTorque() << '\n';
    }
    file << "  </DataArray>\n";
    
    //check if this type of Interaction has extra fields
    if (handler_.getSize() != 0)
    {
        for (unsigned i = 0; i < handler_.getLastObject()->getNumberOfFieldsVTK(); i++)
        {
            file << "  <DataArray type=\"" << handler_.getLastObject()->getTypeVTK(i) << "\" Name=\""
                 << handler_.getLastObject()->getNameVTK(i) << "\" format=\"ascii\">\n";
            // Add species type
            for (const auto& p: handler_)
            {
                for (auto f : p->getFieldVTK(i))
                    file << '\t' << f << '\n';
            }
            file << "  </DataArray>\n";
        }
    }
}
