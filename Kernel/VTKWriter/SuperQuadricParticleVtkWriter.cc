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


#include "VTKWriter/SuperQuadricParticleVtkWriter.h"
#include "DPMBase.h"

void SuperQuadricParticleVtkWriter::writeVTK() const
{
    std::fstream file = makeVTKFileWithHeader();
    file << "<Piece NumberOfPoints=\"" << handler_.getSize() << "\" NumberOfCells=\"" << 0 << "\">\n";
    writeVTKPositions(file);
    file << "<PointData Scalars=\"speciesType\" Vectors=\"axesScales\" Tensors=\"orientationTensor\">\n";
    writeVTKIndSpecies(file);
    writeVTKSuperquadricGeometry(file);
    writeVTKOrientation(file);
    file << "</PointData>\n";
    writeVTKFooterAndClose(file);
}

void SuperQuadricParticleVtkWriter::writeVTKOrientation(std::fstream& file) const
{
    file << "  <DataArray type=\"Float32\" Name=\"orientationTensor\" NumberOfComponents=\"9\" format=\"ascii\">\n";
    // Add orientationTensor
    for (const auto& p: handler_)
    {
        SmallMatrix<3, 3> A;
        p->getOrientation().getRotationMatrix(A);
        file << '\t';
        for (unsigned i = 0; i < 3; ++i)
        {
            for (unsigned j = 0; j < 3; ++j)
            {
                file << A(i, j) << " ";
            }
        }
        file << '\n';
    }
    file << "  </DataArray>\n";
}

void SuperQuadricParticleVtkWriter::writeVTKSuperquadricGeometry(std::fstream& file) const
{
    file << "  <DataArray type=\"Float32\" Name=\"phiAndTheta\" NumberOfComponents=\"2\" format=\"ascii\">\n";
    // Add exponents eps1 and eps2. Note, they are defined in reverse order in Paraview!
    for (const auto& p: handler_)
    {
        file << '\t' << p->getExponentEps2() << " " <<  p->getExponentEps1()  << '\n';
    }
    file << "  </DataArray>\n";
    file << "  <DataArray type=\"Float32\" Name=\"axesScales\" NumberOfComponents=\"3\" format=\"ascii\">\n";
    // Add axesScales
    for (const auto& p: handler_)
    {
        file << '\t' << 2 * p->getAxes() << '\n';
    }
    file << "  </DataArray>\n";
}
