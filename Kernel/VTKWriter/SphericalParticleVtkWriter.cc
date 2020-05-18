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

#include "VTKWriter/SphericalParticleVtkWriter.h"
#include "DPMBase.h"

/*!
 * \details Writes all points ans cells to a file in the VTK format. The
 * filename is hard-coded in this method, and is based on the name of the
 * DPMBase and has a unique counter in it to ensure there are no two files
 * with the same name.
 */
void SphericalParticleVtkWriter::writeVTK() const
{
    std::fstream file = makeVTKFileWithHeader();
    long numberOfParticlesToWrite = std::count_if(handler_.begin(), handler_.end(),
                                                  [this](BaseParticle* p){return particleMustBeWritten(p);});
    file << "<Piece NumberOfPoints=\"" << numberOfParticlesToWrite
         << "\" NumberOfCells=\"" << 0 << "\">\n";
    writeVTKPositions(file);
    file << "<PointData  Vectors=\"vector\">\n";
    writeVTKVelocity(file);
    writeVTKRadius(file);
    writeVTKIndSpecies(file);
    writeExtraFields(file);

    file << "</PointData>\n";
    writeVTKFooterAndClose(file);
}

void SphericalParticleVtkWriter::writeVTKVelocity(std::fstream& file) const
{
    file << "  <DataArray type=\"Float32\" Name=\"Velocity\" NumberOfComponents=\"3\" format=\"ascii\">\n";
    // Add velocity
    for (const auto& p: handler_) {
#ifdef MERCURY_USE_MPI
        if (particleMustBeWritten(p))
#endif
        {
            file << ' ' << (float) p->getVelocity().X
                 << ' ' << (float) p->getVelocity().Y
                 << ' ' << (float) p->getVelocity().Z << '\n';
        }
    }
    file << "  </DataArray>\n";
}

/*!
 * Notice that we write GrainRadius in the file, since there is a bug in Paraview that defaults to the first
 * scalar-value in lexicographic order. We therefore need a description for the radius which starts with a letter before
 * N.
 * \param file The filestream to which the  radius must be written.
 */
void SphericalParticleVtkWriter::writeVTKRadius(std::fstream& file) const
{
    file << "  <DataArray type=\"Float32\" Name=\"Radius\" format=\"ascii\">\n";
    // Add radius
    for (const auto& p: handler_)
    {
#ifdef MERCURY_USE_MPI
        if (particleMustBeWritten(p))
      {
        file << '\t' << p->getRadius() << '\n';
      }
#else
        file << '\t' << p->getRadius() << '\n';//Radius
#endif
    }
    file << "  </DataArray>\n";
}
