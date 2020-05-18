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

#include "VTKWriter/ParticleVtkWriter.h"
#include "DPMBase.h"

void ParticleVtkWriter::writeVTKPositions(std::fstream& file) const
{
    file << "<Points>\n";
    file << "  <DataArray type=\"Float32\" Name=\"Position\" NumberOfComponents=\"3\" format=\"ascii\">\n";
    for (const auto& p: handler_)
    {
#ifdef MERCURY_USE_MPI
        if (particleMustBeWritten(p))
        {
            file << '\t' << p->getPosition() << '\n';
        }
#else
        file << '\t' << p->getPosition() << '\n';
#endif
    
    }
    file << "  </DataArray>\n";
    file << "</Points>\n";
}

void ParticleVtkWriter::writeVTKIndSpecies(std::fstream& file) const
{
    file << "  <DataArray type=\"Float32\" Name=\"speciesType\" format=\"ascii\">\n";
    // Add species type
    for (const auto& p: handler_)
    {
#ifdef MERCURY_USE_MPI
        if (particleMustBeWritten(p))
        {
            file << '\t' << p->getIndSpecies() << '\n';
        }
#else
        file << '\t' << p->getIndSpecies() << '\n';
#endif
    }
    file << "  </DataArray>\n";
}

void ParticleVtkWriter::writeExtraFields(std::fstream& file) const
{
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
                for (auto f : p->getFieldVTK(i)) {
#ifdef MERCURY_USE_MPI
                    if (particleMustBeWritten(p))
                    {
                        file << '\t' << f << '\n';
                    }
#else
                    file << '\t' << f << '\n';
#endif
                }
            }
            file << "  </DataArray>\n";
        }
    }
}
