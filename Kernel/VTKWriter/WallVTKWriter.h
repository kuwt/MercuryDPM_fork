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


#ifndef WALL_VTKWRITER_H
#define WALL_VTKWRITER_H

#include "VTKWriter/BaseVTKWriter.h"
#include "WallHandler.h"

class WallVTKWriter final : public BaseVTKWriter<WallHandler>
{

public:
    
    /**
     * Non-default constructor; sets the handler and fileCounter
     */
    explicit WallVTKWriter(WallHandler& wallHandler) : BaseVTKWriter(wallHandler)
    {}
    
    /**
     * Default copy constructor
     */
    WallVTKWriter(const WallVTKWriter&) = default;
    
    /**
     * extracts vtk data from the wallHandler and stores it in a VTKContainer
     */
    void getVTKData(VTKContainer& vtk) const;
    
    /**
     * writes a vtk file
     */
    void writeVTK() const override;
    
    /**
     * the name of the class in the restart file
     */
    std::string getName() const
    { return "WallVTKWriter"; }

    void setWriteWallSurfaceAreaVTK(bool writeWallSurfaceAreaVTK);
    bool getWriteWallSurfaceAreaVTK() const;

protected:
    
    void write(std::fstream& file, std::string name, std::function<double(BaseWall*)> f) const;
    
    /**
     * writes the point data to the vtu file (i.e. the vertices of the mesh displayed in paraview)
     */
    void writeVTKPoints(std::fstream& file, VTKContainer& vtk) const;
    
    /**
     * writes the cell data to the vtu file (i.e. the faces of the mesh displayed in paraview)
     */
    void writeVTKCells(std::fstream& file, VTKContainer& vtk) const;
    
    /**
     * writes the cell data to the vtu file (i.e. surface area)
     */
    void writeVTKCellData(std::fstream& file, VTKContainer& vtk) const;

    /*!
     * \brief Calculates and writes the surface areas of the cells to the vtu file.
     */
    void writeVTKSurfaceArea(std::fstream& file, VTKContainer& vtk) const;

    bool writeWallSurfaceAreaVTK_ { false };
};


#endif
