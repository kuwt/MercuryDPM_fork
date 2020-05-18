//Copyright (c) 2013-2020, The MercuryDPM Developers Team. All rights reserved.
//For the list of developers, see <http://MercuryDPM.org/Team>.
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

/** \brief Test of the STL reader. The files read is an STL file containing 12 triangles that forms a unit cube, created in autocad
**/

#include<fstream>
#include<iostream>
#include<vector>
#include <string>
#include<BinaryReader.h>
#include<Logger.h>
#include<Math/Vector.h>
#include<CMakeDefinitions.h>
#include <Mercury3D.h>
#include <Species/LinearViscoelasticSpecies.h>


int main()
{
    std::string directory = getMercurySourceDir();
    std::string filename = directory+"/Drivers/ImportTools/ExampleSTLFiles/Box1x1x1.stl";
    Mercury3D dpm;
    auto species = dpm.speciesHandler.copyAndAddObject(LinearViscoelasticSpecies());
    dpm.wallHandler.readTriangleWall(filename,species);
    helpers::check(dpm.wallHandler.getNumberOfObjects(),12,0,"Number of walls read");
    dpm.write(std::cout,true);
    return 0;
}
