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

/* 
 * File:   mercury.cpp
 * Author: dducks
 *
 * Created on April 11, 2014, 1:30 PM
 */

#include <cstdlib>
#include <iostream>
#include <fstream>

#include <Mercury3D.h>
#include <Species/LinearViscoelasticSpecies.h>


#include "Serialization.h"
#include <cereal/archives/json.hpp>

using namespace std;



/*
 * 
 */
int main(int argc, char** argv) {
    
    std::ifstream file("../../../../sim.json");
    cereal::JSONInputArchive ar(file);
    
    SerializedProblem p;
    p.setSystemDimensions(3);
    p.setXMax(1.0);
    p.setYMax(1.0);
    p.setZMax(1.0);
    
    auto species = p.speciesHandler.copyAndAddObject(LinearViscoelasticSpecies());
    species->setDensity(2500.0); // sets the species type-0 density
    species->setStiffness(258.5);// sets the spring stiffness
    species->setDissipation(0.0);// sets the dissipation
    
    //p.setSaveCount(10);
    p.dataFile.setFileType(FileType::ONE_FILE);
    p.restartFile.setFileType(FileType::ONE_FILE);
    
    p.setXBallsAdditionalArguments("-solidf -v0");
    
    load(ar, p);
    
    p.solve(argc,argv);
    
    file.close();
    
    return 0;
}

