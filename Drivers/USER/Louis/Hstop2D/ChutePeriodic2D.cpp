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

#include <iostream>
#include "ChutePeriodic2D.h"

int main(int argc, char* argv[])
{

    // read arguments
    if (argc < 5)
    {
        std::cout << "Please specify Height, Angle, BottomType, SizeRatio, and Spacing." << std::endl;
        exit(-1);
    }
    Mdouble h=strtod(argv[1],nullptr);
    Mdouble a=strtod(argv[2],nullptr);
    auto bt= static_cast<int>(strtod(argv[3], nullptr));
    Mdouble srb=strtod(argv[4],nullptr);
    Mdouble eps=strtod(argv[5],nullptr);
    std::string btShort;
    BottomType bottomType;
    switch(bt) {
        case 0 : btShort = "F"; bottomType = FLAT; break;
        case 1 : btShort = "O"; bottomType = ORDERED; break;
        case 2 : default:  btShort = "D"; bottomType = DISORDERED; break;
    }

    // name
    std::stringstream name;
    name << "H" << h
         << "A" << a
         << "_" << btShort
         << "_R" << srb
         << "E" << eps;
    std::cout << "Case is named as: " << name.str().c_str() << std::endl;

    // solve problem
    ChutePeriodic2D problem;
    problem.setName(name.str().c_str());
    problem.setInflowHeight(h);
    problem.setChuteAngle(a);
    problem.setBottomType(bottomType);
    problem.setFixedParticleRadiusRatio(srb);
    problem.setFixedParticleSpacing(eps);
    problem.setTimeMax(3000);
    problem.cgHandler.getObject(0)->setY(problem.getYMin(),problem.getYMax());
    problem.cgHandler.getObject(1)->setY(problem.getYMin(),problem.getYMax());
    problem.cgHandler.getObject(1)->setTimeMin(2000);
    problem.solve();

    return 0;
}