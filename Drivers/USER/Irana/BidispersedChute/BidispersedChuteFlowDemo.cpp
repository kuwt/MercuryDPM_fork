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

#include <StatisticsVector.h>
#include <CG/TimeAveragedCG.h>
#include "CG/CG.h"
#include "Logger.h"
#include "BidispersedChute.h"


int main(int argc, char* argv[])
{
    auto parameters = BidispersedChuteParameters(15, 22, 0.5);
    //parameters.setLargeParticleRadius(1.0 /3.0 * (2 + std::cbrt(2) - std::cbrt(4)));
    parameters.setLargeParticleRadius(0.75);
    BidispersedChute problem(parameters);
    problem.setChuteLength(20);
    problem.setChuteWidth(10);
    problem.setTimeMax(500);
    problem.restartFile.setFileType(FileType::MULTIPLE_FILES);
    problem.setRoughBottomType(RoughBottomType::MULTILAYER);
    problem.setName("BidispersedChuteH15A22Phi05R075");
    problem.setSaveCount(10000);
    /*auto cg0 = problem.cgHandler.copyAndAddObject(TimeAveragedCG<CGCoordinates::Z,CGFunctions::Lucy>());
    cg0->statFile.setSaveCount(10);
    cg0->statFile.setName(problem.getName() + ".LucyZ.stat");
    cg0->setNZ(200);*/
    problem.solve();
    return 0;
}
