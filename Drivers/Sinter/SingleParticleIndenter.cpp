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

#include "SphericalIndenter.h"
#include "Species/SinterSpecies.h"
#include "Walls/InfiniteWall.h"
#include "Boundaries/PeriodicBoundary.h"
#include "Walls/AxisymmetricIntersectionOfWalls.h"
using std::cout;
using std::endl;

/// Single particle, indented slowly by spherical indenter.

class SingleParticleIndenter : public SphericalIndenter
{
public:

    SingleParticleIndenter(Mdouble indenterDiameter, Mdouble indentationVelocity, Mdouble indentationForce)
        : SphericalIndenter(indenterDiameter, indentationVelocity, indentationForce)
    {
        setName("SingleParticleConstantTemperature");
        readRestartFile();
        setRestarted(false);
        setTime(0);
        setName("SingleParticleIndenter");
        std::cout << "Restarting SingleParticleConstantTemperature N=" << particleHandler.getNumberOfObjects() << " \n";
        write(std::cout);
    }

};

int main(int argc UNUSED, char *argv[] UNUSED)
{
    Mdouble timeMax = 2e-5;
    Mdouble indenterDiameter = 127e-6;
    Mdouble indentationDepth = 1e-6; //4mN, 1um
    Mdouble indentationVelocity = indentationDepth/timeMax*2.0;
    Mdouble indentationForce = 4e-3;

    SingleParticleIndenter sp(indenterDiameter, indentationVelocity, indentationForce);
    sp.setFileType(FileType::ONE_FILE);
    sp.setXBallsAdditionalArguments(" -v0 -solidf ");
    sp.setSaveCount(100);
    sp.setTimeMax(2.0*timeMax);
    sp.solve();

    std::cout << "Execute 'gnuplot SingleParticleIndenter.gnu' to view output" << std::endl;
    helpers::writeToFile("SingleParticleIndenter.gnu",
                         "set xlabel 'displacement [um]'\n"
                         "set ylabel 'force [mN]'\n"
                         "p 'SingleParticleIndenter.ene' u (-$2*1e6):($3*1e3) w lp\n"
                         );
}
