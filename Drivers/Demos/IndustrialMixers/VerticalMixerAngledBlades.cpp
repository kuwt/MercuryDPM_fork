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

#include "VerticalMixer.h"

int main(int argc, char *argv[])
{
    std::string name = helpers::readFromCommandLine(argc, argv, "-name", std::string("VerticalMixerAngledBlades"));

    //first, do the real run
    VerticalMixerAngledBlades mixer(argc, argv);
    mixer.setName(name);
    mixer.removeOldFiles();
    // use straight blades
    mixer.bladeAngle_ = helpers::readFromCommandLine(argc, argv, "-bladeAngle", 0.0)*constants::pi;
    mixer.particleNumber_ = helpers::readFromCommandLine(argc, argv, "-particleNumber", 5000);
    mixer.setTimeMax(helpers::readFromCommandLine(argc, argv, "-timeMax", 30.0));
    mixer.setParticlesWriteVTK(true);
    mixer.solve();

    //then create pretty wall data
    VerticalMixerAngledBlades walls(argc, argv);
    walls.setName(name + "PrettyWalls");
    walls.removeOldFiles();
    walls.bladeAngle_ = mixer.bladeAngle_;
    walls.prettyWalls_ = true;
    walls.setTimeMax(mixer.getTimeMax());
    walls.solve();

    //finally, create pretty blade data
    VerticalMixerAngledBlades blades(argc, argv);
    blades.setName(name + "PrettyBlades");
    blades.removeOldFiles();
    blades.bladeAngle_ = mixer.bladeAngle_;
    blades.prettyWalls_ = true;
    blades.haveOuterWalls = false;
    blades.setTimeMax(mixer.getTimeMax());
    blades.solve();
    return 0;
}
