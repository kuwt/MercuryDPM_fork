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

#include "Sinter2.h"
#include <iostream>

class two_particle_collision : public Sinter {
public:
	void setupInitialConditions() override {	}
	
	void add_particles() {}
	
	void actionsBeforeTimeStep() override {
		if (getTime()<0.7*getTimeMax())
            setTemperature(1.3);
		else
            setTemperature(0.0);
	}		
};

int main(int argc UNUSED, char *argv[] UNUSED)
{
	two_particle_collision md;
	md.setName("Bed");
    if (argc>=2) {
        std::stringstream s;
        s << "BedH" << atoi(argv[1]);
		md.setName(s.str().c_str());
	}
    md.readRestartFile();
 	md.setRestarted(false);
 	md.setName("Sinter");
    if (argc>=2) {
        std::stringstream s;
        s << "SinterH" << atoi(argv[1]);
		md.setName(s.str().c_str());
	}
 	md.NeckGrowthSpeed = 400e-6;	//a^2/2, delta = a sqrt(t)
 	md.setXBallsAdditionalArguments("-v0 -solidf");
 	md.setTimeMax(0.25*md.getTimeMax());
 	md.setSaveCount(md.getTimeMax()/ md.getTimeStep()/100);
 	md.restartFile.setFileType(FileType::MULTIPLE_FILES);
 	md.fStatFile.setFileType(FileType::ONE_FILE);
    md.solve();
}
//to see how much it melted:
//p 'Sinter.ene' u 1:($8*1e6)
