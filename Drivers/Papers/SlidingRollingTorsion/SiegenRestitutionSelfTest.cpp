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

#include "Siegen.h"
#include<fstream>
#include<iomanip>

///We let a particle slide along a plane with a constant tangential velocity 
///and no rotation 
class Slide : public Siegen {
public:

	Slide() : Siegen()
	{
		setName("SiegenRestitutionSelfTest");
		
		// set size of loop
		double Radius = particleHandler.getObject(0)->getRadius();
		LoopTime = sqrt(2.*100.*0.1*Radius/10.);
		std::cout << "LoopTime=" << LoopTime << std::endl;

		//time stepping
		setTimeMax(LoopTime);
		setSaveCount(helpers::getSaveCountFromNumberOfSavesAndTimeMaxAndTimeStep(200, getTimeMax(), getTimeStep()));

        //set wall
		InfiniteWall w;
        w.set(Vec3D(0., -1., 0.), Vec3D(0,0,0));
        wallHandler.copyAndAddObject(w);

		//set_Particle
		particleHandler.getObject(0)->setPosition(Vec3D(0,1.1*Radius,0));
		particleHandler.getObject(0)->setVelocity(Vec3D(0,0,0));

        setGravity(Vec3D(0,-10,0));

	}
};

int main(int argc UNUSED, char *argv[] UNUSED)
{
	Slide md;
	md.solve(argc, argv);
}
