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


///This code is a example on how to write a restartable mercury code
#include "DPMBase.h"
#include "Walls/InfiniteWall.h"
#include "Species/LinearViscoelasticSpecies.h"
class FreeFall : public DPMBase{
public:

    /// Set species, particles, boundaries, walls here
    void setupInitialConditions() override
	{
        LinearViscoelasticSpecies s;
        s.setStiffness(1e5);
        s.setDissipation(0.01);
        s.setDensity(2e3);
        auto sp = speciesHandler.copyAndAddObject(s);

        InfiniteWall w0;
        w0.setSpecies(speciesHandler.getObject(0));
        w0.set(Vec3D(0,0,-1), Vec3D(0, 0, 0));
		wallHandler.copyAndAddObject(w0);

		SphericalParticle p0;
        p0.setSpecies(speciesHandler.getObject(0));
		p0.setPosition(Vec3D(0,0,getZMax()-1e-3));
		p0.setRadius(1e-3);
		particleHandler.copyAndAddObject(p0);
    }
    
    void actionsOnRestart() override
    {
        //this function replaces setupInitialConditions in case of a restart
        //if you have problem-specific variables, set them here
    }

};

int main(int argc, char *argv[])
{
	///Start off my solving the default problem
	FreeFall dpm;

    //set FreeFall-specific parameters here
    dpm.setName("FreeFallRestart");
    dpm.setSaveCount(1000);
    dpm.setTimeStep(1e-6);
    dpm.setTimeMax(0.5);
    dpm.setDomain(Vec3D(-1e-3,-1e-3,0e-3),Vec3D(1e-3,1e-3,10e-3));
    dpm.setGravity(Vec3D(0,0,-9.8));

    //restart and run until final time 1, using command line arguments:
    // ./freeFallRestart -r -tmax 1
    dpm.solve(argc,argv);
    return 0;
}
