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



#include "Mercury2D.h"
#include "Walls/InfiniteWall.h"
#include <Species/LinearViscoelasticSpecies.h>

/// This case does a single elastic particle falling on an infinite plane.
/// The k is chosen so that the maximum overlap with the wall is around 2% of the particles diameter;
/// whereas, the time step must be taken to ensure 50 steps with a collision.
class FreeFallSelfTest : public Mercury2D
{
public:
	
	void setupInitialConditions() override
	{
		setMax({0.01,0.1,0.0});
		setGravity({0.0,-9.8,0.0});
		setTimeMax(1.0);

		InfiniteWall w0;
        w0.setSpecies(speciesHandler.getObject(0));
		w0.set(Vec3D(0,-1,0), Vec3D(0, getYMin(), 0));
		wallHandler.copyAndAddObject(w0);
		
		SphericalParticle p0;
        p0.setSpecies(speciesHandler.getObject(0));
		p0.setPosition(Vec3D(getXMax()/2,getYMax()*0.95,0.0));
		p0.setVelocity(Vec3D(0.0,0.0,0.0));
		p0.setRadius(0.005);
		particleHandler.copyAndAddObject(p0);
	}

};

int main(int argc, char* argv[])
{
	logger(INFO, "Single particle bouncing vertically on the bottom plate");
	// Make the problem and set the name
	FreeFallSelfTest freeFallSelfTestProblem;
	freeFallSelfTestProblem.setName("FreeFallSelfTest");

	//Set the species of the particle and wall, and its properties
	LinearViscoelasticSpecies species;
    species.setDensity(2000.0);
    freeFallSelfTestProblem.setParticleDimensions(3);
    species.setStiffness(8000000.0);
	freeFallSelfTestProblem.speciesHandler.copyAndAddObject(species);

	//set the parameters for the solver
	freeFallSelfTestProblem.setSaveCount(500);
	freeFallSelfTestProblem.fStatFile.setFileType(FileType::NO_FILE);
	freeFallSelfTestProblem.setWallsWriteVTK(FileType::ONE_FILE);
    freeFallSelfTestProblem.setTimeStep(1e-6);

	//solve the system, the single particle will now bounce on the plate
    freeFallSelfTestProblem.solve(argc, argv);
}
