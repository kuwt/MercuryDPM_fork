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

#include "Mercury2D.h"
#include "Walls/InfiniteWall.h"
#include <iostream>
#include "Species/Species.h"
#include "Species/LinearViscoelasticSpecies.h"

/// In this file two particles are symmetrically placed in a bi-axial box are allowed to jump around under gravity. It tests walls gravity and symmetry.

class TwoParticleElasticCollisionInteraction : public Mercury2D {

	void setupInitialConditions() override {
		setMax({0.01,0.01,0.0});
		setGravity({0.0,-9.8,0.0});

		SphericalParticle P0,P1;
		P0.setSpecies(speciesHandler.getObject(0));
		P1.setSpecies(speciesHandler.getObject(0));
		P0.setPosition(Vec3D(0.006,0.005,0.0));
		P1.setPosition(Vec3D(0.004,0.005,0.0));
	
		P0.setVelocity(Vec3D(-0.1,0.0,0.0));
		P1.setVelocity(Vec3D( 0.1,0.0,0.0));
	
		P0.setRadius(0.0005);
		P1.setRadius(0.0005);
		particleHandler.copyAndAddObject(P0);
		particleHandler.copyAndAddObject(P1);
		
		wallHandler.clear();
		InfiniteWall w0;
		w0.setSpecies(speciesHandler.getObject(0));
		w0.set(Vec3D(-1, 0, 0), getMin());
		wallHandler.copyAndAddObject(w0);
		w0.set(Vec3D( 1, 0, 0), getMax());
		wallHandler.copyAndAddObject(w0);
		w0.set(Vec3D( 0,-1, 0), getMin());
		wallHandler.copyAndAddObject(w0);
		w0.set(Vec3D( 0, 1, 0), getMax());
		wallHandler.copyAndAddObject(w0);		
	}

};

int main(int argc UNUSED, char *argv[] UNUSED)
{
	TwoParticleElasticCollisionInteraction twoParticleElasticCollisionInteractionProblem;
    twoParticleElasticCollisionInteractionProblem.setName("TwoParticleElasticCollisionInteractionSelfTest");

    auto species = new LinearViscoelasticSpecies;
    twoParticleElasticCollisionInteractionProblem.speciesHandler.addObject(species);
    species->setDensity(2000);
    species->setStiffness(1e4);

	twoParticleElasticCollisionInteractionProblem.setTimeMax(0.25);
	twoParticleElasticCollisionInteractionProblem.setSaveCount(10);
    twoParticleElasticCollisionInteractionProblem.setTimeStep(2e-5);
	twoParticleElasticCollisionInteractionProblem.fStatFile.setFileType(FileType::NO_FILE);
	twoParticleElasticCollisionInteractionProblem.getInteractionFile().setFileType(FileType::ONE_FILE);

	twoParticleElasticCollisionInteractionProblem.solve();

}
