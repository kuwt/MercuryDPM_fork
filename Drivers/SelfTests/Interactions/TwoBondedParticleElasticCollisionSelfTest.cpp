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

#include "Mercury3D.h"
#include "Walls/InfiniteWall.h"
#include <iostream>
#include "Species/Species.h"
#include "Species/LinearViscoelasticSlidingFrictionBondedSpecies.h"

/**
 * In this file two bonded particle pairs are placed in a box, and are allowed to jump around under gravity.
 */
class TwoBondedParticleElasticCollision : public Mercury3D {

	void setupInitialConditions() override {
		//set system dimension
		setMax(Vec3D(1,1,1)*0.01);
		setMin(Vec3D(0,0,0));

		//set particle properties
		SphericalParticle p;
        p.setSpecies(speciesHandler.getObject(0));
		p.setPosition(Vec3D(0.006,0.0051,0.0));
		p.setVelocity(Vec3D(-0.1,0.0,0.0));
		p.setRadius(0.0005);
		particleHandler.copyAndAddObject(p);
		p.setVelocity(Vec3D(-0.11,0.0,0.0));
		p.setPosition(Vec3D(0.006,0.006,0.0));
        particleHandler.copyAndAddObject(p);
        
		p.setPosition(Vec3D(0.004,0.005,0.0));
        p.setVelocity(Vec3D( 0.1,0.0,0.0));
        particleHandler.copyAndAddObject(p);
        p.setPosition(Vec3D(0.004,0.0059,0.0));
        particleHandler.copyAndAddObject(p);

		//bond certain particles
		auto i = dynamic_cast<BondedInteraction*>(interactionHandler.getInteraction(particleHandler.getObject(0),particleHandler.getObject(1),0));
		i->bond();
		i = dynamic_cast<BondedInteraction*>(interactionHandler.getInteraction(particleHandler.getObject(2),particleHandler.getObject(3),0));
		i->bond();

		//set wall properties
		wallHandler.clear();
		InfiniteWall w0;
        w0.setSpecies(speciesHandler.getObject(0));
		w0.set(Vec3D(-1, 0, 0), Vec3D(getXMin(), 0, 0));
		wallHandler.copyAndAddObject(w0);
		w0.set(Vec3D( 1, 0, 0), Vec3D(getXMax(), 0, 0));
		wallHandler.copyAndAddObject(w0);
		w0.set(Vec3D( 0,-1, 0), Vec3D(0, getYMin(), 0));
		wallHandler.copyAndAddObject(w0);
		w0.set(Vec3D( 0, 1, 0), Vec3D(0, getYMax(), 0));
		wallHandler.copyAndAddObject(w0);		
	}

};

int main(int argc UNUSED, char *argv[] UNUSED)
{
	//create an instance of the DPM class
	TwoBondedParticleElasticCollision twoParticleElasticCollisionProblem;
    twoParticleElasticCollisionProblem.setName("TwoBondedParticleElasticCollisionSelfTest");

	//set species properties
	auto species = new LinearViscoelasticSlidingFrictionBondedSpecies;
    twoParticleElasticCollisionProblem.speciesHandler.addObject(species);
    species->setDensity(2000);
	Mdouble radius = 0.0005;
	Mdouble mass = species->getMassFromRadius(radius);
	Mdouble tc = 1e-3;
	Mdouble restitution = 0.8;
    species->setCollisionTimeAndNormalAndTangentialRestitutionCoefficient(tc,restitution,restitution,mass);
    species->setSlidingFrictionCoefficient(0.05);
	//set bond strength such that particles have 0.5r overlap under equilibrium conditions (no external forces)
    species->setBondForceMax(species->getStiffness()*radius*0.5);

	//set global DPM properties
	twoParticleElasticCollisionProblem.setTimeMax(0.3);
	twoParticleElasticCollisionProblem.setSaveCount(20);
    twoParticleElasticCollisionProblem.setTimeStep(0.02*tc);
	twoParticleElasticCollisionProblem.fStatFile.setFileType(FileType::NO_FILE);
	twoParticleElasticCollisionProblem.setXBallsAdditionalArguments("-cmode 7 -v0 -solid -3dturn 1");
	twoParticleElasticCollisionProblem.solve();
}
