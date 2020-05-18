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

#include "Chute.h"
#include "StatisticsVector.h"
#include <iostream>
#include <Species/LinearViscoelasticSlidingFrictionSpecies.h>

class Cstatic3D : public Chute {
public:
	
	void set_particle_properties() {
        auto species = speciesHandler.copyAndAddObject(LinearViscoelasticSlidingFrictionSpecies());
        //define species
	        species->setDensity(6./constants::pi);
		setGravity(Vec3D(sin(20.*constants::pi/180.),0,-cos(20.*constants::pi/180.)));
		setInflowParticleRadius(0.5);
		setFixedParticleRadius(0.5);
		//double mass = 4. / 3. * constants::pi * pow(getInflowParticleRadius(), 3.0) * S->getDensity();
                species->setCollisionTimeAndRestitutionCoefficient(0.005,0.88, species->getMassFromRadius(getInflowParticleRadius()));
		//~ setCollisionTimeAndRestitutionCoefficient(0.005,1);
		species->setSlidingStiffness(2.0/7.0* species->getStiffness());
		species->setSlidingDissipation(species->getDissipation());
		species->setSlidingFrictionCoefficient(0.5);
		setTimeStep(0.02 * species->getCollisionTime(species->getMassFromRadius(0.5 * (getMinInflowParticleRadius() + getMaxInflowParticleRadius()))));
		std::cout << getTimeStep() << std::endl;

		//		setDensity(6./pi);
		//setGravity(Vec3D(sin(20.*pi/180.),0,-cos(20.*pi/180.)));
		//setInflowParticleRadius(0.5);
		//setFixedParticleRadius(0.5);
		//setCollisionTimeAndRestitutionCoefficient(0.005,0.88);
		//~ setCollisionTimeAndRestitutionCoefficient(0.005,1);
		//setSlidingStiffness(2.0/7.0*getStiffness());
		//setSlidingDissipation(getDissipation());
		//setSlidingFrictionCoefficient(0.5);
		//setTimeStep();
		//cout << getTimeStep() << endl;
		//fix_hgrid();
	}
	
	void actionsBeforeTimeStep(){};
	
	void setupInitialConditions()
	{
		if (false) {
			// a single moving particle
			//set_NWall(0);
			setZMax(2);
			setYMax(1);
			setXMax(2);
			//set_N(2);
			SphericalParticle p0;
			p0.setPosition(Vec3D(1,0.5,0.5));
                        particleHandler.copyAndAddObject(p0);
			p0.setPosition(Vec3D(1,0.5,0.5)-Vec3D(sin(20.*constants::pi/180.),0,-cos(20.*constants::pi/180.)));
			particleHandler.copyAndAddObject(p0);
			particleHandler.getObject(0)->setVelocity(Vec3D(0.,0.,0.));
			particleHandler.getObject(1)->setVelocity(Vec3D(0.,0.,0.));
			particleHandler.getObject(0)->setRadius(0.5);
			particleHandler.getObject(1)->setRadius(0.5);
			particleHandler.getObject(0)->fixParticle();
		} else if (false) {
			// a two particles on two fixed ones
			//set_NWall(0);
			setZMax(2);
			setYMax(1);
			setXMax(2);
			//set_N(4);
			SphericalParticle p0;
		        p0.setPosition(Vec3D(0.5,0.5,0.5));
                        particleHandler.copyAndAddObject(p0);                       
		        p0.setPosition(Vec3D(1.5,0.5,0.5));
                        particleHandler.copyAndAddObject(p0);
 		        p0.setPosition(Vec3D(1,0.5,0.5+sqrt(0.75)));
                        particleHandler.copyAndAddObject(p0);
			p0.setPosition(Vec3D(1,0.5,0.5+sqrt(0.75))-Vec3D(sin(20.*constants::pi/180.),0,-cos(20.*constants::pi/180.)));
                        particleHandler.copyAndAddObject(p0);
			particleHandler.getObject(0)->setVelocity(Vec3D(0.,0.,0.));
                        particleHandler.getObject(1)->setVelocity(Vec3D(0.,0.,0.));
                        particleHandler.getObject(2)->setVelocity(Vec3D(0.,0.,0.));
                        particleHandler.getObject(3)->setVelocity(Vec3D(0.,0.,0.));
			//particleHandler.getObject(1)->getVelocity().set_zero();
			//Particles[2].Velocity.set_zero();
			//Particles[3].Velocity.set_zero();
			particleHandler.getObject(0)->setRadius(getInflowParticleRadius());
			particleHandler.getObject(1)->setRadius(getInflowParticleRadius());
                        particleHandler.getObject(2)->setRadius(getInflowParticleRadius());
                        particleHandler.getObject(3)->setRadius(getInflowParticleRadius());
			particleHandler.getObject(0)->fixParticle();
			particleHandler.getObject(1)->fixParticle();
		} else {
			readRestartFile("/storage/usr/people/weinhartt/DRIVERS/FlowRulePaper/run/full_runs/restart_flowrule/H10A24L1M0.5B0.5.restart"); 
			setName("static3d"); 
                   	dataFile.setFileType(FileType::ONE_FILE);
                       	restartFile.setFileType(FileType::ONE_FILE);
                	fStatFile.setFileType(FileType::ONE_FILE);
	                eneFile.setFileType(FileType::ONE_FILE);
                                   
			setTime(0.);
			setTimeMax(.5);
		}
		if (particleHandler.getNumberOfObjects()<10) write(std::cout,true);
		else write(std::cout,false);
	}

};

int main(int argc UNUSED, char *argv[] UNUSED)
{
	///Start off my solving the default problem
	Cstatic3D static3d;
	static3d.setName("static3d");
	static3d.set_particle_properties();
	static3d.setSaveCount(12000);
	static3d.setTimeMax(1);
	static3d.solve();
	//~ static3d.write(cout,true);
}
