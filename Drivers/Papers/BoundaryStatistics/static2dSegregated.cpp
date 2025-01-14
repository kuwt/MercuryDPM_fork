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

// This produces the picture for the Powders and Grains 2013 conference talk by Thomas Weinhart.
#include "Chute.h"
#include "StatisticsVector.h"
#include "Walls/InfiniteWall.h"
#include <iostream>
#include <Species/LinearViscoelasticSlidingFrictionSpecies.h>

class Cstatic2d : public Chute {
public:
	
	void set_particle_properties() {
		//define species
        auto species = speciesHandler.copyAndAddObject(LinearViscoelasticSlidingFrictionSpecies());
        species->setDensity(6. / constants::pi);
        setGravity(Vec3D(0, 0, -1));
        setInflowParticleRadius(0.5);
        setFixedParticleRadius(0.5);
        //double mass = 4. / 3. * constants::pi * pow(getInflowParticleRadius(), 3.0) * S->getDensity();
        species->setCollisionTimeAndRestitutionCoefficient(0.005, 0.88,
                                                           species->getMassFromRadius(getInflowParticleRadius()));
        //~ setCollisionTimeAndRestitutionCoefficient(0.005,1);
        species->setSlidingStiffness(2.0 / 7.0 * species->getStiffness());
        species->setSlidingDissipation(species->getDissipation());
        species->setSlidingFrictionCoefficient(0.5);
        setTimeStep(0.02 * species->getCollisionTime(
                species->getMassFromRadius(0.5 * (getMinInflowParticleRadius() + getMaxInflowParticleRadius()))));
        logger(INFO, "%", getTimeStep());
        
        //		setDensity(6./constants::pi);
        //setGravity(Vec3D(0,0,-1));
        //setInflowParticleRadius(0.5);
        //setFixedParticleRadius(0.5);
        //
        //setCollisionTimeAndRestitutionCoefficient(0.005,0.88);
        //~ setCollisionTimeAndRestitutionCoefficient(0.005,1);
        //setSlidingStiffness(2.0/7.0*getStiffness());
        //setSlidingDissipation(getDissipation());
		//setSlidingFrictionCoefficient(0.5);
		//setTimeStep();
		//std::cout << getTimeStep() << std::endl;
		//fix_hgrid();
	}
	
	void actionsBeforeTimeStep() override {};
	
	void setupInitialConditions() override
	{
		if (true) {
			//set domain and walls
			setZMax(5);
			setYMax(1);
			setXMax(5);
						
			InfiniteWall w0;
			w0.set(Vec3D( 0.0, 0,-1.0), getMin());
			wallHandler.copyAndAddObject(w0);
			w0.set(Vec3D( 0.0, 0, 1.0), getMax());
			wallHandler.copyAndAddObject(w0);
			//w0.set(Vec3D(-1.0, 0, 0.0), -getXMin());
			//get_WallHandler().copyAndAddWall(w0);
			//w0.set(Vec3D( 1.0, 0, 0.0),  getXMax());
			//get_WallHandler().copyAndAddWall(w0);
			
			PeriodicBoundary b0;
			b0.set(Vec3D( 1.0, 0.0, 0.0), getXMin(), getXMax());
			boundaryHandler.copyAndAddObject(b0);

			
			//set particles
			SphericalParticle p0;
			p0.setRadius(getFixedParticleRadius());
			p0.fixParticle();

			p0.setPosition(Vec3D(0.5,0.5,0.1));
			particleHandler.copyAndAddObject(p0);
			p0.setPosition(Vec3D(1.5,0.5,-.2));
			particleHandler.copyAndAddObject(p0);
			p0.setPosition(Vec3D(2.5,0.5,0));
			particleHandler.copyAndAddObject(p0);
			p0.setPosition(Vec3D(3.5,0.5,0.1));
			particleHandler.copyAndAddObject(p0);
			p0.setPosition(Vec3D(4.5,0.5,0.05));
			particleHandler.copyAndAddObject(p0);
			
			SphericalParticle p1;
			p1.setRadius(getInflowParticleRadius());
			
			p1.setPosition(Vec3D(1,0.5,1));
			p1.setVelocity(Vec3D(random.getRandomNumber(-.1,.1),0,random.getRandomNumber(-.1,.1)));
			particleHandler.copyAndAddObject(p1);
			p1.setPosition(Vec3D(2,0.5,1));
			p1.setVelocity(Vec3D(random.getRandomNumber(-.1,.1),0,random.getRandomNumber(-.1,.1)));
			particleHandler.copyAndAddObject(p1);
			p1.setPosition(Vec3D(3,0.5,1));
			p1.setVelocity(Vec3D(random.getRandomNumber(-.1,.1),0,random.getRandomNumber(-.1,.1)));
			particleHandler.copyAndAddObject(p1);
			p1.setPosition(Vec3D(4,0.5,1));
			p1.setVelocity(Vec3D(random.getRandomNumber(-.1,.1),0,random.getRandomNumber(-.1,.1)));
			particleHandler.copyAndAddObject(p1);
			p1.setPosition(Vec3D(1.5,0.5,2));
			p1.setVelocity(Vec3D(random.getRandomNumber(-.1,.1),0,random.getRandomNumber(-.1,.1)));
			particleHandler.copyAndAddObject(p1);

			p1.setPosition(Vec3D(1,0.5,2));
			p1.setVelocity(Vec3D(random.getRandomNumber(-.1,.1),0,random.getRandomNumber(-.1,.1)));
			particleHandler.copyAndAddObject(p1);
			p1.setPosition(Vec3D(2,0.5,2));
			p1.setVelocity(Vec3D(random.getRandomNumber(-.1,.1),0,random.getRandomNumber(-.1,.1)));
			particleHandler.copyAndAddObject(p1);
			p1.setPosition(Vec3D(3,0.5,2));
			p1.setVelocity(Vec3D(random.getRandomNumber(-.1,.1),0,random.getRandomNumber(-.1,.1)));
			particleHandler.copyAndAddObject(p1);
			p1.setPosition(Vec3D(4,0.5,2));
			p1.setVelocity(Vec3D(random.getRandomNumber(-.1,.1),0,random.getRandomNumber(-.1,.1)));
			particleHandler.copyAndAddObject(p1);
			p1.setPosition(Vec3D(1.5,0.5,3));
			p1.setVelocity(Vec3D(random.getRandomNumber(-.1,.1),0,random.getRandomNumber(-.1,.1)));
			particleHandler.copyAndAddObject(p1);
            
            p1.setRadius(getInflowParticleRadius() * 1.5);
            
            p1.setPosition(Vec3D(1, 0.5, 3.25));
            p1.setVelocity(Vec3D(random.getRandomNumber(-.1, .1), 0, random.getRandomNumber(-.1, .1)));
            particleHandler.copyAndAddObject(p1);
            p1.setPosition(Vec3D(2, 0.5, 3.25));
            p1.setVelocity(Vec3D(random.getRandomNumber(-.1, .1), 0, random.getRandomNumber(-.1, .1)));
            particleHandler.copyAndAddObject(p1);
            p1.setPosition(Vec3D(3, 0.5, 3.25));
            p1.setVelocity(Vec3D(random.getRandomNumber(-.1, .1), 0, random.getRandomNumber(-.1, .1)));
            particleHandler.copyAndAddObject(p1);
            p1.setPosition(Vec3D(4, 0.5, 3.25));
            p1.setVelocity(Vec3D(random.getRandomNumber(-.1, .1), 0, random.getRandomNumber(-.1, .1)));
            particleHandler.copyAndAddObject(p1);
            p1.setPosition(Vec3D(1.5, 0.5, 4));
            p1.setVelocity(Vec3D(random.getRandomNumber(-.1, .1), 0, random.getRandomNumber(-.1, .1)));
            particleHandler.copyAndAddObject(p1);
            
        }
        if (particleHandler.getNumberOfObjects() < 10)
        {
            write(std::cout, true);
        }
        else
        {
            write(std::cout, false);
        }
    }

};

int main(int argc UNUSED, char *argv[] UNUSED)
{
	///Start off my solving the default problem
	Cstatic2d static2d;
	static2d.setName("static2d");
	static2d.set_particle_properties();
	static2d.setSaveCount(12000);
	static2d.setTimeMax(400);
	static2d.solve();
	//~ static2d.write(cout,true);
}
