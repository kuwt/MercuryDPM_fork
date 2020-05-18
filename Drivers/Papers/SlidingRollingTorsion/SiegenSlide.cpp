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

#include <iomanip>
#include "Siegen.h"

//We let a particle slide along a plane with a constant tangential velocity 
//and no rotation 
class Slide : public Siegen{
public:

	Slide(double NormalForce) : Siegen( NormalForce) {
		// set problem name
		setName("slide");

		// fixed tangential velocity of the particle
		TangentialVelocity = 1e-6*10; //m/s
		
		// set size of loop
		double RadiusPrime = particleHandler.getObject(0)->getRadius()/(1+relOverlap);
		RelLoopSize = 1e-6/RadiusPrime; //wrt Radius
		LoopTime = (4*RelLoopSize*RadiusPrime)/TangentialVelocity;
		std::cout << "LoopTime=" << LoopTime << std::endl;

		//time stepping
		setTimeMax(1.15*LoopTime);
		setSaveCount(helpers::getSaveCountFromNumberOfSavesAndTimeMaxAndTimeStep(4000,getTimeMax(),getTimeStep()));
		//eneFile.setSaveCount(1000000000);
		//dataFile.setSaveCount(1000000000);
		write(std::cout,true);
		
		//set wall
		InfiniteWall w0;
		w0.set(Vec3D(0.0,-1.0,0.0), Vec3D(0,0,0));
		wallHandler.copyAndAddObject(w0);
	
		//set_Particle
		particleHandler.getObject(0)->setPosition(Vec3D(0,RadiusPrime,0));
		particleHandler.getObject(0)->setVelocity(Vec3D(TangentialVelocity,0,0));

		//setNormal_force
		Mdouble Mass = particleHandler.getObject(0)->getMass();
		setGravity(Vec3D(0,-NormalForce/Mass,0));
	}
	
	void actionsBeforeTimeLoop() {
		//don't allow rotations
		particleHandler.getObject(0)->setInfiniteInertia();
	}
		
	void actionsBeforeTimeStep(){
		//fix tangential velocity
		Mdouble T1 = getTime()/LoopTime-floor(getTime()/LoopTime);
		if (T1>0.5) {
			Vec3D V = particleHandler.getObject(0)->getVelocity();
			Vec3D P = particleHandler.getObject(0)->getPosition();
			particleHandler.getObject(0)->setVelocity(Vec3D(-TangentialVelocity,V.Y,V.Z));
			particleHandler.getObject(0)->setPosition(Vec3D( TangentialVelocity*(1.-T1)*LoopTime,P.Y,P.Z));
		} else {
			Vec3D V = particleHandler.getObject(0)->getVelocity();
			Vec3D P = particleHandler.getObject(0)->getPosition();
			particleHandler.getObject(0)->setVelocity(Vec3D(-TangentialVelocity,V.Y,V.Z));
			particleHandler.getObject(0)->setPosition(Vec3D( TangentialVelocity*T1*LoopTime,P.Y,P.Z));
		}
	}

	void writeEneTimeStep(std::ostream& os) const{
		//MD::writeToEne();
		os
		<< " " << std::setw(12) << particleHandler.getObject(0)->getForce().X
		<< " " << std::setw(12) << particleHandler.getObject(0)->getPosition().X 
		<< std::endl;
	}

	void create_rough_wall(double Radius){
		setName("sliderough");
		//reset spring constants
		double NOverlap = particleHandler.getObject(0)->getRadius() 
			*sqrt(mathsFunc::square(1+relOverlap)-1)/Radius;
		std::cout << "N" << NOverlap << std::endl;
		//~ Mdouble Mass = particleHandler.getObject(0)->getMass();
		//~ setStiffnessAndRestitutionCoefficient(getStiffness()/NOverlap, 0.88, Mass);
        species->setStiffness(species->getStiffness() / NOverlap);
        species->setDissipation(species->getDissipation()/sqrt(NOverlap));
        species->setSlidingStiffness(species->getSlidingStiffness()/NOverlap);
        species->setSlidingDissipation(species->getSlidingDissipation()/sqrt(NOverlap));

		wallHandler.clear();
		double HalfWidth = particleHandler.getObject(0)->getRadius() 
			*(RelLoopSize + sqrt(mathsFunc::square(1+relOverlap)-1));
		//insert particle and get mass
		SphericalParticle p0;
		p0.setRadius(Radius);
		p0.fixParticle();
		double PositionX = -HalfWidth;
		while (PositionX<HalfWidth) {
			double PositionZ = -sqrt(mathsFunc::square(HalfWidth)+mathsFunc::square(PositionX));
			while (PositionZ<sqrt(mathsFunc::square(HalfWidth)+mathsFunc::square(PositionX))) {
				p0.setPosition(Vec3D(PositionX,0,PositionZ)+0.5*Radius*Vec3D(random.getRandomNumber(-1,1),random.getRandomNumber(-1,1),random.getRandomNumber(-1,1)));
				particleHandler.copyAndAddObject(p0);
				//std::cout << PositionX << "," << PositionY << std::endl;
				PositionZ += 2.0*Radius;
			}
			PositionX += 2.0*Radius;
		}
		std::cout << "N" << particleHandler.getNumberOfObjects() << std::endl;
	}

};

int main(int argc UNUSED, char *argv[] UNUSED)
{
	///Start off my solving the default problem
	Slide md(100e-6);
	//md.create_rough_wall(.4e-7);
	//md.speciesHandler.getObject()[0].setSlidingFrictionCoefficientStatic(.33);
	md.solve(argc, argv);
	std::cout << " " << std::setw(12) << md.Angle*180./constants::pi
		<< " " << std::setw(12) << md.species->getRollingFrictionCoefficient()
		<< " " << std::setw(12) << md.species->getTorsionFrictionCoefficient()
		<< " " << std::setw(12) <<-md.wallHandler.getObject(0)->getForce().X/md.wallHandler.getObject(0)->getForce().Y
		<< " " << std::setw(12) << md.wallHandler.getObject(0)->getForce().Y
		<< " " << std::setw(12) << md.particleHandler.getObject(0)->getVelocity().X;
	std::cout << std::endl;

	Slide md2(100e-6);
    //species->setForceType(ForceType::HERTZ_MINDLIN_DERESIEWICZ);
    //md2.species.setForceType(ForceType::HERTZ_MINDLIN_DERESIEWICZ);
	md2.setName("slideHMD");	
	md2.solve(argc, argv);

}
