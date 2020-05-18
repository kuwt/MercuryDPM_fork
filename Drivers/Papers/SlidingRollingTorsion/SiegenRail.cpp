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

///We let a particle slide between two planes with a constant tangential 
///velocity of the upper plane and no rotation 
class Slide : public Siegen{
public:

	Slide(Mdouble AngleDeg=0) : Siegen() {
		// set problem name
		setName("rail");
		// fixed tangential velocity of the particle
		TangentialVelocity = 1e-6*100; //m/s
		Angle = AngleDeg/180.*constants::pi;
		
		// set size of loop
		double RadiusPrime = particleHandler.getObject(0)->getRadius()/(1+relOverlap);
		RelLoopSize = .01; //wrt Radius
		LoopTime = (4*RelLoopSize*RadiusPrime)/TangentialVelocity;

		//time stepping
		//setTimeMax(1.3*LoopTime);
		//setSaveCount(helpers::getSaveCountFromNumberOfSavesAndTimeMaxAndTimeStep(2000,getTimeMax(),getTimeStep()));
		setTimeMax(0.24*LoopTime);
        setSaveCount(helpers::getSaveCountFromNumberOfSavesAndTimeMaxAndTimeStep(8, getTimeMax(), getTimeStep()));
		//eneFile.setSaveCount(1000000000);
		//dataFile.setSaveCount(1000000000);
		//write(std::cout,true);
		
		//set wall
		set_Walls(Angle);

		//setNormal_force
		setGravity(Vec3D(0,0,0));
		
	}
	
	void set_Walls(Mdouble Angle) {
		std::cout << "Walls " << std::endl;
		Mdouble c = cos(Angle);
		Mdouble s = sin(Angle);
		Mdouble Radius = particleHandler.getObject(0)->getRadius()/(1+relOverlap);
		std::cout << "Angle " << Angle*180/constants::pi;
		std::cout << "Radius " << particleHandler.getObject(0)->getRadius();
		std::cout << "RadiusPrime " << Radius << std::endl;
		setYMin(0);
		setYMax(Radius*(1.+1./c));
		setZMin(-2.*Radius);
		setZMax(+2.*Radius);
		if (Angle==0.0) {
            InfiniteWall w;
			w.set(Vec3D(0., 1., 0.), getMax());
            wallHandler.copyAndAddObject(w);
			w.set(Vec3D(0.,-1., 0.), Vec3D(0,0,0));
            wallHandler.copyAndAddObject(w);
		} else {
            InfiniteWall w;
			w.set(Vec3D(0., 1., 0.), getMax());
            wallHandler.copyAndAddObject(w);
			w.set(Vec3D(0., -c, -s), Vec3D(0,0,0));
            wallHandler.copyAndAddObject(w);
			w.set(Vec3D(0., -c,  s), Vec3D(0,0,0));
            wallHandler.copyAndAddObject(w);
		}
		particleHandler.getObject(0)->setPosition(Vec3D(0,Radius/c,0));
		particleHandler.getObject(0)->setVelocity(Vec3D(0,0,0));
	}

	void actionsBeforeTimeStep(){
        static Mdouble mus, mur, mut;
		if (getTime()==0) {
       		mus= species->getSlidingFrictionCoefficient();
        	mur= species->getRollingFrictionCoefficient();
            mut= species->getTorsionFrictionCoefficient();
            species->setSlidingFrictionCoefficient(0);
            species->setRollingFrictionCoefficient(0);
            species->setTorsionFrictionCoefficient(0);
        } else if (getTime()<0.5*getTimeMax()) {
            wallHandler.getObject(0)->setPosition(Vec3D(0.0,5e-6*(wallHandler.getObject(0)->getForce().Y-NormalForce),0.0));
            ///todo{DK: this is the old version of moving walls, the new version has to be tested}
			//wallHandler.getObject(0)->move(5e-6*(wallHandler.getObject(0)->getForce().Y-NormalForce));
		} else if (getTime()-getTimeStep()<0.5*getTimeMax()) {
//            TangentialSpringParticle* tsp0 = static_cast<TangentialSpringParticle*>(particleHandler.getObject(0));
//            std::cout << tsp0->get_TangentialSprings() << std::endl;
//            tsp0->get_TangentialSprings().reset();
//			std::cout << tsp0->get_TangentialSprings() << std::endl;
            wallHandler.getObject(0)->setVelocity(Vec3D(TangentialVelocity,0,0));
            ///todo{DK: this is the old version of moving walls, the new version has to be tested}
            //static_cast<InfiniteWall*>(wallHandler.getObject(0))->move(Vec3D(TangentialVelocity,0,0),0);
       		species->setSlidingFrictionCoefficient(mus);
            species->setRollingFrictionCoefficient(mur);
            species->setTorsionFrictionCoefficient(mut);

		}
// 		
// 		Mdouble T1 = getTime()/LoopTime-floor(getTime()/LoopTime);
// 		Mdouble T0 = (getTime()-getTimeStep())/LoopTime-floor(getTime()/LoopTime);
// 		if ((T1>0.25 && T0<0.25) ||
// 			(T1>0.75 && T0<0.75)) {
// 			wallHandler.getObject(0)->setVelocity(-wallHandler.getObject(0)->getVelocity());
// 		}
	}

	void writeEneTimeStep(std::ostream& os) const{
		//MD::writeToEne();
		os
		<< " " << std::setw(12) << wallHandler.getObject(0)->getForce().X
		<< " " << std::setw(12) << particleHandler.getObject(0)->getPosition().X
		<< " " << std::setw(12) << wallHandler.getObject(0)->getForce().Y
		<< " " << std::setw(12) << particleHandler.getObject(0)->getPosition().Y
		<< " " << std::setw(12) << Angle*180./constants::pi
		<< std::endl;
	}
	
};

int main(int argc UNUSED, char *argv[] UNUSED)
{
	FILE *counter_file=fopen("record", "w");
    char buffer [50];
	for (Mdouble a=0; a<71; a+=8) {
 		Slide md(a);
        sprintf (buffer, "rail%fR", a); md.setName(buffer);
        md.species->setTorsionFrictionCoefficient(0);
 		md.solve(argc, argv);
		fprintf(counter_file, "%f %.10e %.10e\n",a,md.wallHandler.getObject(0)->getForce().X,md.wallHandler.getObject(0)->getForce().Y);
 	}
	for (Mdouble a=0; a<71; a+=8) {
		Slide md(a);
        sprintf (buffer, "rail%f", a); md.setName(buffer);
        std::cout << buffer << std::endl;
		md.solve(argc, argv);
		fprintf(counter_file, "%f %.10e %.10e\n",a,md.wallHandler.getObject(0)->getForce().X,md.wallHandler.getObject(0)->getForce().Y);
	}
	fclose(counter_file);
}
