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

#include <sstream>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <sys/types.h>
#include <sys/stat.h>
#include "Walls/InfiniteWall.h"
#include "Chute.h"
#include "Species/LinearViscoelasticFrictionSpecies.h"

//nohup nice -19 sc/quick_run BegerSetup > BegerSetup.out &

//used to simulate Dirk Begers LawinenBox
class LawinenBox : public Chute {
public:
	
	LawinenBox() {
        species  = speciesHandler.copyAndAddObject(LinearViscoelasticFrictionSpecies());
		random.randomise();
		//load Beger's bottom
		readDataFile("../ini/Bottom.data");
		
		//setXMax(getXMax()/1.);
		setYMax(getYMax()/10.);
		for (unsigned int i=0;i<particleHandler.getNumberOfObjects();) {
			if (particleHandler.getObject(i)->getPosition().X>getXMax()||
			    particleHandler.getObject(i)->getPosition().X<getXMin()||
			    particleHandler.getObject(i)->getPosition().Y>getYMax()||
			    particleHandler.getObject(i)->getPosition().Y<getYMin())
			{
				particleHandler.removeObject(i);
			}	else i++;
		}

		setName("Beger");
		setXBallsAdditionalArguments("-v0 -solidf");
		setFixedParticleRadius(particleHandler.getObject(0)->getRadius());
		setInflowParticleRadius(particleHandler.getObject(0)->getRadius());
		species->setDensity(2500);
		setGravity(Vec3D(0,0,-9.8));
		//collision time 5 ms
		Mdouble tc = 50e-4;
		Mdouble eps = 0.97; //0.97;
		Mdouble beta = 0.44; //0.44;
		Mdouble mass = species->getDensity()*mathsFunc::cubic(10e-3)*constants::pi/6.;
		species->setCollisionTimeAndNormalAndTangentialRestitutionCoefficient(tc, eps, beta, mass);
		setTimeStep(tc/50.);
		species->setSlidingFrictionCoefficient(0.1);

		species->setRollingStiffness(2./5.*species->getStiffness());
		species->setRollingDissipation(2./5.*species->getDissipation());
		species->setRollingFrictionCoefficient(0.);//getSlidingFrictionCoefficient()*0.2);
		//setSlidingFrictionCoefficient(getSlidingFrictionCoefficient()-getRollingFrictionCoefficient());

		for (unsigned int i=0; i<particleHandler.getNumberOfObjects(); i++) particleHandler.getObject(i)->fixParticle();
	}

	void actionsBeforeTimeStep(){
		Mdouble DeltaAngle = 0.025; //in degree
		Mdouble CheckInterval = 0.1; //in seconds
		static Mdouble NextCheck = CheckInterval;
		static Mdouble NextIncrease = 0;
		static Mdouble ChuteAngle = 0;
		
		//move chute while you are in the increase interval
		if (getTime()<NextIncrease) {
			//http://www.comp-fu.com/2012/01/nukes-smooth-ramp-functions/
			// smooth: traditional smoothstep
			double x=(getTime()-(NextIncrease-CheckInterval))/CheckInterval;
			setChuteAngle(ChuteAngle+x*x*(3 - 2*x)*DeltaAngle);
			if (getChuteAngleDegrees()>80) setTimeMax(getTime());
			//~ cout
			//~ << " t=" << x << getTime()
			//~ << " theta=" << getChuteAngleDegrees()
			//~ << endl;
		}	
		
		//move chute while you are in the check interval
		if (getTime()>NextCheck) {
			NextCheck = CheckInterval + getTime();
		
			Mdouble ene_kin = 0;
			for (unsigned int i=0;i<particleHandler.getNumberOfObjects();i++) if (!particleHandler.getObject(i)->isFixed())
				ene_kin += .5 * particleHandler.getObject(i)->getMass() * particleHandler.getObject(i)->getVelocity().getLengthSquared();
			if (ene_kin<1e-12) {
				NextIncrease = CheckInterval + getTime();
				ChuteAngle = getChuteAngleDegrees();
			}
			std::cout
			<< " t=" << getTime()
			<< " theta=" << getChuteAngleDegrees()
			<< " Ene_kin=" << ene_kin  
			<< std::endl;
		}
	}

	void writeEneTimeStep(std::ostream& os) const
	{
		Mdouble ene_kin = 0, ene_rot = 0, ene_gra = 0, mass_sum= 0, x_masslength=0, y_masslength=0, z_masslength=0;

		for (unsigned int i=0;i<particleHandler.getNumberOfObjects();i++) if (!particleHandler.getObject(i)->isFixed())
		{
			ene_kin += .5 * particleHandler.getObject(i)->getMass() * particleHandler.getObject(i)->getVelocity().getLengthSquared();
			ene_rot += particleHandler.getObject(i)->getRotationalEnergy();
			ene_gra -= particleHandler.getObject(i)->getMass() * Vec3D::dot(getGravity(),particleHandler.getObject(i)->getPosition());
			mass_sum +=particleHandler.getObject(i)->getMass();
			x_masslength +=particleHandler.getObject(i)->getMass()*particleHandler.getObject(i)->getPosition().X;
			y_masslength +=particleHandler.getObject(i)->getMass()*particleHandler.getObject(i)->getPosition().Y;
			z_masslength +=particleHandler.getObject(i)->getMass()*particleHandler.getObject(i)->getPosition().Z;
		} //end for loop over Particles

		///todo{Why is there a +6 here?}
		static int width = os.precision() + 6;
		os  << std::setw(width) << getTime() 
					   << " " << std::setw(width) << ene_gra
					   << " " << std::setw(width) << ene_kin
					   << " " << std::setw(width) << ene_rot
					   << " " << std::setw(width) << getElasticEnergy()
					   << " " << std::setw(width) << (mass_sum?x_masslength/mass_sum:NAN)
					   << " " << std::setw(width) << (mass_sum?y_masslength/mass_sum:NAN) 
					   << " " << std::setw(width) << (mass_sum?z_masslength/mass_sum:NAN)
					   << " " << std::setw(width) << getChuteAngleDegrees()
					   << std::endl;
		
		//resetElasticEnergy();
	}

	void create_inflow_particle()
	{
		P0.setRadius(random.getRandomNumber(getMinInflowParticleRadius(),getMaxInflowParticleRadius()));
		P0.setVelocity(Vec3D(0,0,0));
		Vec3D position;
		position.X = random.getRandomNumber(getXMin()+P0.getRadius(), getXMax()-P0.getRadius());
		position.Y = random.getRandomNumber(getYMin()+P0.getRadius(), getYMax()-P0.getRadius());
		position.Z = random.getRandomNumber(getZMin()+P0.getRadius(), getZMax()-P0.getRadius());
        P0.setPosition(position);
	}

	void printTime() const {
			//~ Mdouble ene_kin = 0;
			//~ for (unsigned int i=0;i<particleHandler.getNumberOfObjects();i++) if (!particleHandler.getObject(i)->isFixed())
				//~ ene_kin += .5 * particleHandler.getObject(i)->get_Mass() * particleHandler.getObject(i)->getVelocity().GetLength2();
			//~ cout << "t=" << setprecision(3) << left << setw(6) << getTime() 
			//~ << ", enekin=" << setprecision(3) << left << setw(6) << ene_kin
			//~ << endl;
	}

	void setupInitialConditions() {

	  boundaryHandler.clear();//WallsPeriodic.resize(0);
	  InfiniteWall w0;//Walls.resize(5);
		w0.set(Vec3D( 0.0, 0.0,-1.0), -getZMin());
                wallHandler.copyAndAddObject(w0);
		w0.set(Vec3D(-1.0, 0.0, 0.0), -getXMin());
                wallHandler.copyAndAddObject(w0);
		w0.set(Vec3D( 1.0, 0.0, 0.0),  getMax());
                wallHandler.copyAndAddObject(w0);
		w0.set(Vec3D( 0.0,-1.0, 0.0), getMin());
                wallHandler.copyAndAddObject(w0);
		w0.set(Vec3D( 0.0, 1.0, 0.0),  getMax());
                wallHandler.copyAndAddObject(w0);

		//number of flowing particles
		///
		int M=1000;
		int num_created = particleHandler.getNumberOfObjects();
                particleHandler.setStorageCapacity(particleHandler.getNumberOfObjects()+M);
		hGridActionsBeforeTimeLoop();
		hGridActionsBeforeTimeStep();
                	
		while (num_created<M){
			create_inflow_particle();
			if (checkParticleForInteraction(P0)) {
				num_created++;
			} else setZMax(getZMax()+0.00001);
		};
	}
  SphericalParticle P0;
    LinearViscoelasticFrictionSpecies* species;
};

int main(int argc, char *argv[])
{
	LawinenBox md;
	md.autoNumber();
	md.setTimeMax(1000);
	md.setSaveCount(1e4);
	md.eneFile.setSaveCount(1e2);
	md.restartFile.setFileType(FileType::ONE_FILE);//restartFile.setFileType(FileType::ONE_FILE);
	md.dataFile.setFileType(FileType::ONE_FILE);//dataFile.setFileType(FileType::ONE_FILE);
	md.solve(argc,argv);
	md.writeRestartFile();
}
