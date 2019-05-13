//Copyright (c) 2013-2018, The MercuryDPM Developers Team. All rights reserved.
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

#include "Sinter.h"
#include <iostream>
#include "Walls/AxisymmetricIntersectionOfWalls.h"

enum IT {Flat, Spherical};
enum StageType {InitialRelaxation, Compression, RelaxationAfterCompression, MoveDownQuick, MoveDown, RelaxationAfterMoveDown, MoveUp};


class Indenter : public Sinter {
public:

	void setupInitialConditions() {	}
	
	Indenter(Mdouble IndenterRadius_, IT IndenterType_, Mdouble IndenterVelocity_, int N_, Mdouble IndentationDepth_, std::string restartfile_) 
	: IndenterRadius(IndenterRadius_), IndenterType(IndenterType_), IndenterVelocity(IndenterVelocity_), N(N_), IndentationDepth(IndentationDepth_), restartfile(restartfile_)
	{
		//Mdouble MinIndenterHeight = 0;

        //print input variables
        std::cout << "Creating Indenter with "
        << "IndenterRadius=" << IndenterRadius << "d, "
        << "IndenterType=" << IndenterType << " " 
        << "IndenterVelocity=" << IndenterVelocity << "sqrt(d/g) " 
        << std::endl;

        //check input variables
        if (IndenterRadius<0 || IndenterVelocity<0) {
            std::cerr << "error in constructor of Indenter; inputs have to be positive" << std::endl;
            exit(-1);
        }

        std::stringstream s, t;
        s << "Sinter" << restartfile << ".restart.0" << (N>=100?"":"0") << (N>=10?"":"0") << N;
       	t << "Indenter" << restartfile << "N" << (N>=100?"":"0") << (N>=10?"":"0") << N << "D" << IndentationDepth*1e6;
       	std::cout << "filenames: " << s.str() << " " << t.str() << std::endl;

	 	//setName("Sinter");
        readRestartFile(s.str().c_str());
		setRestarted(false);
		
		//reset parameters
		species->setDensity(2000);
		species->setSlidingFrictionCoefficient(0);
		Mdouble d = 4.5e-6;
		Mdouble m = 4.0/3.0*constants::pi*speciesHandler.getObject(0)->getDensity()*d*d*d;
		Mdouble tc = constants::pi/sqrt(1344/(m/2.0));
		Mdouble r = 0.88;
		//Mdouble g = 9.8e3;
	    Mdouble k2_k1_ratio = 5;
	    Mdouble kc_k1_ratio = 5;
	    Mdouble depth = 0.05; //delta0max=d/10=4.5e-7
		species->setCollisionTimeAndNormalAndTangentialRestitutionCoefficient(tc, r, r, m);
	    species->setPlasticParameters(species->getLoadingStiffness(), k2_k1_ratio* species->getLoadingStiffness(), kc_k1_ratio* species->getLoadingStiffness(), depth);
		//NeckGrowthExponent = 2.5e4;
		//Temperature = 4.0;
		NumSteps=10;
		std::cout << *species << std::endl;
		
		setName("Indenter");
		
		if (IndenterType==Flat) {
			setName("flatIndenter");
			Vec3D Axis(0, 0, 1);
			Vec3D PointOnAxis = Vec3D(0.5*(getXMax()+getXMin()), 0.5*(getYMax()+getYMin()), getBedHeight());
			Vec3D A(-1.4*IndenterRadius ,0.0, IndenterRadius);
			Vec3D B( 1.4*IndenterRadius ,0.0, IndenterRadius);
			Vec3D C( IndenterRadius, 0.0, 0.0);
			Vec3D D(-IndenterRadius, 0.0, 0.0);
			Vec3D E(0,1,0); //Periodic direction of the prism
			//FlatIndenter.addObject(Vec3D::Cross(A-B,E),A);
			FlatIndenter.addObject(Vec3D::cross(B-C,E),B);
			FlatIndenter.addObject(Vec3D::cross(C-D,E),C);
			FlatIndenter.addObject(Vec3D::cross(D-A,E),D);
            FlatIndenter.setPosition(PointOnAxis);
            FlatIndenter.setAxis(Axis);
	        wallHandler.addObject(&FlatIndenter);
		} else if (IndenterType==Spherical) {
			setName(t.str().c_str());
	        setFixedParticleRadius(IndenterRadius);
			SphericalIndenter.setPosition(Vec3D(0.5*(getXMax()+getXMin()), 0.5*(getYMax()+getYMin()), getBedHeight()+IndenterRadius));
			SphericalIndenter.fixParticle();
			//SphericalIndenter.set_IndSpecies(1);
			SphericalIndenter.setRadius(IndenterRadius);
			SphericalIndenter.setIndex(particleHandler.getNumberOfObjects());
	        //particleHandler.addObject(&SphericalIndenter);
		} else {
			exit(-1);
		}

		Mdouble dtDown = IndentationDepth/IndenterVelocity;
		setTimeMax(NumSteps*dtDown);
	 	setSaveCount(getTimeMax()/ getTimeStep()/200.0);
	 	dataFile.setFileType(FileType::ONE_FILE);
	 	restartFile.setFileType(FileType::MULTIPLE_FILES);
		setXBallsAdditionalArguments("-v0 -solidf");// -moh 250 -s 0.2 -oh 1300");
		write(std::cout,false);
	}

	void computeExternalForces(BaseParticle* CI)
	{
		if (!CI->isFixed()) {
			/// Now add on gravity
			CI->addForce(getGravity() * CI->getMass());
			///Finally walls
			//computeForcesDueToWalls(CI);
			//add spherical indenter as wall (quicker due to hGrid problems)
			if (IndenterType==Spherical)
				computeInternalForce(CI, &SphericalIndenter);
		}
	}

	void outputXBallsData(std::ostream& os) const
	{
		os  << particleHandler.getNumberOfObjects()+1 << " " << getTime() << " "
			<< getXMin() << " " << getYMin() << " " << getZMin() << " "
			<< getXMax() << " " << getYMax() << " " << getZMax() << " " << std::endl;
		// This outputs the particle data
		for (unsigned int i = 0;i<particleHandler.getNumberOfObjects();i++) 
			outputXBallsDataParticle(i,14,os);

		if (IndenterType==Flat) {
		    os
					<< FlatIndenter.getPosition().X << " " 
					<< FlatIndenter.getPosition().Y << " " 
					<< FlatIndenter.getPosition().Z+IndenterRadius << " " 
					<< FlatIndenter.getVelocity() << " "
					<< IndenterRadius << " 0 0 0 0 0 0 0"
					<< std::endl;
		} else if (IndenterType==Spherical) {
		    os
					<< SphericalIndenter.getPosition() << " " 
					<< SphericalIndenter.getVelocity() << " "
					<< SphericalIndenter.getRadius() << " "
					<< SphericalIndenter.getOrientation() << " "
					<< SphericalIndenter.getAngularVelocity() << " "
					<< getInfo(SphericalIndenter) <<std::endl;
		} else {
			exit(-1);
		}

	}

	void writeEneHeader(std::ostream& os) const
	{
		static int width = os.precision() + 6;
		os
		<< std::setw(width) << "t" << " "
		<< std::setw(width) << "ene_gra" << " "
		<< std::setw(width) << "ene_kin" << " "
		<< std::setw(width) << "ene_rot" << " "
		<< std::setw(width) << "ene_ela" << " "
		<< std::setw(width) << "X_COM" << " "
		<< std::setw(width) << "Y_COM" << " "
		<< std::setw(width) << "Z_COM" << " "
		<< std::setw(width) << "forceOnIndenter" << " "
		<< std::setw(width) << "positionOfIndenter" << " "
		<< std::endl;
	}

	void writeEneTimeStep(std::ostream& os) const
	{
		Mdouble ene_kin = 0, ene_rot = 0, ene_gra = 0, mass_sum= 0, x_masslength=0, y_masslength=0, z_masslength=0;

		for (std::vector<BaseParticle*>::const_iterator it = particleHandler.begin(); it!=particleHandler.end(); ++it) if (!(*it)->isFixed())
		{
			ene_kin += .5 * (*it)->getMass() * (*it)->getVelocity().getLengthSquared();
			ene_rot += particleHandler.getRotationalEnergy();
			ene_gra -= (*it)->getMass() * Vec3D::dot(getGravity(),(*it)->getPosition());
			mass_sum+= (*it)->getMass();
			x_masslength +=(*it)->getMass()*(*it)->getPosition().X;
			y_masslength +=(*it)->getMass()*(*it)->getPosition().Y;
			z_masslength +=(*it)->getMass()*(*it)->getPosition().Z;
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
			<< " " << std::setw(width) << getForceOnIndenter()
			<< " " << std::setw(width) << getIndenterHeight()-MinIndenterHeight
		<< std::endl;
	}

	Mdouble getIndenterHeight() const {
		if (IndenterType==Flat) {
			return FlatIndenter.getPosition().Z;
	    } else {
			return SphericalIndenter.getPosition().Z-SphericalIndenter.getRadius();
	    }
	}

	Mdouble getForceOnIndenter() const {
		if (IndenterType==Flat) {
			return FlatIndenter.getForce().Z ;
	    } else {
			return SphericalIndenter.getForce().Z;
	    }
	}

	void actionsAfterTimeStep(){
		if (getTime()<getTimeMax()*0.4*(1.0-3.5/NumSteps)) {
			
		} else if (getTime()+ getTimeStep()<getTimeMax()*(1.0-3.5/NumSteps)) {
            species->setNeckGrowthRate(0.0);
		} else if (getTime()<getTimeMax()*(1.0-3.5/NumSteps)) {
	 		std::cout<< "starting to indent" << std::endl;
	 		adjustIndenterHeight();
	 		MinIndenterHeight = getIndenterHeight();
		} else if (getTime()>getTimeMax()*(1.0-3.5/NumSteps)&& getTime()<getTimeMax()*(1.0-2.5/NumSteps)) {
			if (IndenterType==Flat) {
			    FlatIndenter.setVelocity(Vec3D(0.0,0.0,-IndenterVelocity));
				///todo{DK: This is the old moving wall functionallity, new one has to be tested}
				//FlatIndenter.move(Vec3D(0,0,-IndenterVelocity), getTimeStep());
			} else if (IndenterType==Spherical) {
			    SphericalIndenter.setVelocity(Vec3D(0,0,-IndenterVelocity));
				///todo{DK: This is the old moving wall functionallity, new one has to be tested}
				//SphericalIndenter.move(Vec3D(0,0,-IndenterVelocity* getTimeStep()));
			}
		} else if (getTime()>getTimeMax()*(1.0-1.5/NumSteps)) {
			if (IndenterType==Flat) {
				//FlatIndenter.move(Vec3D(0,0,IndenterVelocity), getTimeStep());
				///todo{DK: This is the old moving wall functionallity, new one has to be tested}
				FlatIndenter.setVelocity(Vec3D(0.0,0.0,IndenterVelocity));
			} else if (IndenterType==Spherical) {
			    SphericalIndenter.setVelocity(Vec3D(0,0,IndenterVelocity));
				///todo{DK: This is the old moving wall functionallity, new one has to be tested}
			    //SphericalIndenter.move(Vec3D(0,0,IndenterVelocity* getTimeStep()));
			}
		}
	}


	void actionsBeforeTimeStep(){
		if (IndenterType==Flat) {
			return FlatIndenter.setForce(Vec3D(0,0,0));
	    } else {
			return SphericalIndenter.setForce(Vec3D(0,0,0));
	    }
	}

	//Move the indenter down untill it touches a particle
	void adjustIndenterHeight() {
		for (Mdouble dz = 0.1e-6; dz>1e-11; dz/=4.0)
			while (!isTouchingBed(dz)) {};
	}


	bool isTouchingBed(Mdouble dz) {
		if (IndenterType==Flat) {
			FlatIndenter.move(Vec3D(0,0,-dz));
			Vec3D normal;
			Mdouble dist;
			for (std::vector<BaseParticle*>::const_iterator it=particleHandler.begin();it!=particleHandler.end();it++) {
	    		if (FlatIndenter.getDistanceAndNormal(*(*it), dist, normal)) {
					FlatIndenter.move(Vec3D(0,0,dz));
	    			return true;
	    		}
			}
		} else if (IndenterType==Spherical) {
			SphericalIndenter.move(Vec3D(0,0,-dz));
			for (std::vector<BaseParticle*>::const_iterator it=particleHandler.begin();it!=particleHandler.end();it++) {
	    		if (Vec3D::getDistanceSquared(SphericalIndenter.getPosition(),(*it)->getPosition())<=mathsFunc::square(SphericalIndenter.getRadius()+(*it)->getRadius())) {
					SphericalIndenter.move(Vec3D(0,0,dz));
	    			return true;
				}
			}
		} else {
			exit(-1);
		}
		return false;
	}

	Mdouble getBedHeight() {
	    Mdouble BedHeight = getZMin();
	    for (std::vector<BaseParticle*>::const_iterator it=particleHandler.begin();it!=particleHandler.end();it++) {
	        if (!(*it)->isFixed()) {
	        	Mdouble ParticleHeight = (*it)->getPosition().Z+(*it)->getRadius();
	        	if (ParticleHeight>BedHeight) BedHeight=ParticleHeight;
	        }
	    }
	    // std::cout << "BedHeight" << BedHeight << std::endl;
	    return BedHeight;
	}

private:
	//Member variables
	Mdouble IndenterRadius;
	IT IndenterType;
	Mdouble IndenterVelocity;
	AxisymmetricIntersectionOfWalls FlatIndenter;
	SphericalParticle SphericalIndenter;
	//StageType Stage;
	Mdouble MinIndenterHeight;
    int N;
    Mdouble IndentationDepth;
    std::string restartfile;
    Mdouble NumSteps;
};

int main(int argc UNUSED, char *argv[] UNUSED)
{
    int N;
    Mdouble IndentationDepth;
	std::cout << "V" << argc << std::endl;
	std::stringstream s("");
	if (argc<3) {
		s << "Sinter";
		N = 100;
	   	IndentationDepth = 0.1e-6;
	} else {
		s << argv[1];
		N = atoi(argv[2]);
	    IndentationDepth = atof(argv[3])*1e-6;
	    std::cout << s.str() << " N" << N << "ID" << IndentationDepth << std::endl;
	}
 
 	Mdouble IndenterRadius = 21e-6;
	Mdouble IndenterVelocity = 3e-3;
	std::cout << "spherical Indenter" << std::endl;	
	Indenter SI(IndenterRadius, Spherical, IndenterVelocity, N, IndentationDepth, s.str());

   SI.solve();
}
//to plot the force curve
//load '../../../Source/DRIVERS/Sintering/plotIndenter.gnu'
