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

#include "Species/Species.h"
#include "Species/NormalForceSpecies/LinearViscoelasticNormalSpecies.h"
#include "Species/NormalForceSpecies/LinearPlasticViscoelasticNormalSpecies.h"
#include "Species/FrictionForceSpecies/FrictionSpecies.h"
#include "Species/AdhesiveForceSpecies/ReversibleAdhesiveSpecies.h"
//#include "Species/AdhesiveForceSpecies/IrreversibleAdhesiveSpecies.h"
#include "DPMBase.h"
#include "Walls/InfiniteWall.h"
#include "Logger.h"


#include "Walls/AxisymmetricIntersectionOfWalls.h"
#include "Mercury3D.h"
#include <string>
#include <sstream>

///This code is written to test short-distance non-contact
///forces such as van-der-Waals or liquid bridge forces.
///A simple reversible adhesion force is combined with a linear-elastic
///contact force.
class CoDrum : public Mercury3D{
public:

    Species<LinearViscoelasticNormalSpecies, FrictionSpecies,ReversibleAdhesiveSpecies>* species1;
    Species<LinearViscoelasticNormalSpecies, FrictionSpecies,ReversibleAdhesiveSpecies>* species2;
    Species<LinearViscoelasticNormalSpecies, FrictionSpecies,ReversibleAdhesiveSpecies>* speciesWall;

    CoDrum()
    {
        //Creating species
        species1 = new Species<LinearViscoelasticNormalSpecies, FrictionSpecies,ReversibleAdhesiveSpecies>;
        speciesHandler.addObject(species1);
        species2 = new Species<LinearViscoelasticNormalSpecies, FrictionSpecies,ReversibleAdhesiveSpecies>;
        speciesHandler.addObject(species2);
        speciesWall = new Species<LinearViscoelasticNormalSpecies, FrictionSpecies,ReversibleAdhesiveSpecies>;
        speciesHandler.addObject(speciesWall);
    }

	void setupInitialConditions() override
	{
        setParticleDimensions(3);
        /////////////////////////////////////////////////////////////////////////////////////////////
        ////////////////////////////////Setting Species Properties///////////////////////////////////
        /////////////////////////////////////////////////////////////////////////////////////////////

		//calculating masses of all species and effective mass of particles 1 and 2 or use in
		//later calculations
		double mass1 = 4.0 / 3.0 * constants::pi * pow(radiusSpecies1, 3.0) * densitySpecies1;
		double mass2 = 4.0 / 3.0 * constants::pi * pow(radiusSpecies2, 3.0) * densitySpecies2;
		double mEff = 2*mass1*mass2 / (mass1 + mass2);
        //******************************SPECIES 1****************************************************
        species1->setDensity(densitySpecies1);
        //TO FIX - GET ACCURATE MASS!
        species1->setStiffnessAndRestitutionCoefficient(k,e11,mass1);
	    species1->setAdhesionStiffness(k_a11);
        species1->setAdhesionForceMax(f_a11*species1->getAdhesionStiffness());

        ///////////////////////////////FRICTIONAL PROPERTIES/////////////////////////////////////////
        species1->setSlidingStiffness(species1->getStiffness()*2./7.);
        species1->setSlidingFrictionCoefficient(mu11);
        species1->setSlidingDissipation(species1->getDissipation()*2./7.);

		species1->setRollingStiffness(species1->getStiffness()*2.0/7.0);
		species1->setRollingFrictionCoefficient(mur11);
		species1->setRollingDissipation(species1->getDissipation()*2./7.);


        //******************************SPECIES 2****************************************************
        species2->setDensity(densitySpecies2);
        //TO FIX - GET ACCURATE MASS!
        species2->setStiffnessAndRestitutionCoefficient(k,e22,mass2);
        species2->setAdhesionStiffness(k_a22);
		species2->setAdhesionForceMax(f_a22*species2->getAdhesionStiffness());

        ///////////////////////////////FRICTIONAL PROPERTIES/////////////////////////////////////////
        species2->setSlidingStiffness(species2->getStiffness()*2./7.);
        species2->setSlidingFrictionCoefficient(mu22);
        species2->setSlidingDissipation(species2->getDissipation()*2./7.);

		species2->setRollingStiffness(species2->getStiffness()*2.0/7.0);
		species2->setRollingFrictionCoefficient(mur22);
		species2->setRollingDissipation(species2->getDissipation()*2./7.);


        //******************************WALL SPECIES****************************************************
        speciesWall->setDensity(2000.0);
        //TO FIX - GET ACCURATE MASS!
        speciesWall->setStiffnessAndRestitutionCoefficient(k,0.01,mEff);

        ///////////////////////////////ADHESIVE PROPERTIES///////////////////////////////////////////
	    speciesWall->setAdhesionStiffness(k_aWW);
	    speciesWall->setAdhesionForceMax(f_aWW*speciesWall->getAdhesionStiffness());

        ///////////////////////////////FRICTIONAL PROPERTIES/////////////////////////////////////////
        speciesWall->setSlidingStiffness(speciesWall->getStiffness()*2./7.);
        speciesWall->setSlidingFrictionCoefficient(0.1);
        speciesWall->setSlidingDissipation(speciesWall->getDissipation()*2./7.);

		speciesWall->setRollingStiffness(speciesWall->getStiffness()*2.0/7.0);
		speciesWall->setRollingFrictionCoefficient(0.01);
		speciesWall->setRollingDissipation(speciesWall->getDissipation()*2./7.);

        /////////////////////////////////////////////////////////////////////////////////////////////
        ////////////////////////////////Setting Mixed Properties/////////////////////////////////////
        /////////////////////////////////////////////////////////////////////////////////////////////

        //******************************SPECIES 1 AND 2***********************************************
		auto species1And2 = speciesHandler.getMixedObject(species1,species2);
		//(Stiffness, CoR, effective mass)
		//Note: for final version, consider switching out to set tc instead of k
		species1And2->setStiffnessAndRestitutionCoefficient(k, e12, mEff);

        ///////////////////////////////ADHESIVE PROPERTIES///////////////////////////////////////////
	    species1And2->setAdhesionStiffness(k_a12);
	    species1And2->setAdhesionForceMax(f_a12*species1And2->getAdhesionStiffness());

        ///////////////////////////////FRICTIONAL PROPERTIES/////////////////////////////////////////
		species1And2->setSlidingDissipation(species1And2->getDissipation()*2./7.);
		species1And2->setSlidingFrictionCoefficient(mu12);
		species1And2->setSlidingStiffness(species1And2->getStiffness()*2.0/7.0);

		species1And2->setRollingStiffness(species1And2->getStiffness()*2.0/7.0);
		species1And2->setRollingFrictionCoefficient(mur12);
		species1And2->setRollingDissipation(species1And2->getDissipation()*2./7.);

        //******************************SPECIES 1 AND WALL***********************************************
		auto species1AndWall = speciesHandler.getMixedObject(species1,speciesWall);
		//(Stiffness, CoR, effective mass)
		//Note: for final version, consider switching out to set tc instead of k
		species1AndWall->setStiffnessAndRestitutionCoefficient(k, e1W, 2*mass1);

        ///////////////////////////////ADHESIVE PROPERTIES///////////////////////////////////////////
	    species1AndWall->setAdhesionStiffness(k_a1W);
	    species1AndWall->setAdhesionForceMax(f_a1W*species1AndWall->getAdhesionStiffness());

        ///////////////////////////////FRICTIONAL PROPERTIES/////////////////////////////////////////
		species1AndWall->setSlidingDissipation(species1AndWall->getDissipation()*2./7.);
		species1AndWall->setSlidingFrictionCoefficient(mu1W);
		species1AndWall->setSlidingStiffness(species1AndWall->getStiffness()*2.0/7.0);

		species1AndWall->setRollingStiffness(species1AndWall->getStiffness()*2.0/7.0);
		species1AndWall->setRollingFrictionCoefficient(mur1W);
		species1AndWall->setRollingDissipation(species1AndWall->getDissipation()*2./7.);

        //******************************SPECIES 2 AND WALL***********************************************
		auto species2AndWall = speciesHandler.getMixedObject(species2,speciesWall);
		//(Stiffness, CoR, effective mass)
		//Note: for final version, consider switching out to set tc instead of k
		species2AndWall->setStiffnessAndRestitutionCoefficient(k, e1W, 2*mass2);

        ///////////////////////////////ADHESIVE PROPERTIES///////////////////////////////////////////
	    species2AndWall->setAdhesionStiffness(k_a2W);
	    species2AndWall->setAdhesionForceMax(f_a2W*species2AndWall->getAdhesionStiffness());

        ///////////////////////////////FRICTIONAL PROPERTIES/////////////////////////////////////////
		species2AndWall->setSlidingDissipation(species2AndWall->getDissipation()*2./7.);
		species2AndWall->setSlidingFrictionCoefficient(mu2W);
		species2AndWall->setSlidingStiffness(species2AndWall->getStiffness()*2.0/7.0);

		species2AndWall->setRollingStiffness(species2AndWall->getStiffness()*2.0/7.0);
		species2AndWall->setRollingFrictionCoefficient(mur2W);
		species2AndWall->setRollingDissipation(species2AndWall->getDissipation()*2./7.);


        /////////////////////////////////////////////////////////////////////////////////////////////
        ////////////////////////////////SETTING SYSTEM PARAMETERS////////////////////////////////////
        /////////////////////////////////////////////////////////////////////////////////////////////
        setSystemDimensions(3);
        setGravity(Vec3D(0,0,-9.81));
        //The initial drum rotation rate - set to 0 to give initially stationary drum
        double RPSInitial = 0.0;
        //The particle radius
        Mdouble particleRadius = 1e-1;

        /////////////////////////////////////////////////////////////////////////////////////////////
        ////////////////////////////////SETTING SYSTEM GEOMETRY (WALLS)//////////////////////////////
        /////////////////////////////////////////////////////////////////////////////////////////////
   		wallHandler.clear();

        setXMin(0.0);
        setYMin(0.0);
        setZMin(0.0);

        setXMax(2. * drumRadius);
        setZMax(2. * drumRadius);
        setYMax(drumLength);

        //Positioning the centre of the drum
		Vec3D drumCenter = {0.5*(getXMin() + getXMax()),
							0.5*(getYMin() + getYMax()),
							0.5*(getZMin() + getZMax())};

		auto drumWall = wallHandler.copyAndAddObject(AxisymmetricIntersectionOfWalls());
		drumWall -> setSpecies(speciesWall);
		drumWall -> setPosition(drumCenter);
		drumWall -> setOrientation(Vec3D(0.0,1.0,0.0));
		drumWall -> addObject(Vec3D(1,0,0), Vec3D(drumRadius,0.0,0.0));
		drumWall -> setAngularVelocity(Vec3D(0.0,RPSInitial * 2.0 * constants::pi,0.0));

		InfiniteWall w0;
		w0.setSpecies(speciesWall);

		w0.set(Vec3D(0.,-1.,0.),Vec3D(drumCenter.X,getYMin(),drumCenter.Z));
		wallHandler.copyAndAddObject(w0);
		w0.set(Vec3D(0.,1.,0.),Vec3D(drumCenter.X,getYMax(),drumCenter.Z));
		wallHandler.copyAndAddObject(w0);

        /////////////////////////////////////////////////////////////////////////////////////////////
        ////////////////////////////////CREATING PARTICLES///////////////////////////////////////////
        /////////////////////////////////////////////////////////////////////////////////////////////
        
        //Declaring a generalised base particle, not yet assigned species or any other properties
        SphericalParticle P;
                
        //Counters to keep track of the number of each species already inserted
        //Initialised at zero as no particles initially inserted
        double nSpecies1_In = 0;
        double nSpecies2_In = 0;
        
        //Randomly inserting particles until the desired number has been inserted
        while ( (nSpecies1_In < nSpecies1) || (nSpecies2_In < nSpecies2)) {
			
			//determining the number left to be inserted of each species
			double nSpecies1_Left = nSpecies1 - nSpecies1_In;
			double nSpecies2_Left = nSpecies2 - nSpecies2_In;
			
			//creating a boolean "flag" to let the code know which species is being inserted
			bool species1_Flag;
			
			//setting the species and radius corresponding to the particle type of which there are
			//currently fewest
			if (nSpecies1_Left >= nSpecies2_Left) {
				//giving the particle the radius and species corresponding to the relevant species
				P.setSpecies(species1);
				//setting particle radius to that of correct species and simultaneously introducing polydispersity
				particleRadius = random.getRandomNumber((1. - polydispersity)*radiusSpecies1,(1. + polydispersity)*radiusSpecies1);
				P.setRadius(particleRadius);
				species1_Flag = true;
			}
			else {
				//giving the particle the radius and species corresponding to the relevant species
				P.setSpecies(species2);
				//setting particle radius to that of correct species and simultaneously introducing polydispersity
				particleRadius = random.getRandomNumber((1. - polydispersity)*radiusSpecies2,(1. + polydispersity)*radiusSpecies2);
				P.setRadius(particleRadius);
				species1_Flag = false;
			}

			//declaring a set of values that can be used to randomly place particles
			//**within** the confines of the drum!
			double r = random.getRandomNumber(1.0*particleRadius,drumRadius - 1.0*particleRadius);
			double theta = random.getRandomNumber(-constants::pi,constants::pi);
			double y = random.getRandomNumber(getYMin()+1.0*particleRadius,getYMax()-1.0*particleRadius);
						
			Vec3D position;
			position.X = drumRadius + r*cos(theta);
			position.Y = y;
			position.Z = drumRadius + r*sin(theta);
			
			std::cout << position << std::endl;
			
			//Checking that the value changes appropriately!
			P.setPosition(position);
			P.setVelocity(Vec3D(0.0,0.0,0.0));
			    
			//checking for overlaps before insertion
			if(checkParticleForInteraction(P)) {
				//putting the particle in the system
				particleHandler.copyAndAddObject(P);
				//incrementing nSpecies_In for the relevant species, so that the code know the 
				//particle has been added
				if (species1_Flag) {
					nSpecies1_In++;
				}
				else {
					nSpecies2_In++;
				}
				//Rebuilding the hGrid for the current set of particles
				hGridRebuild();
			}
			else{
				std::cout << "Error - collision - particle not inserted" << std::endl;
			}
		}
	}

    void actionsBeforeTimeStep() override
	{

		if (step==2)
		{
			if (getTime() > checkTime)
			{
				std::cout << "Current KE" << getKineticEnergy() << std::endl;
				if (getKineticEnergy() < startingKE)
				{
					step = 3;
					double drumStartTime = getTime();
					std::cout << "\n \n \n";
					std::cout << "-------------------------" << std::endl;
					std::cout << "Starting drum rotation" << std::endl;
					std::cout << "-------------------------" << std::endl;
					std::cout << "\n\n\n";
					// rotate the drum
					wallHandler.getObject(0)->setAngularVelocity(Vec3D(0.0,revsPerSecond*constants::pi*2.0,0.0));
					wallHandler.getObject(1)->setAngularVelocity(Vec3D(0.0,revsPerSecond*constants::pi*2.0,0.0));
					wallHandler.getObject(2)->setAngularVelocity(Vec3D(0.0,revsPerSecond*constants::pi*2.0,0.0));
				}
				else
				{
					checkTime = getTime() + 1.0;
				}
			}
		}
	}
    private:
    //TODO: Take out the initialisation at 2 - this is simply because other stuff is not yet added in!
    //Also, explain what this does!
    int step = 2;
    //The time at which we wish thee rotation to begin, all other factors being equal.
    double checkTime = 2.0;
    //The maximum kinetic energy (KE) at which it is acceptable to begin rotation of the drum
    //This ensures that the drum is 'settled' before motion commences.
    double startingKE = 10000000000;
    //The speed of rotation in revolutions per second
    double revsPerSecond;

	//*****************************************Particle Properties**********************************************
	//The sizes (radii) of the different species of particle used
	double radiusSpecies1;
	double radiusSpecies2;
	
	//The material densities of the different species of particle used
	double densitySpecies1;
	double densitySpecies2;
	
	//The (fractional) polydispersity of particles
	double polydispersity;
	
	//The (sliding) frictional coefficients for all species and mixed species
	//species1-species1 interactions
	double mu11;
	//species2-species2 interactions
	double mu22;
	//species1-species2 interactions
	double mu12;
	//species1-wall interactions
	double mu1W;
	//species12-wall interactions
	double mu2W;
	
	//The (rolling) frictional coefficients for all species and mixed species
	//species1-species1 interactions
	double mur11;
	//species2-species2 interactions
	double mur22;
	//species1-species2 interactions
	double mur12;
	//species1-wall interactions
	double mur1W;
	//species12-wall interactions
	double mur2W;
	
	//the elastic coefficients (coefficients of restitution) for all species and mixed species
	//species1-species1 interactions
	double e11;
	//species2-species2 interactions
	double e22;
	//species1-species2 interactions
	double e12;
	//species1-wall interactions
	double e1W;
	//species12-wall interactions
	double e2W;
	
	//The stiffness of particles (for simplicity and simulation accuracy, using same stiffness for all particles
	//and walls
	double k;
	
	//The adhesion stiffness coefficients for all interactions
	//species1-species1 interactions
	double k_a11;
	//species2-species2 interactions
	double k_a22;
	//species1-species2 interactions
	double k_a12;
	//species1-wall interactions
	double k_a1W;
	//species12-wall interactions
	double k_a2W;
	//Wall
	double k_aWW;
	
	//The adhesion maximal force values for all interactions
	//Note that f_a is a *MULTIPLIER*, not an absolute value - 
	//The force is equal to f_a * adhesionStiffness!
	//species1-species1 interactions
	double f_a11;
	//species2-species2 interactions
	double f_a22;
	//species1-species2 interactions
	double f_a12;
	//species1-wall interactions
	double f_a1W;
	//species12-wall interactions
	double f_a2W;
	//Wall
	double f_aWW;
	//*****************************************System Properties**********************************************
	//Setting the radius of the drum
    double drumRadius;
	//Setting the length of the drum
    double drumLength;
    //The number of particles of each species requested by the user
    int nSpecies1 = 30;
    int nSpecies2 = 30;
	public:
	//**********************************************************************************************************
	//*****************************************Getters and Setters**********************************************
	//**********************************************************************************************************
	//A function to set the radii of both particle species
	void setParticleRadii(double s1, double s2) {
		radiusSpecies1 = s1;
		radiusSpecies2 = s2;
	}
	//A function to set the material densities of both particle species
	void setParticleDensities(double d1, double d2) {
		densitySpecies1 = d1;
		densitySpecies2 = d2;
	}
	//Setting the (fractional) polydispersity of particles
	void setPolydispersity(double p) {
		polydispersity = p;
	}
	//Setting the sliding frictional coefficients for all interactions - individual and mixed
	void setSlidingFrictionCoefficients(double s11, double s22, double s12, double s1W, double s2W) {
		mu11 = s11;
		mu22 = s22;
		mu12 = s12;
		mu1W = s1W;
		mu2W = s2W;
	}
	//Setting the rolling frictional coefficients for all interactions - individual and mixed
	void setRollingFrictionCoefficients(double s11, double s22, double s12, double s1W, double s2W) {
		mur11 = s11;
		mur22 = s22;
		mur12 = s12;
		mur1W = s1W;
		mur2W = s2W;
	}
	//Setting the restitution coefficients for all interactions - individual and mixed
	void setRestitutionCoefficients(double s11, double s22, double s12, double s1W, double s2W) {
		e11 = s11;
		e22 = s22;
		e12 = s12;
		e1W = s1W;
		e2W = s2W;
	}
	//Setting the stiffness of the system components
	void setAllStiffnesses(double K) {
		k = K;
	}
	//Setting the adhesion stiffness for all species and mixed species
	void setAdhesionStiffnesses(double s11, double s22, double s12, double s1W, double s2W, double sWW) {
		k_a11 = s11;
		k_a22 = s22;
		k_a12 = s12;
		k_a1W = s1W;
		k_a2W = s2W;
		k_aWW = sWW;
	}
	//Setting the maximum adhesion force for all species and mixed species
	void setAdhesionMaximalForces(double s11, double s22, double s12, double s1W, double s2W, double sWW) {
		f_a11 = s11;
		f_a22 = s22;
		f_a12 = s12;
		f_a1W = s1W;
		f_a2W = s2W;
		f_aWW = sWW;
	}
    //setting the size of the drum
    void setDrumLengthAndRadius(double l, double r) {
		drumLength = l;
		drumRadius = r;
	} 
	//setting the rotation rate of the drum 
	void setRPM(double rpm) {
		revsPerSecond = rpm / 60.0;
	}
	//setting the (absolute) filling fraction and the relative fraction of each species
	void setFillFractionAndSpecies1VolumeFraction(double ff, double vf) {
		
		//sanity check
		if ( (vf > 1.0) || (ff > 1.0) ) {
			std::cerr << "Error: Fill fraction and Volume Fraction must be smaller less than or equal to 1"
						<< std::endl;
			exit(EXIT_FAILURE);
		}
		
		//calculating the total volume of the system
		double systemVolume = constants::pi * drumRadius * drumRadius * drumLength;

		//calculating the volumes of an individual particle of each species
		double particle1Volume = (4.0/3.0) * constants::pi * radiusSpecies1 * radiusSpecies1 * radiusSpecies1; 
		double particle2Volume = (4.0/3.0) * constants::pi * radiusSpecies1 * radiusSpecies1 * radiusSpecies1; 

		//calculating the required volume of each individual species
		double volumeNeeded1 = ff * vf * systemVolume;
		double volumeNeeded2 = ff * (1.0 - vf) * systemVolume;
		
		//determining the required number of particles to give the required volume
		nSpecies1 = 1 + volumeNeeded1 / particle1Volume;
		nSpecies2 = 1 + volumeNeeded2 / particle2Volume;
		
		std::cout << "number of particles requested: " << nSpecies1 << '\t' << nSpecies2 << std::endl;
		//exit(EXIT_SUCCESS);
		
		
	}
};

int main(int argc UNUSED, char *argv[] UNUSED)
{
	//Instantiating the CoDrum class
	CoDrum coDrum;

	//**********************************************************************************************************
	//************************************Setting Name and Key Parameters***************************************
	//**********************************************************************************************************
	//a stringstream to hold the full name
	std::stringstream nameStream;






	
	//**********************************************************************************************************
	//************************************Setting Particle Properties*******************************************
	//**********************************************************************************************************
	//Choosing the radii of particles
	coDrum.setParticleRadii(5e-3,5e-3); //(species1 radius, species2 radius)
	//Choosing the densities of particles
	coDrum.setParticleDensities(2000,2000); //(species1 density, species2 density)
	//Choosing the degree of polydispersity (0 -> monodisperse, 1 -> radius varies from zero to 2*r)
	coDrum.setPolydispersity(0.05);
	//Setting the sliding and rolling friction parameters for all species and mixed-species interactions
	coDrum.setSlidingFrictionCoefficients(0.1,0.1,0.1,0.6,0.6); //For species combinations (1-1, 2-2, 1-2, 1-wall, 2-wall)
	coDrum.setRollingFrictionCoefficients(0.01,0.01,0.01,0.3,0.3); //For species combinations (1-1, 2-2, 1-2, 1-wall, 2-wall)
	//setting the elastic properties (resitution coefficient) for all species and mixed-species interactions
	coDrum.setRestitutionCoefficients(0.5, 0.5, 0.5, 0.1, 0.1); //For species combinations (1-1, 2-2, 1-2, 1-wall, 2-wall)
	//setting the adhesion stiffnesses and adhesion maximal forces for all species and combinations thereof
	//Note that if we want to keep the *EXTENT* of the force constant, we need to ensure f_a/k_a remains constant
	coDrum.setAdhesionStiffnesses(1.0e-20,1.0e-20,1.0e-20,1.0e-20,1.0e-20, 1.0e-20); //(1-1, 2-2, 1-2, 1-wall, 2-wall, wall)
	//Note: this sets the multiplier value (see above), i.e. the value here actually gives the *LENGTH* of the interaction
	coDrum.setAdhesionMaximalForces(1.0e-5,1.0e-20,1.0e-20,1.0e-20,1.0e-20, 1.0e-20); //(1-1, 2-2, 1-2, 1-wall, 2-wall, wall) 
	//setting the ('normal', normal) stiffnesses of all particles and walls to the same value
	coDrum.setAllStiffnesses(1.0e4);
	//**********************************************************************************************************
	//************************************Setting System Parameters*********************************************
	//**********************************************************************************************************
	//setting the radius of the rotating drum
	coDrum.setDrumLengthAndRadius(.05,.15);
	//setting the strength and direction of the gravitational acceleration
	coDrum.setGravity(Vec3D(0.,0.,-9.81));
	//setting the (absolute) filling fraction and the relative fraction of each particle species
	coDrum.setFillFractionAndSpecies1VolumeFraction(0.25,0.5);
	//setting the drum's rotation rate in RPM
	coDrum.setRPM(10.0);
	//**********************************************************************************************************
	//************************************Setting Simulation Parameters*********************************************
	//**********************************************************************************************************
	coDrum.setName("1adhesive-k_a1e-20-f_a1e-5");
	coDrum.setTimeStep(1.0 / 8000.0);
	coDrum.setTimeMax(500.0);
	coDrum.setFileType(FileType::ONE_FILE);
	coDrum.setSaveCount(100);
	
	coDrum.setXBallsAdditionalArguments("-cmode 8 -solidf -v0");
	coDrum.solve();

	return 0;
}
