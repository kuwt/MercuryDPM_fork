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


#include "Mercury3D.h"
#include "Species/LinearViscoelasticFrictionSpecies.h"
#include "Walls/AxisymmetricIntersectionOfWalls.h"
#include "Walls/InfiniteWall.h"
#include "Boundaries/PeriodicBoundary.h"
#include <string>
#include <sstream>

class CoDrum : public Mercury3D{
public:

	void setupInitialConditions() override
	{

	    auto species1 = speciesHandler.copyAndAddObject(LinearViscoelasticFrictionSpecies());	
	    auto species2 =  speciesHandler.copyAndAddObject(LinearViscoelasticFrictionSpecies());
	    auto speciesWall =  speciesHandler.copyAndAddObject(LinearViscoelasticFrictionSpecies());
  

	  setParticleDimensions(3);
        /////////////////////////////////////////////////////////////////////////////////////////////
        ////////////////////////////////Setting Species Properties///////////////////////////////////
        /////////////////////////////////////////////////////////////////////////////////////////////

	  //calculating masses of all species and effective mass of particles 1 and 2 or use in
	  //later calculations
	  double mass1 = 4.0 / 3.0 * constants::pi * pow(radiusSpecies1, 3.0) * densitySpecies1;
	  double mass2 = 4.0 / 3.0 * constants::pi * pow(radiusSpecies2, 3.0) * densitySpecies2;
	  double mEff = 2*mass1*mass2 / (mass1 + mass2);

	  //Letting user know masses of particles - allows to check that if equal to unity for dimensionless systems
	  std::cout << "Mass 1: " << mass1 << std::endl;
	  std::cout << "Mass 2: " << mass2 << std::endl;
	  std::cout << "Effective mass: " << mEff << std::endl;
	  std::cout << "Average mass: " << (mass1 + mass2) / 2  << std::endl;

	  
        //******************************SPECIES 1****************************************************
        species1->setDensity(densitySpecies1);
        //TO FIX - GET ACCURATE MASS!
        species1->setStiffnessAndRestitutionCoefficient(k,e11,mass1);
	   

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
       

        ///////////////////////////////FRICTIONAL PROPERTIES/////////////////////////////////////////
        species2->setSlidingStiffness(species2->getStiffness()*2./7.);
        species2->setSlidingFrictionCoefficient(mu22);
        species2->setSlidingDissipation(species2->getDissipation()*2./7.);

		species2->setRollingStiffness(species2->getStiffness()*2.0/7.0);
		species2->setRollingFrictionCoefficient(mur22);
		species2->setRollingDissipation(species2->getDissipation()*2./7.);


        //******************************WALL SPECIES****************************************************
        speciesWall->setDensity(densitySpecies1);
        //TO FIX - GET ACCURATE MASS!
        speciesWall->setStiffnessAndRestitutionCoefficient(k,e11,mEff);

       
        ///////////////////////////////FRICTIONAL PROPERTIES/////////////////////////////////////////
        speciesWall->setSlidingStiffness(speciesWall->getStiffness()*2./7.);
        speciesWall->setSlidingFrictionCoefficient(mu22);
        speciesWall->setSlidingDissipation(speciesWall->getDissipation()*2./7.);

		speciesWall->setRollingStiffness(speciesWall->getStiffness()*2.0/7.0);
		speciesWall->setRollingFrictionCoefficient(mur22);
		speciesWall->setRollingDissipation(speciesWall->getDissipation()*2./7.);

        /////////////////////////////////////////////////////////////////////////////////////////////
        ////////////////////////////////Setting Mixed Properties/////////////////////////////////////
        /////////////////////////////////////////////////////////////////////////////////////////////

        //******************************SPECIES 1 AND 2***********************************************
		auto species1And2 = speciesHandler.getMixedObject(species1,species2);
		//(Stiffness, CoR, effective mass)
		//Note: for final version, consider switching out to set tc instead of k
		species1And2->setStiffnessAndRestitutionCoefficient(k, e12, mEff);


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
		species1AndWall->setStiffnessAndRestitutionCoefficient(k, e1W, mass1);


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
		species2AndWall->setStiffnessAndRestitutionCoefficient(k, e1W, mass2);

 
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
        setGravity(Vec3D(0,0,-1));
        //The initial drum rotation rate - set to 0 to give initially stationary drum
        double RPSInitial = 0.0;
        //The particle radius
        Mdouble particleRadius = 5e-1;

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

		PeriodicBoundary b0;
		b0.set(Vec3D(0.0,1.0,0.0),getYMin(),getYMax());
		boundaryHandler.copyAndAddObject(b0);
		
		/*
		InfiniteWall w0;
		w0.setSpecies(speciesWall);

		w0.set(Vec3D(0.,-1.,0.),Vec3D(drumCenter.X,getYMin(),drumCenter.Z));
		wallHandler.copyAndAddObject(w0);
		w0.set(Vec3D(0.,1.,0.),Vec3D(drumCenter.X,getYMax(),drumCenter.Z));
		wallHandler.copyAndAddObject(w0);
		*/

        /////////////////////////////////////////////////////////////////////////////////////////////
        ////////////////////////////////CREATING PARTICLES///////////////////////////////////////////
        /////////////////////////////////////////////////////////////////////////////////////////////
        
        //Declaring a generalised base particle, not yet assigned species or any other properties
        SphericalParticle P;
                
        //Counters to keep track of the number of each species already inserted
        //Initialised at zero as no particles initially inserted
        double nSpecies1_In = 0;
        double nSpecies2_In = 0;

	//Counter to track the number of failed attempts to insert a particle
	int failCounter = 0;
	
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
				//if successfully inserted, reset fail counter
				failCounter = 0;
			}
			else{
			  //Counting the number of failed attempts
			  failCounter++;
			  if (failCounter > 10000) {
			    std::cerr << "Warning - difficulties inserting particle - attempt "
				      << failCounter << std::endl;
			  }
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
					//wallHandler.getObject(1)->setAngularVelocity(Vec3D(0.0,revsPerSecond*constants::pi*2.0,0.0));
					//wallHandler.getObject(2)->setAngularVelocity(Vec3D(0.0,revsPerSecond*constants::pi*2.0,0.0));
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
    double checkTime = 5.0;
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
  //A function to automatically set the sizes and densities of particles based on the particle
  //such that the mean mass and mean particle diameter are equal to unity. Requires as input
  //only the desired particle size ratio and the factor by which particle density should be increased
  void setNondimensional(double ratio, double dFactor) {
    //First, ensuring that the mean sizes of the particle species equal unity
    radiusSpecies1 = ratio / (1 + ratio);
    radiusSpecies2 = 1 / (1 + ratio);

    //Outputting to allow user to check that mean diameter is correct
    std::cout << "Size species 1 = " << 2 * radiusSpecies1 << std::endl;
    std::cout << "Size species 2 = " << 2 * radiusSpecies2 << std::endl;
    std::cout << "Mean particle diameter = " << radiusSpecies1 + radiusSpecies2 << std::endl;

    //Next, setting the particle density such that the mean particle mass is equal to unity
    densitySpecies1 = dFactor * 3.0 / (2 * constants::pi *
			   (radiusSpecies1*radiusSpecies1*radiusSpecies1 +
			    radiusSpecies2*radiusSpecies2*radiusSpecies2 ) );
    
    densitySpecies2 = densitySpecies1;

    //Outputting particle masses to check that their mean equals unity
    double m1 = (4.0/3.0) * constants::pi * radiusSpecies1*radiusSpecies1*radiusSpecies1 * densitySpecies1;
    double m2 = (4.0/3.0) * constants::pi * radiusSpecies2*radiusSpecies2*radiusSpecies2 * densitySpecies2;

    std::cout << "Mass species 1 = " << m1 << std::endl;
    std::cout << "Mass species 2 = " << m2 << std::endl;
    std::cout << "Mean particle mass = " << (m1 + m2) / 2.0 << std::endl;
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
    double particle2Volume = (4.0/3.0) * constants::pi * radiusSpecies2 * radiusSpecies2 * radiusSpecies2; 
    
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
	//The basis of the name of the file
	std::string nameBase = "frictionTest";
	//the multiplier by which particle density is increased
        double dens = 1.0; 	
	//The drum size
	double drumDiam = 35.0;
	//The relevant frictional properties
	double mu_p_s = 0.1; //the sliding friction coefficient for particle-particle interactions (Silbert glass = 0.1)
	double mu_p_r = 1.0e100; //rolling friction for p-p interactions (Silbert glass = 0.01)
	double mu_w_s = mu_p_s; //sliding friction for particle-wall interactions
	double mu_w_r = mu_p_r; //rolling friction for p-w interactions
	// The desired Froude number
	double Fr = 2.0;
	//Calculating the relevant angular velocity from the desired Froude number
	double angV = sqrt(Fr * 2 / drumDiam);
	//converting to RPM
	double rotRate = angV * 60 / (2*constants::pi);
	//creating a full name via stringstream
	nameStream << nameBase << "-Fr" << Fr << "-mu_w_s" << mu_w_s
		   << "-mu_w_r" << mu_w_r;
	//setting the appropriate name of the file
	coDrum.setName(nameStream.str());

	
	//**********************************************************************************************************
	//************************************Setting Particle Properties*******************************************
	//**********************************************************************************************************
	//Setting particle size and density automatically to give mean mass of one and mean diameter of one at the
	//desired particle diameter ratio (given as the argument)
	coDrum.setNondimensional(1.0,dens);
	//Choosing the degree of polydispersity (0 -> monodisperse, 1 -> radius varies from zero to 2*r)
	coDrum.setPolydispersity(0.05);
	//Setting the sliding and rolling friction parameters for all species and mixed-species interactions
	coDrum.setSlidingFrictionCoefficients(mu_p_s,mu_p_s,mu_p_s,mu_w_s,mu_w_s); //For species combinations (1-1, 2-2, 1-2, 1-wall, 2-wall)
	coDrum.setRollingFrictionCoefficients(mu_p_r,mu_p_r,mu_p_r,mu_w_r,mu_w_r); //For species combinations (1-1, 2-2, 1-2, 1-wall, 2-wall)
	//setting the elastic properties (resitution coefficient) for all species and mixed-species interactions
	coDrum.setRestitutionCoefficients(0.91, 0.91, 0.91, 0.91, 0.91); //For species combinations (1-1, 2-2, 1-2, 1-wall, 2-wall) 
	//setting the ('normal', normal) stiffnesses of all particles and walls to the same value
	coDrum.setAllStiffnesses(1e5);
	//**********************************************************************************************************
	//************************************Setting System Parameters*********************************************
	//**********************************************************************************************************
	//setting the radius of the rotating drum
	coDrum.setDrumLengthAndRadius(5.0,drumDiam/2.0);
	//setting the strength and direction of the gravitational acceleration
	coDrum.setGravity(Vec3D(0.,0.,-1.0));
	//setting the (absolute) filling fraction and the relative fraction of each particle species
	coDrum.setFillFractionAndSpecies1VolumeFraction(0.2,0.5); //normally 0.2, 0.5
	//setting the drum's rotation rate in RPM
	coDrum.setRPM(rotRate);
	//**********************************************************************************************************
	//************************************Setting Simulation Parameters*********************************************
	//**********************************************************************************************************
	coDrum.setTimeStep(1.0 / 5000.0);
	//setting the time to give a set number of revolutions
	coDrum.setTimeMax((20.0*60.0 / rotRate)+5.0);
	coDrum.setFileType(FileType::ONE_FILE);
	//Setting to record a fixed number (5000) of individual time steps - 1000 per complete rotation
	coDrum.setSaveCount(coDrum.getTimeMax() / (5000 * coDrum.getTimeStep()));	
	coDrum.setXBallsAdditionalArguments("-cmode 0 -solidf -v0");
	coDrum.solve();

	return 0;
}
