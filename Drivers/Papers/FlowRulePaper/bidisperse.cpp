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
#include <cstring>

#include "SilbertPeriodic.h"

class FlowRule : public SilbertPeriodic {
public:

	//save only every 100th restart data in a new file
  //	void writeRestartFile() {
  //		static int counter=0;
  //		if (!(counter%100)) {
  //			//writeRestartFile_counter--;
  //			MD::writeRestartFile();
  //		}
  //		counter++;
  //	}

	void setName() {
	  std::stringstream name;
		name << "H" << getInflowHeight() 
		<< "A" << getChuteAngleDegrees() 
		<< "P" << getPolydispersity()
		<< "D" << getDensityVariation()
		<< "N" << NumberFraction;
		dataFile.setName(name.str());
	}
	
	void run(int argc, char *argv[])
	{
		//Set up a parameter study
		setSaveCount(10000);
		setTimeMax(2000);
		restartFile.setFileType(FileType::MULTIPLE_FILES);
		dataFile.setFileType(FileType::NO_FILE);
		fStatFile.setFileType(FileType::NO_FILE);
                eneFile.setFileType(FileType::ONE_FILE);
		//restartFile.setFileType(FileType::MULTIPLE_FILES_PADDED);
		//dataFile.setFileType(FileType::NO_FILE);
		//fStatFile.setFileType(FileType::NO_FILE);
		//eneFile.setFileType(FileType::ONE_FILE);

		//1 is the base species
		NumberFraction = 0.5;
		createBaseSpecies();
		setPolydispersity(1);
		setDensityVariation(1);
		
		readArguments(argc, argv);
		setName();
	
		//Save the info to disk
		writeRestartFile();

		//Check if the run has been done before
		if (helpers::fileExists(dataFile.getName())) {
			//If it has move on to the next run immedently
		  std::cout << "Run " << getName() << " has already been done " << std::endl;
		} else {
			//launch this code
			solve();
		}
	}
	
	///Sets variable values for particles that are created at the inflow
	///Chooses between small species 0 and large species 1 according to the number fraction
	void create_inflow_particle()
	{
	  if (random.getRandomNumber(0,1)<NumberFraction)
	    inflowParticle_.setSpecies(speciesHandler.getObject(1));
	  else
	    inflowParticle_.setSpecies(speciesHandler.getObject(0));
	  if (inflowParticle_.getIndSpecies()==0)
	    inflowParticle_.setRadius(getMinInflowParticleRadius());
	  else
	    inflowParticle_.setRadius(getMaxInflowParticleRadius());
	  //P0.indSpecies_+ = random(0,1)<NumberFraction;
	  //P0.getRadius() = P0.indSpecies_+ ? getMaxInflowParticleRadius() : getMinInflowParticleRadius();
	  //inflowParticle_.computeMass();

	  Vec3D position;
	  position.X = random.getRandomNumber(getXMin()+2.0* inflowParticle_.getRadius(),getXMax());
	  position.Y = random.getRandomNumber(getYMin()+2.0* inflowParticle_.getRadius(),getYMax());
	  position.Z = random.getRandomNumber(getZMin()+2.0* inflowParticle_.getRadius(),getInflowHeight());
		
	  inflowParticle_.setPosition(position);
	  inflowParticle_.setVelocity(Vec3D(0.0,0.0,0.0));
	}
	
	/// Sets radii such that the total volume remains constant, and rmax/rmin = Polydispersity
	/// Also sets FixedParticleRadius to MinInflowParticleRadius to avoid Particles falling through
	void setPolydispersity(double Polydispersity) {
		if (Polydispersity>=1.) {
		  setMinInflowParticleRadius(getFixedParticleRadius()*std::pow((NumberFraction+mathsFunc::cubic(Polydispersity)*(1-NumberFraction)),-1./3.));
			setMaxInflowParticleRadius(getMinInflowParticleRadius()*Polydispersity);
			setFixedParticleRadius(getMinInflowParticleRadius());
			std::cout 
			<< "r0=" << getMinInflowParticleRadius()
			<< "r1=" << getMaxInflowParticleRadius()
			<< "Numberfraction=" << NumberFraction
			<< std::endl;
	} else {
		  std::cerr<<"Error: polydispersity " << Polydispersity << " needs to be >=1"<<std::endl; exit(1);
		}
	}

	/// Changes density of (small) species 0
	void setDensityVariation(double Densityvariation) {
	  if (Densityvariation>0) speciesHandler.getObject(0)->setDensity(speciesHandler.getObject(1)->getDensity()*Densityvariation);
		else {std::cerr<<"Error: densityvariation needs to be positive"<<std::endl; exit(-1);}
		std::cout 
		  << "rho0=" << speciesHandler.getObject(0)->getDensity()
		  << "rho1=" << speciesHandler.getObject(1)->getDensity()
		<< std::endl;
	}

        double getDensityVariation(){return speciesHandler.getObject(0)->getDensity()/speciesHandler.getObject(1)->getDensity();}
	
	double getPolydispersity(){return getMaxInflowParticleRadius()/getMinInflowParticleRadius();}

private:

	/// allows input from the command line      
	bool readNextArgument(int& i, int argc, char *argv[]) override {
		if (!strcmp(argv[i],"-polydispersity")) {
			setPolydispersity(atof(argv[i+1]));
		} else if (!strcmp(argv[i],"-densityvariation")) {
			setDensityVariation(atof(argv[i+1]));
		} else return SilbertPeriodic::readNextArgument(i, argc, argv); 
		return true; //returns true if argv[i] is found here
	}
		
	/// number ratio of (large) species 1
	double NumberFraction;
	
	// if you want to give VolumeFraction instead of Number fraction, use thefollowing conversion:	
	// double NumberFraction = 1.-1./(1.+VolumeFraction/(1.-VolumeFraction)*mathsFunc::cubic(get_Polydispersity()));

};

int main(int argc, char *argv[])
{
	FlowRule problem;
	problem.run(argc, argv);
}
