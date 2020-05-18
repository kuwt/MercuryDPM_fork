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

#include "SilbertPeriodic.h"
using namespace std;

class FlowRule : public SilbertPeriodic {
public:

	//save only every 100th restart data in a new file
	void writeRestartFile() {
		static int counter=0;
		if (!(counter%100)) {
			//writeRestartFile_counter--;
			DPMBase::writeRestartFile();
		}
		counter++;
	}

	void setname() {
		stringstream name;
		name << "H" << getInflowHeight() 
		<< "A" << getChuteAngleDegrees() 
		<< "P" << get_Polydispersity()
		<< "D" << get_Densityvariation()
		<< "N" << NumberFraction;
		setName(name.str().c_str());
		set_data_filename();
	}
	
	void run(int argc, char *argv[])
	{
		//Set up a parameter study
		setSaveCount(1e4);
		setTimeMax(2000);
		restartFile.setFileType(FileType::MULTIPLE_FILES_PADDED);
		dataFile.setFileType(FileType::NO_FILE);
		fStatFile.setFileType(FileType::NO_FILE);
		eneFile.setFileType(FileType::ONE_FILE);

		//1 is the base species
		NumberFraction = 0.5;
		createBaseSpecies();
		set_Polydispersity(1);
		set_Densityvariation(1);
		
		readArguments(argc, argv);
		setname();
	
		//Save the info to disk
		writeRestartFile();

		//Check if the run has been done before
		if (FileExists(data_filename.str())) {
			//If it has move on to the next run immedently
			cout << "Run " << getName() << " has already been done " << endl;
		} else {
			//launch this code
			solve();
		}
	}
	
	///Sets variable values for particles that are created at the inflow
	///Chooses between small species 0 and large species 1 according to the number fraction
	void create_inflow_particle()
	{
		P0.indSpecies_+ = random.get_RN(0,1)<NumberFraction;
		P0.Radius = P0.indSpecies_+ ? getMaxInflowParticleRadius() : getMinInflowParticleRadius();
		P0.computeMass(Species);

		P0.Position.X = random.get_RN(getXMin()+2.0*P0.Radius,getXMax());
		P0.Position.Y = random.get_RN(getYMin()+2.0*P0.Radius,getYMax());
		P0.Position.Z = random.get_RN(getZMin()+2.0*P0.Radius,getInflowHeight());
		P0.Velocity = Vec3D(0.0,0.0,0.0);
	}
	
	/// Sets radii such that the total volume remains constant, and rmax/rmin = Polydispersity
	/// Also sets FixedParticleRadius to MinInflowParticleRadius to avoid Particles falling through
	void set_Polydispersity(Mdouble Polydispersity) {
		if (Polydispersity>=1.) {
			setMinInflowParticleRadius(getFixedParticleRadius()*pow((NumberFraction+mathsFunc::cubic(Polydispersity)*(1-NumberFraction)),-1./3.));
			setMaxInflowParticleRadius(getMinInflowParticleRadius()*Polydispersity);
			setFixedParticleRadius(getMinInflowParticleRadius());
			cout 
			<< "r0=" << getMinInflowParticleRadius()
			<< "r1=" << getMaxInflowParticleRadius()
			<< "Numberfraction=" << NumberFraction
			<< endl;
	} else {
			cerr<<"Error: polydispersity " << Polydispersity << " needs to be >=1"<<endl; exit(1);
		}
	}

	/// Changes density of (small) species 0
	void set_Densityvariation(Mdouble Densityvariation) {
		if (Densityvariation>0) setDensity(getDensity(1)*Densityvariation,0);
		else {cerr<<"Error: densityvariation needs to be positive"<<endl; exit(-1);}
		cout 
		<< "rho0=" << getDensity(0)
		<< "rho1=" << getDensity(1)
		<< endl;
	}

	Mdouble get_Densityvariation(){return getDensity(0)/getDensity(1);}
	
	Mdouble get_Polydispersity(){return getMaxInflowParticleRadius()/getMinInflowParticleRadius();}

private:

	/// allows input from the command line
	int readNextArgument(unsigned int& i, unsigned int argc, char *argv[]) {
		if (!strcmp(argv[i],"-polydispersity")) {
			set_Polydispersity(atof(argv[i+1]));
		} else if (!strcmp(argv[i],"-densityvariation")) {
			set_Densityvariation(atof(argv[i+1]));
		} else return SilbertPeriodic::readNextArgument(i, argc, argv); 
		return true; //returns true if argv[i] is found here
	}
		
	/// number ratio of (large) species 1
	Mdouble NumberFraction;
	
	// if you want to give VolumeFraction instead of Number fraction, use thefollowing conversion:	
	// Mdouble NumberFraction = 1.-1./(1.+VolumeFraction/(1.-VolumeFraction)*mathsFunc::cubic(get_Polydispersity()));

};

int main(int argc, char *argv[])
{
	FlowRule problem;
	problem.run(argc, argv);
}
