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

#include<iostream>
#include "scr/DPMBase.h"
#include "scr/Mercury3D.h"
#include "scr/Chute.h"
using namespace std;

double ratio;
double ChuteAngleInc;
int savecount;

class ChutePeriodic : public Chute{
public:

	void actionsBeforeTimeStep(){};
		
	void setupInitialConditions()
	{
		stringstream com("");
		com << getName() << ".ini";
		readDataFile(com.str().c_str());

		t=time;
		write(std::cout,false);
	}
	
	void writeEneTimeStep(std::ostream& os) const {
		double ene_kin = 0, ene_rot = 0, ene_gra = 0;

		for (unsigned int i=0;i<Particles.size();i++) if (!Particles[i].is_fixed())
		{
			ene_kin += .5 * Particles[i].get_mass() * Particles[i].Velocity.GetLength2();
			ene_rot += Particles[i].getRotationalEnergy();
			ene_gra -= Particles[i].get_mass() * Vec3D::Dot(getGravity(),Particles[i].Position);
		} //end for loop over Particles

		///todo{Why is there a +6 here?}
		static int width = ene_file.precision() + 6;
		ene_file  << setw(width) << getTime() 
			<< " " << setw(width) << ene_gra
			<< " " << setw(width) << ene_kin
			<< " " << setw(width) << ene_rot
			<< " " << setw(width) << getElasticEnergy()
			<< " " << setw(width) << getChuteAngleDegrees()
			<< endl;
		
		static double t_change = t;
		if (ene_kin/getElasticEnergy()<ratio) {
			t_change = t;
			cout << "ChuteAngle: " << getChuteAngleDegrees();
			setChuteAngle(getChuteAngleDegrees()+ChuteAngleInc);
			cout << " to " << getChuteAngleDegrees() << ", ene_rat: " << ene_kin/getElasticEnergy() << endl;
		} else if (t-t_change>300&&t>1000) exit(0);
		resetElasticEnergy();
	}

	double time;
};

int main(int argc, char *argv[])
{
	if (argc<2) {
		cout << "Please enter the name of the simulation you want to restart and, optionally, the name of the simulation afted restart" << endl;
		exit(-1);
	} else {
		cout << "restart data: " << argv[1] << ".restart" << endl;
	}
	
	savecount = atoi(argv[2]);
	ChuteAngleInc = atof(argv[3]);
	ratio = atof(argv[4]);
	cout << "savecount " << savecount
		<< ", ChuteAngleInc " << ChuteAngleInc
		<< ", ratio " << ratio
		<< endl;
	
 	ChutePeriodic problem;
 	problem.setName(argv[1]);
	problem.load_restart_data();
	problem.readArguments(argc-4, argv+4);
	problem.time = problem.getTime();
	problem.setSaveCount(savecount);
	problem.setTimeMax(1e20);
	problem.dataFile.setFileType(FileType::NO_FILE);
	problem.fStatFile.setFileType(FileType::NO_FILE);

	problem.solve();
	
	problem.write(std::cout,false);
	//problem.writeRestartFile();
}
