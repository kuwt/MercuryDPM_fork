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

#include "../Papers/FlowRulePaper/SilbertPeriodic.h"
#include "ChuteBottom.h"
using namespace std;

class FlowRule : public SilbertPeriodic {
public:

	void launch(bool startfast=false) {
		launchNewRun("flowRule",startfast);
	}

	//save only every 100th restart data in a new file
  //	void save_restart_data() {
  //		static int counter=0;
  //		if (!(counter%100)) {
  //			//cout << endl << counter << endl;
  //			//save_restart_data_counter--;
  //			MD::save_restart_data();
  //		}
  //		counter++;
  //	}

  void run(std::vector<Mdouble> study_num, int argc, char *argv[])
	{
		//Set up a parameter study
		setInflowHeight(study_num[1]);
		setChuteAngle(study_num[2]);
		set_study(study_num[0]);
		//
		restartFile.setFileType(FileType::ONE_FILE);//restartFile.setFileType(FileType::ONE_FILE);
		dataFile.setFileType(FileType::ONE_FILE);
		fStatFile.setFileType(FileType::ONE_FILE);
		eneFile.setFileType(FileType::ONE_FILE);
		readArguments(argc, argv);
		//set_study();	
		
		//Save the info to disk
		writeRestartFile();

		//Check if the run has been done before
		// if (FileExists(data_filename.str())) {
// 			//If it has move on to teh next run immedently
// 			cout << "Run " << getName() << " has already been done " << endl;
// 		} else 
		{
			//launch this code
		  std::stringstream com("");
			com << "echo started \tstudy_num \t" 
				 << study_num[0] << " " 
				 << study_num[1] << " " 
				 << study_num[2] << " \tname \t" 
				 << getName() << " &>>ReportFlowRule";
			std::cout << system(com.str().c_str()) << std::endl;
			std::cout << "started study_num " 
				 << study_num[0] << " " 
				 << study_num[1] << " " 
				 << study_num[2] << ", name " 
			     << getName() << std::endl;
		//new commands to make it 2D
		setChuteWidth(2);
		setChuteLength(5);
		setInflowParticleRadius(0.47,0.53);
		setSaveCount(2e3);
		setTimeMax(30);
		//set_NWallPeriodic(1);
		//setRoughBottomType(MONOLAYER_ORDERED);
			solve();

			com.str("");
			com << "echo finished \tstudy_num \t" 
				 << study_num[0] << " " 
				 << study_num[1] << " " 
				 << study_num[2] << " \tname \t" 
				 << getName() << " &>>ReportFlowRule";
			std::cout << system(com.str().c_str()) << std::endl;
		}
		

	}
	
// 	void createBottom() {
// 		if (abs(getFixedParticleRadius())<1e-12) // smooth bottom
// 		{
// 			//bottom wall 
// 			Walls.resize(Walls.size()+1);
// 			Walls.back().set(Vec3D(0.0, 0.0, -1.0), -get_zmin());
// 		}
// 		else //rough bottom
// 		{
// 			//this pointer is the current MD class, so the bottom is create with the particles properties from the MD class
// 			ChuteBottom bottom(*this);
// 			bottom.setThickness(3);
// 			bottom.setInflowParticleRadius(getFixedParticleRadius());
// 			bottom.makeRoughBottom(getParticleHandler());
// 			cout<<"Starting to destruct ChuteBottom"<<endl;
// 			for (vector<Particle*>::iterator it= getParticleHandler().begin(); it!=getParticleHandler().end(); ++it) 
// 				(*it)->fixParticle();
// 		}
// 	}
	

};

int main(int argc, char *argv[])
{
	FlowRule problem;
	problem.restartFile.setFileType(FileType::MULTIPLE_FILES);//restartFile.setFileType(FileType::MULTIPLE_FILES_PADDED);
	problem.setSaveCount(1e4);
	problem.restartFile.setSaveCount(100.0*problem.restartFile.getSaveCount());
	problem.setTimeMax(2000);
	//~ problem.setRoughBottomType(MONOLAYER_DISORDERED);

	//set case, height, angle to given or default values
	vector<Mdouble> study_num;
	study_num.resize(3);

	//this line is needed for the code to work with demoparams
	//if (argc<4) exit(-1);
	
	if (argc>1) study_num[0]=atof(argv[1]);
	else {      study_num[0]=4;  std::cout << "No study number given; using lambda=1" << std::endl;} 
	if (argc>2) study_num[1]=atof(argv[2]);
	else {      study_num[1]=5; std::cout << "No height given; using h=10 d" << std::endl;} 
	if (argc>3) study_num[2]=atof(argv[3]);
	else {      study_num[2]=28; std::cout << "No angle given; using theta=24 deg" << std::endl;} 
	
	//call run
	if (argc>3) problem.run(study_num,argc-3,argv+3);
	else problem.run(study_num,1,argv);
}
