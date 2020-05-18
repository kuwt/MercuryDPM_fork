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

	void launch(bool startfast=false) {
		launchNewRun("flowRule",startfast);
	}

	void run(vector<int> study_num, int argc, char *argv[])
	{
		//if counter == -1 we do a parameter study
		bool dostudy=false;
		if (study_num.size()!=3) {
			dostudy=true;
			autoNumber();
			//Set up a 2 parameter study
			study_num=get2DParametersFromRunNumber(4/*#Heights*/,9/*#Angles*/);
		}

		set_study(study_num);
		readArguments(argc, argv);
		
		//If the study is not complete save the data to disk and move on
		writeRestartFile();
		
		//This can be used to exclude certain numbers from running
		bool launchthiscode = true;

		//start next code
		if (dostudy) {
			if (launchthiscode) {
				launch();
			} else {
				cout << "Code " << get_counter() << " was NOT launched" << endl;
				launch(true);
			}
		}		

		//Check if the run has been done before
		if (FileExists(data_filename.str()))
		{
			//If it has move on to teh next run immedently
			cout << "Run " << getName() << " has already been done " << endl;
			exit(0);
		}
		
		//launch this code
		stringstream com("");
		com << "echo started \tstudy_num \t" 
		     << study_num[0] << " " 
		     << study_num[1] << " " 
		     << study_num[2] << " \tname \t" 
		     << getName() << " &>>ReportFlowRule";
		int sysret;
		sysret = system(com.str().c_str());
		cout << "started study_num " 
		     << study_num[0] << " " 
		     << study_num[1] << " " 
		     << study_num[2] << ", name " 
		     << getName() << endl;

		if (launchthiscode) solve();

		com.str("");
		com << "echo finished \tstudy_num \t" 
		     << study_num[0] << " " 
		     << study_num[1] << " " 
		     << study_num[2] << " \tname \t" 
		     << getName() << " &>>ReportFlowRule";
		sysret = system(com.str().c_str());
	}
};

int main(int argc, char *argv[])
{
	FlowRule problem;
	problem.restartFile.setFileType(FileType::MULTIPLE_FILES_PADDED);
	//~ problem.setRoughBottomType(MONOLAYER_DISORDERED);
	//~ problem.setTimeMax(1e-1);
	//problem.readArguments(argc, argv);
	vector<int> study_num;
	if (argc>3) {
		study_num.resize(3);
		study_num[0]=atoi(argv[1]);
		study_num[1]=atoi(argv[2]);
		study_num[2]=atoi(argv[3]);
	}
	problem.run(study_num,argc-3,argv+3);
}
