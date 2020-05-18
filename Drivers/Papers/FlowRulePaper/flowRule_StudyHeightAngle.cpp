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

#include "SilbertPeriodic.h"

class FlowRule : public SilbertPeriodic {
public:   
    
	void launch(bool startfast=false) {
		launchNewRun("flowRule",startfast);
	}

    void run(std::vector<Mdouble> studyNumber, int argc, char *argv[])
	{
		//Set up a parameter study
		setInflowHeight(studyNumber[1]);
		setChuteAngle(studyNumber[2]);
		set_study(studyNumber[0]);          
		readArguments(argc, argv);
		//set_study();	
		
		//Save the info to disk
		writeRestartFile();

		//Check if the run has been done before
		//if (helpers::fileExists(dataFile.getName())) {
		//	//If it has move on to teh next run immedently
        //    std::cout << "Run " << getName() << " has already been done " << std::endl;
		//} else 
        {
			//launch this code
            std::stringstream com("");
			com << "echo started \tstudyNumber \t" 
				 << studyNumber[0] << " " 
				 << studyNumber[1] << " " 
				 << studyNumber[2] << " \tname \t" 
				 << getName() << " &>>ReportFlowRule";
			//std::cout << system(com.str().c_str()) << std::endl;
			std::cout << "started studyNumber " 
				 << studyNumber[0] << " " 
				 << studyNumber[1] << " " 
				 << studyNumber[2] << ", name " 
				  << getName() << std::endl;

			solve();

			com.str("");
			com << "echo finished \tstudyNumber \t" 
				 << studyNumber[0] << " " 
				 << studyNumber[1] << " " 
				 << studyNumber[2] << " \tname \t" 
				 << getName() << " &>>ReportFlowRule";
			//std::cout << system(com.str().c_str()) << std::endl;
		}
		
		//~ //launch next code
		//~ stringstream com("");
		//~ com.str("");
		//~ com << "./flowRule.exe " 
		     //~ << studyNumber[0] << " " 
		     //~ << studyNumber[1]+1 << " " 
		     //~ << studyNumber[2] << " &";
		//~ cout << studyNumber[1] << com.str() << endl;
		//~ if (studyNumber[1]<4) system(com.str().c_str());
	}
};

int main(int argc, char *argv[])
{
	FlowRule problem;
	//set case, height, angle to given or default values
	std::vector<Mdouble> studyNumber;
	studyNumber.resize(3);

	//this line is needed for the code to work with demoparams
	//if (argc<4) exit(-1);
	
	if (argc>3) 
    {
        studyNumber[0]=atof(argv[1]);
        studyNumber[1]=atof(argv[2]);
        studyNumber[2]=atof(argv[3]);
		problem.run(studyNumber,argc-3,argv+3);
    } 
    else
    {
        std::cout << "Not enough input arguments given (./flowRule_StudyHeightAngle $study $height $angle); " << std::endl
            << "using demo values (equivalent to ./flowRule_StudyHeightAngle 5 10 24 -tmax 0.01)" << std::endl;
        studyNumber[0]=5;
        studyNumber[1]=2;
        studyNumber[2]=24;
        problem.setTimeMax(0.01);
        problem.dataFile.setFileType(FileType::ONE_FILE);
        problem.setChuteLength(5); ///\todo make selftest smaller (currently 8 sec)); also CoilSelfTest (8s), LeesEdwardsSelfTest (13s), MaserSelfTest (14s)
        problem.setChuteWidth(5);
        //problem.setRoughBottomType(MULTILAYER);
        problem.setRoughBottomType(MONOLAYER_DISORDERED);
        //problem.setTimeMax(10);
        problem.run(studyNumber,1,argv);
        problem.setName("flowRuleSelfTest");
        problem.writeRestartFile();
    }
}
