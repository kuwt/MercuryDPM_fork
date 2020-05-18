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

#include<iostream>
#include "Chute.h"

class ChutePeriodic : public Chute{
public:

	void actionsBeforeTimeStep(){}
		
	void setupInitialConditions()
	{
		setTime(getTimeStep());
		setTimeMax(getTimeStep());
                dataFile.setFileType(FileType::ONE_FILE);
                restartFile.setFileType(FileType::NO_FILE);
                eneFile.setFileType(FileType::NO_FILE);
                fStatFile.setFileType(FileType::NO_FILE);
                write(std::cout,false);
	}

	void writeToLocalFolder() {
		//keep file name but create files in the local directory, i.e. remove folder
	  std::string name = getName();
		size_t found=name.find_last_of("/\\");
		setName(name.substr(found+1).c_str());
		std::cout << "new name: " << getName() << std::endl;
	}

};

int main(int argc, char *argv[])
{
	if (argc<2) {
	  std::cout << "snapshot.exe problemname [args]" << std::endl;
		exit(-1);
	} else {
	  std::cout << "restart data: " << argv[1] << ".restart" << std::endl;
	}
 	ChutePeriodic problem;
 	problem.setName(argv[1]);
	problem.readRestartFile();
	problem.writeToLocalFolder();
	problem.setXBallsScale(1.7/problem.getZMax());
	//problem.setXBallsColourMode(10);
	problem.setXBallsAdditionalArguments("-v0 -solidf -h 900 -w 1450 -o 120 -noborder 4 -rgbb 60 -rgbg 60 -rgbr 60 -rgbs 90");
	problem.readArguments(argc-1, argv+1);
	problem.solve();
	
	std::stringstream command("");
	command << "./" << problem.getName() << ".disp -of " << problem.getName() << ".pdf -die";
	int out = system (command.str().c_str());
	std::cout<<"The return value by system ="<<out<<std::endl;

}
