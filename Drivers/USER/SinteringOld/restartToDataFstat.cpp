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

#include<iostream>
#include "Sinter2.h"
using namespace std;

class Restart : public Sinter
{
public:

	void actionsBeforeTimeStep() override {};

	double time;

	void setupInitialConditions() override {
        setTime(time);}
};

int main(int argc, char *argv[])
{
	if (argc<2) {
		cout << "Please enter the name of the simulation you want to restart and, optionally, the name of the simulation after restart" << endl;
		exit(-1);
	} else {
		cout << "restart data: " << argv[1] << endl;
	}
 	Restart problem;
    problem.readRestartFile(argv[1]);
	problem.setRestarted(false);
	problem.restartFile.setFileType(FileType::NO_FILE);
	problem.eneFile.setFileType(FileType::NO_FILE);
	problem.dataFile.setFileType(FileType::MULTIPLE_FILES);
	problem.fStatFile.setFileType(FileType::MULTIPLE_FILES);
	problem.time = problem.getTime();
	problem.setTimeMax(problem.getTime());
	problem.write(std::cout,false);
	problem.solve(argc-1, argv+1);
}
