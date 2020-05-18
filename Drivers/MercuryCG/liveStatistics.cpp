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
#include <CG/CG.h>
#include "Mercury3D.h"
using namespace std;

class Mercury3DRestart : public Mercury3D
{
	void actionsBeforeTimeStep() override {};
	void setupInitialConditions() override {setTime(time);}
public:
	double time;
};

int main(int argc, char *argv[])
{
//	if (argc<2) {
//		cout << "Please enter the name of the simulation you want to restart and, optionally, the name of the simulation after restart" << endl;
//		exit(-1);
//	} else {
//		cout << "restart data: " << argv[1] << ".restart" << endl;
//	}
 	Mercury3DRestart problem;
 	problem.setName("ShearCell3D_vol_75nl_surf7_30degree");
	problem.readRestartFile();
	problem.setName("Out");
	problem.setRestarted(false);
	problem.time = problem.getTime();
	problem.write(std::cout,false);

//	auto cg2 = problem.cgHandler.copyAndAddObject(CG<CGCoordinates::R,CGFunctions::Gauss>());
//	cg2->setN(100);
//	cg2->setWidth(0.001);
//	cg2->setTimeMin(problem.getTimeMax());


	problem.solve(argc-1, argv+1);
	for (BaseInteraction* p : problem.interactionHandler) {
		if (p->getOverlap()>0) logger(INFO,"overlap %",p->getOverlap());
	}
	logger(INFO,"mean overlap %",problem.interactionHandler.getMeanOverlap());
	logger(INFO,"mean radius  %",problem.particleHandler.getMeanRadius());
}
