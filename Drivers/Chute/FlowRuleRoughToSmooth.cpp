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

#include "scr/DPMBase.h"
#include <iostream>
#include <vector>
using namespace std;

class flowrule : public DPMBase {
public:
	void run()
	{
		//~ //Check if the run has been done before
		//~ if (FileExists(data_filename.str()))
		//~ {
			//~ //If it has move on to teh next run immedently
			//~ cout << "This run has been done " << endl;
			//~ launch_new("parm_demo",true);
			//~ exit(0);
		//~ }
		//~ 
		//~ //Set up a 3 by 1 study
		//~ vector<int> study_num=get_numbers(3,1);
		//~ 
		//~ //If study 0 is complete quit
		//~ if (study_num[0] > 0)
			//~ {
				//~ cout << "Study is complete " << endl;
				//~ exit(0);
			//~ }
		//~ else
		//~ //If the study is not complete save the data to disk and move on
			//~ {
				//~ save_info_to_disk();
				//~ launch_new("FlowRuleRoughToSmooth");
			//~ }
		//~ 
	}
};

int main(int argc UNUSED, char *argv[] UNUSED)
{
	///Start off by solving the default problem
	flowrule problem;
	
	///Autonumber turns on file numbering
	problem.auto_number();
	problem.setName("flowrule");
	//problem.run();
}
