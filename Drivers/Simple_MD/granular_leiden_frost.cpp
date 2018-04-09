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


#include "scr/Mercury2D.h"


//#define DEBUG_OUTPUT



class my_problem : public Mercury2D{
public:
	
void setupInitialConditions()
{



for (int i=0;i<get_N();i++)
	{


	double x=random(getXMin(),getXMax());
	double y=random(getYMin(),getYMax());
	

	particleHandler.getObject(i)->Position=Vec3D(x,y,0.0);
	particleHandler.getObject(i)->setVelocity(Vec3D(0.00,0.00,0.0));
	particleHandler.getObject(i)->getRadius()=0.0001;

	}	

set_NWall(2);
Walls[1].set(Vec3D(0,1,0), getYMax());
Walls[0].set(Vec3D(0, -1,0),  -getYMin());
set_NWallPeriodic(1);
WallsPeriodic[0].set(Vec3D(1,0,0), getXMin(), getXMax());
}

protected:

void actionsBeforeTimeStep()
{
double time=getTime()-0.1;
if (time>0.0)
	{
	Walls[0].move(getYMin()-0.0001*sin(time*400.0));
	set_dissipation(1e-2);
	}

}


		
};




/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int main(int argc UNUSED, char *argv[] UNUSED)
{



	///Start off my solving the default problem
 	my_problem problem;
 	problem.set_N(1000);
 	problem.set_name("leidenfrost");
 	problem.set_dissipation(1);
	problem.setSaveCount(100);
 	problem.setTimeMax(2);
	problem.setTimeStep(1e-5);
 	problem.solve();


	
	
	
	
 
		
}
