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

#include "scr/Time.h"
#include "scr/Mercury3D.h"

#include <cmath>
#include <iostream>
#include <iomanip>


class my_problem : public Mercury3D {

public:

	void setupInitialConditions()
	{
		int N = round(pow(Particles.size(),1.0/3.0));
		set_N(mathsFunc::cubic(N));
		double r = 0.001;
		
		double sqr3 = sqrt(3), sqr2 = sqrt(2);
		vector<CParticle>::iterator it = Particles.begin();
		for (int i=0;i<N;i++)
		for (int j=0;j<N;j++)
		for (int k=0;k<N;k++)
		{
			it->Radius = r;
			it->Position = Vec3D(r+2.0*r*i+r*(j%2)*(k%2?-1.0:1.0)+r*(k%2), r+sqr3*r*j+2.0/sqr3*r*(k%2), r+2.0*sqr2/sqr3*r*k);
			it->setVelocity(Vec3D(0.0,0.0,0.0));
			it++;
		}	
		
		setXMax(2.0*r*N+r);
		setYMax(sqr3*r*N+(2.0-sqr3)*r+2.0/sqr3*r);
		setZMax(2.0*sqr2/sqr3*r*N+(2.0-2.0*sqr2/sqr3)*r);
		
		set_NWall(6);
		Walls[0].set(Vec3D( 1.0,0.0,0.0), getXMax());
		Walls[1].set(Vec3D(-1.0,0.0,0.0),-getXMin());
		Walls[2].set(Vec3D(0.0, 1.0,0.0), getYMax());
		Walls[3].set(Vec3D(0.0,-1.0,0.0),-getYMin());
		Walls[4].set(Vec3D(0.0,0.0, 1.0), getZMax());
		Walls[5].set(Vec3D(0.0,0.0,-1.0),-getZMin());
		set_NWallPeriodic(0);
	}
	
};


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int main(int argc UNUSED, char *argv[] UNUSED)
{
	///Start off my solving the default problem
	my_problem problem;
	problem.set_N(mathsFunc::cubic(1));
	problem.setSystemDimensions(3);
	problem.set_dissipation(1.0);
	problem.set_name("hgrid");
	problem.setupInitialConditions();
	problem.setTimeStep();
	problem.setTimeMax(problem.getTimeStep()*100.0);
	problem.setSaveCount(helpers::getSaveCountFromNumberOfSavesAndTimeMaxAndTimeStep(0,getTimeMax(),getTimeStep()));
	problem.solve();
	//problem.writeRestartFile();
	//std::cout << problem;

	Time time;
	double t = 0, t_previous = 0;
	int n_previous = 0;
	for (int n=5; n<=100; n+=5) {
		time.tic();
		problem.set_N(mathsFunc::cubic(n));
		problem.solve();
		t_previous = t;
		t = time.toc();
		std::cout << "N=" << setw(5) << problem.get_N() << ", " 
			<< "time=" << setw(5) << t << ", ";
		if (t_previous) std::cout << "slope log(time)/log(N) " << setprecision(3) << (log(t)-log(t_previous))/(log(mathsFunc::cubic(n))-log(mathsFunc::cubic(n_previous))) << std::endl;
		n_previous = n;
	}
	std::cout << std::endl;
	
}
