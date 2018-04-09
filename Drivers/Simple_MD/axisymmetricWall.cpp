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

#include "Mercury3D.h"
#include "Matrix.h"
#include "AxisymmetricIntersectionOfWalls.h"


class axisymmetricWall : public Mercury3D {

public:

	void setupInitialConditions()
	{
		setXMin(0);
		setYMin(0);
		setZMin(0);
		setXMax(1);
		setYMax(1);
		setZMax(1);

		int N = 10;
		set_N(N*N*N);
		
		double Radius = (getXMax()-getXMin())/N/2.;
		
		//set Particles' position, radius, velocity and bounding box
		vector<CParticle>::iterator it = Particles.begin();
		for (int i=0;i<N;i++)
		for (int j=0;j<N;j++)
		for (int k=0;k<N;k++)
		{
			it->Radius = Radius;
			it->setVelocity(Vec3D(0.0,0.0,0.0));
			it->Position = Radius*Vec3D(1.+2.*i, 1.+2.*j, 1.+2.*k);
			it->Position.Z+=1;
			it++;
		}	
		
		//set walls
		set_NWall(2);
		//for a prism wall, define points ABC (or more points) 
		//as the corners of the prism base (ordered in clockwise 
		//direction) and D as the 'periodic direction' of the prism 
		//(upwards: Vec3D::Cross(A-B,D) has to point into the wall)
		Vec3D A(0.1,0.0,0.0), B(0.4,0.0,0.8), C(0.4,0.0,0.0);
		Vec3D D(0.0,1.0,0.0); //Periodic direction of the prism
		Walls[0].addObject(Vec3D::cross(A-B,D),A);
		Walls[0].addObject(Vec3D::cross(B-C,D),B);
		Walls[0].addObject(Vec3D::cross(C-A,D),C);
		//~ Walls[0].set(Vec3D(1.0,0.0,1.0), .25*getXMax());
		Walls[0].setPosition(Vec3D(0.5,0.5,0),Vec3D(0,0,1));

		Walls[1].set(Vec3D(0,0,-1),2*Radius);
		set_NWallPeriodic(0);
		
		//create particle properties
		particleHandler.getObject(0)->computeMass(Species);
		double tc = .1;
		setCollisionTimeAndRestitutionCoefficient(tc, .05, particleHandler.getObject(0)->getMass());
		setTimeStep(tc/50.);
		setTimeMax(200.*tc);
		setSaveCount(50);
		setGravity(Vec3D(0,0,-1));
		write(std::cout,false);
	}
	
};

int main(int argc UNUSED, char *argv[] UNUSED)
{
	axisymmetricWall md;
	md.set_name("axisymmetricWall");
	//~ std::cout << md.getTimeStep() << "t" << md.getTimeMax() << std::endl;
	md.solve(argc,argv);
}

