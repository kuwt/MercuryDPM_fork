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

#include "scr/Mercury3D.h"

#include <iostream>

class finite_walls_hopper : public Mercury3D{
public:
	
  //void computeInternalForces(CParticle &PI, CParticle &PJ) {}
	
  	
  void setupInitialConditions()
  {
    double radius = 0.001; //glass beds jaeger's
    double dist = 2.0 * radius;
    int N = floor(getXMax()/dist);

        
    int id = 0;
    for (int i=1; i<=N; i++) {
      for (int j=1; j<=N; j++) {
	for (int k=1; k<=10; k++) {
	  if (sqrt(mathsFunc::square(i*dist - getXMax()/2.0)+mathsFunc::square(j*dist-getYMax()/2.0)) < (50-1)*radius){
	    id++;
	  }
	}
      }
    }
    set_N(id);
    std::cout << "Number of particles: " << get_N() << std::endl;
    id = 0;

    

   
    //Particles[id].set_mass(1);
    
    for (int i=1; i<=N; i++) {
      for (int j=1; j<=N; j++) {
	for (int k=1; k<=10; k++) {
	  if (sqrt(mathsFunc::square(i*dist - getXMax()/2.0)+mathsFunc::square(j*dist-getYMax()/2.0)) < (50-1)*radius){
	    Particles[id].Radius=radius;
	    Particles[id].Position=Vec3D(i*dist - getXMax()/2.0,j*dist  - getYMax()/2.0,2*radius + (k-1)*dist);
	    Particles[id].Velocity=Vec3D(1*random(-.1,.1),1*random(-.1,.1),random(-.01,.01));
	    //Particles[id].set_mass(1e-5);
	    id++;
	    }
	}
      }
    }


    set_NWall(1);
    Walls[0].set_is_hopper();
    //Walls[0].set(Vec3D( 0.0,0.0, -1.0), 0 );
  }
	
};


int main(int argc UNUSED, char *argv[] UNUSED)
{
	
  ///Start off my solving the default problem
  finite_walls_hopper problem;
  problem.setSystemDimensions(3);
  problem.set_name("finite_walls_hopper");
  problem.setTimeMax(.2);
  problem.setTimeStep(5e-6);
  //problem.setSaveCount(helpers::getSaveCountFromNumberOfSavesAndTimeMaxAndTimeStep(1000,problem.getTimeMax(),problem.getTimeStep()));
  problem.setSaveCount(25);
	
  
  problem.setXMax(.1);
  problem.setYMax(.1);
  problem.setZMax(.1);
  problem.setXMin(-.1);
  problem.setYMin(-.1);
  problem.setZMin(-1);
  problem.setGravity(Vec3D(0,0,-9.8));
  

  problem.set_dissipation(1);
  // problem.speciesHandler.getObject(0)->setStiffness(500);
  problem.setStiffnessAndRestitutionCoefficient(500,0.001,8e-7);
  std::cout <<"Collision time "<< problem.getCollisionTime(8e-7) << std::endl;
  
  problem.solve();
  problem.save_info_to_disk();
 
  std::cout << problem << std::endl;
  
  //now run ./finite_walls/finite_walls.disp -rmult .05
}
