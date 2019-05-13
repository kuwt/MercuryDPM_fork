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
#include "InfiniteWallWithHole.h"
#include "Screw.h"
 
class Screw_test : public Mercury3D{
	private:
	
	void setupInitialConditions()
	{
        wallHandler.clear();
		InfiniteWall w0;
		w0.set(Vec3D(-1, 0, 0), -getXMin());
		wallHandler.copyAndAddObject(w0);
		w0.set(Vec3D( 1, 0, 0),  getXMax());
		wallHandler.copyAndAddObject(w0);
		w0.set(Vec3D( 0,-1, 0), -getYMin());
		wallHandler.copyAndAddObject(w0);
		w0.set(Vec3D( 0, 1, 0),  getYMax());
		wallHandler.copyAndAddObject(w0);
		w0.set(Vec3D( 0, 0,-1), -getZMin());
		wallHandler.copyAndAddObject(w0);
		
		InfiniteWallWithHole* infiniteWallWithHole=new InfiniteWallWithHole(Vec3D(0,0,1),getZMax(),1.0);
		wallHandler.addObject(infiniteWallWithHole);
		
        ///Here screw properties are set
		///Screw(Start position,Length,Radius,Number of turns,Rotation speed,Thickness)
		screw=new Screw(Vec3D(0,0,0),1.0,1.0-ParticleRadius,1.0,-1.0,0.5*ParticleRadius);
		wallHandler.addObject(screw);
		
		particleHandler.clear();
		SphericalParticle p0;
		p0.setVelocity(Vec3D(0.0,0.0,0.0));
		p0.setRadius(BaseParticleRadius);
		p0.computeMass(Species);

        
        //Single test case
        /*double distance;
		Vec3D normal;
        p0.setPosition(Vec3D(1.0,0.0,0.0));
        if(screw->getDistance_and_normal(p0, distance, normal))
			std::cout<<"Collision, distance screw="<<distance<<std::endl;
		else
			std::cout<<"No collision, distance screw="<<distance<<std::endl;
        */

        //Simple run settings
        //Nx*Ny*Nz particles are created evenly spaced between [xmin,xmax]*[ymin,ymax]*[zmin,zmax] and checked for contact with the screw
		int Nx=floor((getXMax()-getXMin())/(2.1*ParticleRadius));
		int Ny=floor((getYMax()-getYMin())/(2.1*ParticleRadius));
		int Nz=floor((getZMax()-getZMin())/(2.1*ParticleRadius));
		Mdouble distance;
		Vec3D normal;		
		for(int i=0;i<Nx;i++)
			for(int j=0;j<Ny;j++)
				for(int k=0;k<Nz;k++)
				{
					p0.setPosition(Vec3D(getXMin()+(getXMax()-getXMin())*(0.5+i)/Nx,getYMin()+(getYMax()-getYMin())*(0.5+j)/Ny,getZMin()+(getZMax()-getZMin())*(0.5+k)/Nz));
                    if(!screw->getDistance_and_normal(p0, distance, normal))
                    {
						particleHandler.copyAndAddObject(p0);
                    }
                    else
                    {
                        //std::cout<<p0.getPosition()<<std::endl;
                    }
				}
	}
	
	void actionsBeforeTimeStep()
    {
	  if (getTime()>1)
			screw->move_time(getTimeStep());
			
	}
	
	public:
		double ParticleRadius;
        Screw* screw;

};

int main(int argc UNUSED, char *argv[] UNUSED)
{
	Screw_test problem;
	
	problem.set_name("Screw_test");
	problem.setDensity(1000);
	problem.setXMax(1.0);
	problem.setYMax(5.0);
	problem.setZMax(2.0);
	problem.setXMin(-1.0);
	problem.setYMin(-1.0);
	problem.setTimeMax(2.0);
	problem.ParticleRadius=0.2;	
	problem.speciesHandler.getObject(0)->setCollisionTimeAndRestitutionCoefficient(0.05,0.8, pow(problem.ParticleRadius,3)*constants::pi*4.0/3.0*1000);
    problem.setTimeStep(0.02*0.05);
	problem.setSaveCount(helpers::getSaveCountFromNumberOfSavesAndTimeMaxAndTimeStep(1000,problem.getTimeMax(),problem.getTimeStep()));
	problem.write(std::cout,false);
	problem.solve();
}
