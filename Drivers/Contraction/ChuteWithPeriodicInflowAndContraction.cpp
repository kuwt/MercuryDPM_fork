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

#include "ChuteWithPeriodicInflow.h"
#include "Walls/IntersectionOfWalls.h"

class ChuteWithPeriodicInflowAndContraction : public ChuteWithPeriodicInflow
{
public:
  ChuteWithPeriodicInflowAndContraction(std::string restart_file) : ChuteWithPeriodicInflow(restart_file) {
		//dimensions 700x130x10
		//contraction at 300-500x13
		
		//corrected to dimensions 260x130x4
		//contraction at 60-260x13
		double xContractionStart =  60;
		double xContractionEnd   = 140;
		double yContractionStart =  50;
		double yContractionEnd   =  25;
		
		//defines how often the inflow has to be copied
		double lengthPeriodicChute = getXMax(); //length of restarted periodic chute (typically 20); will also be length of inflow
		double widthPeriodicChute = getYMax(); //width of restarted periodic chute (typically 10);
		int numRepetitions = round(xContractionEnd/getXMax()); //length 13x20=260
		int numRepetitionsInWidth = round(yContractionStart/getYMax()); //width 13x10=130
		
		//creates enough particles to fill contraction
		//set_Nmax(get_N()*numRepetitionsInWidth*(numRepetitions+2.));
		int N=particleHandler.getNumberOfObjects();
		//extend inflow (species 0) in width
		for (int k=1; k<numRepetitionsInWidth; k++) {
			for (int i=0; i<N; i++) {
			  particleHandler.addObject(particleHandler.getObject(i));//Particles.push_back(Particles[i]);
			  particleHandler.getLastObject()->move(Vec3D(0.0,widthPeriodicChute*k,0.0));//Particles.back().Position.Y+=widthPeriodicChute*k;
			}
		}
		//extend contraction (species 1) in width and height
		for (int k=0; k<numRepetitionsInWidth; k++) {
			for (int j=1; j<numRepetitions; j++) {
				for (int i=0; i<N; i++) {
					if (inflowParticle_.isFixed()) {
					  particleHandler.addObject(particleHandler.getObject(i));//Particles.push_back(Particles[i]);
					  particleHandler.getLastObject()->setSpecies(speciesHandler.getObject(1));
					  particleHandler.getLastObject()->move(Vec3D(lengthPeriodicChute*j,0.,0.));//Particles.back().Position.X+=lengthPeriodicChute*j;
					  particleHandler.getLastObject()->move(Vec3D(0.,widthPeriodicChute*k,0.));//Particles.back().Position.Y+=widthPeriodicChute*k;
					}
				}
			}
		}
		std::cout << "Particles added: " << particleHandler.getNumberOfObjects() << std::endl;

		setXMax((numRepetitions+5)      *lengthPeriodicChute);
		setYMax(numRepetitionsInWidth* widthPeriodicChute);
		PeriodicBoundary* perw = static_cast<PeriodicBoundary*> (boundaryHandler.getObject(1));
		perw->set(Vec3D( 0.0, 1.0, 0.0), getYMin(), getYMax());
		//
		//fix_hgrid();
		//set_HGRID_num_buckets_to_power(get_Nmax()); // automated already
		//~ Walls.resize(1);
		//~ Walls.back().set(Vec3D(0,0,-1),getInflowParticleRadius());

		set_symmetric_contraction(xContractionStart,xContractionEnd,(yContractionStart-yContractionEnd)/2.);
		// remove particles in contraction wall
		double dist;
		Vec3D normal;
		bool touch;
		int counter=0;
		for (int i=0; i<particleHandler.getNumberOfObjects(); i++) {
		  for (int k=wallHandler.getNumberOfObjects()-2; k<wallHandler.getNumberOfObjects(); k++) {
		    touch = wallHandler.getObject(k)->getDistanceAndNormal(*particleHandler.getObject(i), dist, normal);
				if (touch) {
					counter++;
					particleHandler.removeObject(i);
					i--;
					break;
				}
			}
		}
		std::cout << "Particles deleted: " << counter << std::endl;
		
		setName("ChuteWithPeriodicInflowAndContraction");
		write(std::cout,false);
	}
	
	void set_symmetric_contraction(double x_min, double x_max, double delta_y) {
	  IntersectionOfWalls wall;
	  //Walls.resize(Walls.size()+1);
		//back wall
		Vec3D normalIntoWall = Vec3D(-1,0,0);
		Vec3D point = Vec3D(x_max,0,0);
		wall.addObject(normalIntoWall, Vec3D::dot(normalIntoWall,point));
		//slanted wall
		double delta_x = x_max-x_min;
		normalIntoWall = Vec3D(delta_y,-delta_x,0)/sqrt(mathsFunc::square(delta_x)+mathsFunc::square(delta_y));
		point = Vec3D(x_min,0,0);
		wall.addObject(normalIntoWall, Vec3D::dot(normalIntoWall,point));

		//Walls.resize(Walls.size()+1);
		//back wall
		normalIntoWall = Vec3D(-1,0,0);
		point = Vec3D(x_max,getChuteWidth(),0);
		wall.addObject(normalIntoWall, Vec3D::dot(normalIntoWall,point));
		//slanted wall
		delta_x = x_max-x_min;
		normalIntoWall = Vec3D(delta_y,delta_x,0)/sqrt(mathsFunc::square(delta_x)+mathsFunc::square(delta_y));
		point = Vec3D(x_min,getChuteWidth(),0);
		wall.addObject(normalIntoWall, Vec3D::dot(normalIntoWall,point));
		wallHandler.copyAndAddObject(wall);
	}
	
};



int main(int argc, char* argv[])
{
	ChuteWithPeriodicInflowAndContraction problem("../ini/H10A28L2M0.5B0.5");
	problem.setTimeMax(1e4);
	//problem.setSaveCount(1./problem.getTimeStep());
	problem.setSaveCount(5e3);
	problem.setXBallsAdditionalArguments("-v0 -solidf -h 600 -w 1400 ");
	//
	problem.dataFile.setFileType(FileType::ONE_FILE);
	problem.restartFile.setFileType(FileType::ONE_FILE);
	problem.fStatFile.setFileType(FileType::ONE_FILE);
	problem.eneFile.setFileType(FileType::ONE_FILE);
	//problem.restartFile.setFileType(FileType::ONE_FILE);
	//problem.dataFile.setFileType(FileType::ONE_FILE);
	//problem.fStatFile.setFileType(FileType::NO_FILE);
	//problem.eneFile.setFileType(FileType::ONE_FILE);
	problem.solve(argc, argv);
}
