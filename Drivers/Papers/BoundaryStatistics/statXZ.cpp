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

#include "Chute.h"
#include "StatisticsVector.h"
bool verbose = false;

template <StatType T> class statistics_while_running : public StatisticsVector<T>, public Chute
{
	public:
	statistics_while_running<T>() : StatisticsVector<T>() , Chute() {}

	void actionsBeforeTimeStep(){};
		
	void setupInitialConditions() {
		if (verbose) StatisticsVector<T>::printStat();
		//~ for (unsigned int i=0; i<particleHandler.getNumberOfObjects(); i++) {
			//~ if (particleHandler.getObject(i)->isFixed() || i>10) {removeParticle(i); i--;}
		//~ }
	};
	
	void printTime() const {
		static Mdouble tint = getTimeMax()-getTime();
		std::cout << "\rProgress: " << std::setprecision(2) << std::setw(10) << (int)100.*(1-(getTimeMax()-getTime())/tint) << "%\r";
		std::cout.flush();
	}

	// allows getZMax() to be set to the height of the highest particle
	void auto_set_domain() {
		if (particleHandler.getNumberOfObjects()) {
		  Vec3D Max= particleHandler.getObject(0)->getRadius()*Vec3D(1,1,1) +particleHandler.getObject(0)->getPosition();//Particles[0].Position+Particles[0].Radius;
		  Vec3D Min=-particleHandler.getObject(0)->getRadius()*Vec3D(1,1,1) +particleHandler.getObject(0)->getPosition();//Particles[0].Position-Particles[0].Radius;
			for (unsigned int i=0; i<particleHandler.getNumberOfObjects(); i++) {
				Max=Vec3D::max(Max,particleHandler.getObject(i)->getPosition()+particleHandler.getObject(i)->getRadius()*Vec3D(1,1,1) );
				Min=Vec3D::min(Min,particleHandler.getObject(i)->getPosition()-particleHandler.getObject(i)->getRadius()*Vec3D(1,1,1) );
			}
			setXMax(Max.X+.05*(Max.X-Min.X));
			setXMin(Min.X-.05*(Max.X-Min.X));
			setYMax(Max.Y+.05*(Max.Y-Min.Y));
			setYMin(Min.Y-.05*(Max.Y-Min.Y));
			setZMax(Max.Z+.05*(Max.Z-Min.Z));
			setZMin(Min.Z-.05*(Max.Z-Min.Z));
		}
	}

	
};

int main(int argc, char *argv[])
{
	//select file name to restart from
  std::string name;
	if (argc>1) {name=argv[1]; argv++; argc--;} else exit(-1);
	//select time interval
	Mdouble tint=.2; if (argc>1) {tint=atof(argv[1]); argv++; argc--; std::cout << "tint=" << tint << std::endl;}
	argv--; argc++;
	
	//load restart data
	std::cout << "obtaining XZ statictics" << std::endl;
 	statistics_while_running<XZ> problem;
	std::cout << "loading restart data: " << name << ".restart" << std::endl;
 	problem.setName(name.c_str());
	problem.readRestartFile();//load_restart_data();
	problem.setFixedParticleRadius(problem.particleHandler.getObject(0)->getRadius());
	problem.setInflowParticleRadius(problem.particleHandler.getObject(0)->getRadius());
	problem.restartFile.setFileType(FileType::ONE_FILE); //problem.restartFile.setFileType(FileType::ONE_FILE);
	problem.writeRestartFile();
	problem.auto_set_domain();
	//keep file name but create files in the local directory, i.e. remove folder
	if (verbose) std::cout << "old name: " << problem.getName() << std::endl;
	size_t found=name.find_last_of("/\\");
	problem.setName(name.substr(found+1).c_str());
	if (verbose) std::cout << "new name: " << problem.getName() << std::endl;
	//set output to minimum
	problem.dataFile.setFileType(FileType::NO_FILE);
	problem.restartFile.setFileType(FileType::NO_FILE);
	problem.fStatFile.setFileType(FileType::NO_FILE);
	problem.eneFile.setFileType(FileType::ONE_FILE);
	//set statistical parameters
        problem.setDoPeriodicWalls(false);
 	problem.setN(50);
 	problem.setCGWidth(.1);
	problem.setSaveCount(25);
	problem.setStressTypeForFixedParticles(0);
 	problem.setCGTimeMin(problem.getTime());
	problem.setTimeMax(problem.getTime()+tint);
	//solve and create live statistics
	problem.readStatArguments(argc, argv);
	problem.solve();

	//~ cout << endl << "Z averaging" << endl;
 	//~ statistics_while_running<Z> problemZ;
 	//~ problemZ.set_name(name.c_str());
	//~ problemZ.load_restart_data();
	//~ problemZ.setFixedParticleRadius(problemZ.getObjects()[0].Radius);
	//~ problemZ.setInflowParticleRadius(problemZ.getObjects()[0].Radius);
	//~ problemZ.restartFile.setFileType(FileType::ONE_FILE);
	//~ problemZ.writeRestartFile();
	//~ problemZ.auto_set_domain();
	//~ //keep file name but create files in the local directory, i.e. remove folder
	//~ cout << "old name: " << problemZ.get_name() << endl;
	//~ stringstream ss;
	//~ ss << name.substr(found+1).c_str() << "Z";
	//~ problemZ.set_name(ss.str().c_str());
	//~ cout << "new name: " << problemZ.get_name() << endl;
	//~ //set output to minimum
	//~ problemZ.dataFile.setFileType(FileType::NO_FILE);
	//~ problemZ.restartFile.setFileType(FileType::ONE_FILE);
	//~ problemZ.fStatFile.setFileType(FileType::NO_FILE);
	//~ problemZ.eneFile.setFileType(FileType::ONE_FILE);
	//~ //set statistical parameters
    //~ problemZ.setDoPeriodicWalls(false);
 	//~ problemZ.setN(50);
 	//~ problemZ.setCGWidth(.1);
	//~ problemZ.setSaveCount(25);
	//~ problemZ.setStressTypeForFixedParticles(3);
 	//~ problemZ.setCGTimeMin(problemZ.get_t());
	//~ problemZ.setTimeMax(problemZ.get_t()+tint);
	//~ //solve and create live statistics
	//~ problemZ.readStatArguments(argc, argv);
	//~ problemZ.solve();
}
