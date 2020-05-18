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

template <StatType T> class statistics_while_running : public StatisticsVector<T>, public Chute
{
	public:
	statistics_while_running<T>() : StatisticsVector<T>() , Chute() {}

	void actionsBeforeTimeStep(){}
		
	void setupInitialConditions() {StatisticsVector<T>::printStat();
		//~ for (unsigned int i=0; i<particleHandler.getNumberOfObjects(); i++) {
			//~ if (particleHandler.getObject(i)->isFixed() || i>10) {removeParticle(i); i--;}
		//~ }
	}
	
	void printTime() const {
		static Mdouble tint = getTimeMax()-getTime();
		std::cout << "\r" << std::setprecision(2) << std::setw(5) << static_cast<int>(100.*(1-(getTimeMax()-getTime())/tint)) << "%\r";
		std::cout.flush();
	}

	// allows getZMax() to be set to the height of the highest particle
	void auto_set_domain() {
		if (particleHandler.getNumberOfObjects()) {
		  Vec3D Max= particleHandler.getObject(0)->getRadius()*Vec3D(1,1,1) + particleHandler.getObject(0)->getPosition();
		  Vec3D Min= - particleHandler.getObject(0)->getRadius()*Vec3D(1,1,1)  + particleHandler.getObject(0)->getPosition();//Particles[0].Position-Particles[0].Radius;
			for (unsigned int i=0; i<particleHandler.getNumberOfObjects(); i++) {
			  Max=Vec3D::max(Max,particleHandler.getObject(i)->getRadius()*Vec3D(1,1,1)  + particleHandler.getObject(0)->getPosition());
			  Min=Vec3D::min(Min,- particleHandler.getObject(i)->getRadius()*Vec3D(1,1,1)  + particleHandler.getObject(0)->getPosition());
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
	std::cout << "restart data: " << name << ".restart" << std::endl;
	
	//load restart data
	std::cout << "XY averaging" << std::endl;
 	statistics_while_running<XY> problem;
 	problem.setName(name.c_str());
	problem.readRestartFile();//load_restart_data();
	problem.setFixedParticleRadius(problem.particleHandler.getObject(0)->getRadius());
	problem.setInflowParticleRadius(problem.particleHandler.getObject(0)->getRadius());
	problem.restartFile.setFileType(FileType::ONE_FILE);//restartFile.setFileType(FileType::ONE_FILE);
	problem.writeRestartFile();
	problem.auto_set_domain();
	//keep file name but create files in the local directory, i.e. remove folder
	std::cout << "old name: " << problem.getName() << std::endl;
	size_t found=name.find_last_of("/\\");
	problem.setName(name.substr(found+1).c_str());
	std::cout << "new name: " << problem.getName() << std::endl;
	//set output to minimum
        problem.dataFile.setFileType(FileType::NO_FILE);
        problem.restartFile.setFileType(FileType::NO_FILE);
        problem.fStatFile.setFileType(FileType::NO_FILE);
        problem.eneFile.setFileType(FileType::ONE_FILE);
	//set statistical parameters
        problem.setDoPeriodicWalls(false);
 	//~ problem.setZMinStat(-1);
 	problem.setN(50);
 	problem.setCGWidth(.1);
	problem.setSaveCount(25);
	problem.setStressTypeForFixedParticles(3);
 	problem.setCGTimeMin(problem.getTime());
	problem.setTimeMax(problem.getTime()+tint);
	//solve and create live statistics
	problem.readStatArguments(argc, argv);
	//~ problem.Chute::write(cout, true);
	//~ cout << endl << problem.StatisticsVector<XY>::write() << endl;
	problem.solve();

	std::cout << std::endl << "Y averaging" << std::endl;
 	statistics_while_running<Y> problemY;
 	problemY.setName(name.c_str());
	problemY.readRestartFile();//load_restart_data();
	problemY.setFixedParticleRadius(problemY.particleHandler.getObject(0)->getRadius());
	problemY.setInflowParticleRadius(problemY.particleHandler.getObject(0)->getRadius());
        problemY.restartFile.setFileType(FileType::ONE_FILE);
	problemY.writeRestartFile();
	problemY.auto_set_domain();
	//keep file name but create files in the local directory, i.e. remove folder
	std::cout << "old name: " << problemY.getName() << std::endl;
	std::stringstream ss;
	ss << name.substr(found+1).c_str() << "Y";
	problemY.setName(ss.str().c_str());
	std::cout << "new name: " << problemY.getName() << std::endl;
	//set output to minimum
        problemY.dataFile.setFileType(FileType::NO_FILE);
        problemY.restartFile.setFileType(FileType::ONE_FILE);
        problemY.fStatFile.setFileType(FileType::NO_FILE);
        problemY.eneFile.setFileType(FileType::ONE_FILE);
	//set statistical parameters
        problemY.setDoPeriodicWalls(false);
 	//~ problemY.setZMinStat(-1);
 	problemY.setN(50);
 	problemY.setCGWidth(.1);
	problemY.setSaveCount(25);
	problemY.setStressTypeForFixedParticles(3);
 	problemY.setCGTimeMin(problemY.getTime());
	problemY.setTimeMax(problemY.getTime()+tint);
	//solve and create live statistics
	problemY.readStatArguments(argc, argv);
	//~ problemY.Chute::write(cout, true);
	//~ cout << endl << problemY.StatisticsVector<Y>::write() << endl;
	problemY.solve();


	std::cout << std::endl << "X averaging" << std::endl;
 	statistics_while_running<X> problemX;
 	problemX.setName(name.c_str());
	problemX.readRestartFile();//load_restart_data();
	problemX.setFixedParticleRadius(problemX.particleHandler.getObject(0)->getRadius());
	problemX.setInflowParticleRadius(problemX.particleHandler.getObject(0)->getRadius());
	problemX.restartFile.setFileType(FileType::ONE_FILE);
	problemX.writeRestartFile();
	problemX.auto_set_domain();
	//keep file name but create files in the local directory, i.e. remove folder
	std::cout << "old name: " << problemX.getName() << std::endl;
	ss << name.substr(found+1).c_str() << "X";
	problemX.setName(ss.str().c_str());
	std::cout << "new name: " << problemX.getName() << std::endl;
	//set output to minimum
        problemX.dataFile.setFileType(FileType::NO_FILE);
        problemX.restartFile.setFileType(FileType::ONE_FILE);
        problemX.fStatFile.setFileType(FileType::NO_FILE);
        problemX.eneFile.setFileType(FileType::ONE_FILE);

	  //problemX.dataFile.setFileType(FileType::NO_FILE);
	  //problemX.restartFile.setFileType(FileType::ONE_FILE);
	  //problemX.fStatFile.setFileType(FileType::NO_FILE);
	  //problemX.eneFile.setFileType(FileType::ONE_FILE);
	//set statistical parameters
        problemX.setDoPeriodicWalls(false);
 	//~ problemX.setZMinStat(-1);
 	problemX.setN(50);
 	problemX.setCGWidth(.1);
	problemX.setSaveCount(25);
	problemX.setStressTypeForFixedParticles(3);
 	problemX.setCGTimeMin(problemX.getTime());
	problemX.setTimeMax(problemX.getTime()+tint);
	//solve and create live statistics
	problemX.readStatArguments(argc, argv);
	//~ problemX.Chute::write(cout, true);
	//~ cout << endl << problemX.StatisticsVector<X>::write() << endl;
	problemX.solve();

}
