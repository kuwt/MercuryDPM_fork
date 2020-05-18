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

template <StatType T> class CLiveStatistics : public StatisticsVector<T>, public Chute
{
	public:
	CLiveStatistics<T>() : StatisticsVector<T>() , Chute() {
		
		}

	CLiveStatistics<T>(unsigned int argc, char *argv[]) : StatisticsVector<T>(argc, argv) , Chute() {
		
		}

	void actionsBeforeTimeStep() override {};
		
	void setupInitialConditions() override {StatisticsVector<T>::printStat();
		write(std::cout,false);
	};
	
	void printTime() const override {
		static double tint = getTimeMax()-getTime();
		std::cout << "\r" << std::setprecision(2) << (int)100.*(1-(getTimeMax()-getTime())/tint) << "%\r";
		std::cout.flush();
	}
	
	void getLiveStatistics() {
	  readRestartFile();//load_restart_data();
		//set output to minimum
          dataFile.setFileType(FileType::NO_FILE);
          restartFile.setFileType(FileType::ONE_FILE);
          fStatFile.setFileType(FileType::NO_FILE);
	  //		dataFile.setFileType(FileType::NO_FILE);
	  //	restartFile.setFileType(FileType::ONE_FILE);
	  //	fStatFile.setFileType(FileType::NO_FILE);
		//set statistical parameters
		//~ setDoPeriodicWalls(false);
		//~ setZMinStat(-1);
		//~ setN(200);
		//~ setCGWidth(.1);
		setSaveCount(5);
		//~ set_infiniteStressForFixedParticles(true);
		//~ setTimeMax(problem.getTime()+.1);
		//solve and create live statistics
		solve();
	}
	
};


int main(int argc, char *argv[]) {	
	
  if (argc>1&&strcmp(argv[1],"-help") != 0) std::cout << std::endl << "Get statistics for " << argv[1] << std::endl;

	//check for '-stattype' option
	StatType T = XYZ;
	for (int i = 2; i<argc; i++) {
		if (!strcmp(argv[i],"-stattype")) {
			if (!strcmp(argv[i+1],"XYZ")) T = XYZ;
			else if (!strcmp(argv[i+1],"XY")) T = XY;
			else if (!strcmp(argv[i+1],"XZ")) T = XZ;
			else if (!strcmp(argv[i+1],"YZ")) T = YZ;
			else if (!strcmp(argv[i+1],"X")) T = X;
			else if (!strcmp(argv[i+1],"Y")) T = Y;
			else if (!strcmp(argv[i+1],"Z")) T = Z;
			else {std::cerr << "stattype unknown" << std::endl; exit(-1);}
		}
	}
	if (T==XY) { // averaging in z-direction
	  std::cout << "averaging in z-direction" << std::endl;
		CLiveStatistics<XY> stats(argc, argv);
		stats.setDoPeriodicWalls(false);
		stats.getLiveStatistics();
	} else if (T==XZ) { // averaging in x-direction
	  std::cout << "averaging in y-direction" << std::endl;
		CLiveStatistics<XZ> stats(argc, argv);
		stats.setDoPeriodicWalls(false);
		stats.getLiveStatistics();
	} else if (T==YZ) { // averaging in x-direction
	  std::cout << "averaging in x-direction" << std::endl;
		CLiveStatistics<YZ> stats(argc, argv);
		stats.setDoPeriodicWalls(false);
		stats.getLiveStatistics();
	} else if (T==X) { // averaging in yz-direction
	  std::cout << "averaging in yz-direction" << std::endl;
		CLiveStatistics<X> stats(argc, argv);
		stats.setDoPeriodicWalls(false);
		stats.getLiveStatistics();
	} else if (T==Y) { // averaging in yz-direction
	  std::cout << "averaging in xz-direction" << std::endl;
		CLiveStatistics<Y> stats(argc, argv);
		stats.setDoPeriodicWalls(false);
		stats.getLiveStatistics();
	} else if (T==Z) { // averaging in yz-direction
	  std::cout << "averaging in xy-direction" << std::endl;
		CLiveStatistics<Z> stats(argc, argv);
		stats.setDoPeriodicWalls(false);
		stats.getLiveStatistics();
	} else { //default: no averaging
		CLiveStatistics<XYZ> stats(argc, argv);
		stats.setDoPeriodicWalls(false);
		stats.getLiveStatistics();
	}
	return 0;
}
