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

#include "scr/Chute.h"
#include "scr/StatisticsVector.h"
using namespace std;

template <StatType T> 
class statistics_while_running : public StatisticsVector<T>, public Chute
{
	public:
	statistics_while_running<T>() : StatisticsVector<T>() , Chute() {}

	void run (int argc, char *argv[]) {
		//select file name to restart from
		string name;
		if (argc>1) {name=argv[1]; argv++; argc--;} else { cout << "please provide a name" << endl; exit(-1); }
		//select time interval
		Mdouble tint=.2; 
		if (argc>1) {tint=atof(argv[1]); argv++; argc--; cout << "tint=" << tint << endl;}
		argv--; argc++;
		//create class and set name
		setName(name.c_str());
		//load restart data
		int i;
		for (i = 1; i<argc; i++) if (!strcmp(argv[i],"-restartfilename")) break;
		if (i<argc) {
			string filename(argv[i+1]);
			cout << "loading restart data:: " << filename << endl;
			load_restart_data(filename);
		} else {
			cout << "loading restart data: " << name << ".restart" << endl;
			load_restart_data();
		}
		//save restart data
		//restartFile.setFileType(FileType::ONE_FILE);
		//writeRestartFile();
		//automatic depth range detection
		auto_set_z();
		//keep file name but create files in the local directory, i.e. remove folder
		cout << "changed name from " << getName();
		size_t found=name.find_last_of("/\\");
		setName(name.substr(found+1).c_str());
		cout << " to " << getName() << endl;
		//save restart data
		//restartFile.setFileType(FileType::ONE_FILE);
		//writeRestartFile();
		//set output to minimum
		dataFile.setFileType(FileType::NO_FILE);
		restartFile.setFileType(FileType::NO_FILE);
		fStatFile.setFileType(FileType::NO_FILE);
		eneFile.setFileType(FileType::NO_FILE);
		//automatic depth range decection
		auto_set_z();
		//set statistical parameters
		this->setDoPeriodicWalls(false);
		//~ setZMinStat(-1);
		this->setN(200);
		setCGWidth(getInflowParticleRadius());
		cout << "by default w=InflowParticleRadius=" << this->getCGWidth() << endl;
		setSaveCount(25);
		//~ setStressTypeForFixedParticles(3); //default 1
		setCGTimeMin(getTime());
		setTimeMax(getTime()+tint);
		//solve and create live statistics
		cout << "read arguments" << endl;
		this->readStatArguments(argc, argv);
		solve();
	}


	void actionsBeforeTimeStep(){};
		
	void setupInitialConditions() {
		write(std::cout,false);
		StatisticsVector<T>::write();
		//~ for (unsigned int i=0; i<particleHandler.getNumberOfObjects(); i++) {
			//~ if (particleHandler.getObject(i)->isFixed() || i>10) {removeParticle(i); i--;}
		//~ }
	};
	
	void printTime() {
		if (this->getVerbosityLevel()<1) return;
		static Mdouble tint = getTimeMax()-getTime();
		cout << "\r" << setprecision(2) << setw(5) << (int)100.*(1-(getTimeMax()-getTime())/tint) << "%\r";
		cout.flush();
	}

	// allows zmax to be set to the height of the highest particle
	void auto_set_z() {
		if (particleHandler.getNumberOfObjects()) {
			Mdouble zmin=particleHandler.getObject(0)->->getPosition().Z+particleHandler.getObject(0)->->getRadius();
			Mdouble zmax=particleHandler.getObject(0)->->getPosition().Z-particleHandler.getObject(0)->->getRadius();
			for (unsigned int i=0; i<particleHandler.getNumberOfObjects(); i++) {
				if (particleHandler.getObject(i)->getPosition().Z+particleHandler.getObject(i)->getRadius()>zmax) zmax=particleHandler.getObject(i)->getPosition().Z+particleHandler.getObject(i)->getRadius();
				if (particleHandler.getObject(i)->getPosition().Z-particleHandler.getObject(i)->getRadius()<zmin) zmin=particleHandler.getObject(i)->getPosition().Z-particleHandler.getObject(i)->getRadius();
			}
			setZMax(zmax+.05*(zmax-zmin));
			setZMin(zmin-.05*(zmax-zmin));
		}
	}

	Mdouble TimeAverageSurface;
	Mdouble TimeAverageSurface2;
	Mdouble TimeAverageBase;
	Mdouble TimeAverageBase2;
	int nTimeAvg;

	void initialiseStatistics()
	{
		StatisticsVector<T>::initialiseStatistics();
		TimeAverageSurface = 0.0;
		TimeAverageSurface2 = 0.0;
		TimeAverageBase = 0.0;
		TimeAverageBase2 = 0.0;
		nTimeAvg = 0;
	}

	void processStatistics()
	{
		StatisticsVector<T>::processStatistics();
		Mdouble Base = 1e20, Surface=-1e20;
		for (unsigned int i=0; i<particleHandler.getNumberOfObjects(); i++) if (!particleHandler.getObject(i)->isFixed()) {
			Base=min(Base,particleHandler.getObject(i)->getPosition().Z);
			Surface=max(Surface,particleHandler.getObject(i)->getPosition().Z);
		}
		TimeAverageSurface += Surface;
		TimeAverageSurface2 += mathsFunc::square(Surface);
		TimeAverageBase += Base;
		TimeAverageBase2 += mathsFunc::square(Base);
		nTimeAvg++;
		//~ cout << "Height " << setprecision(14) << Surface << "_" << Base << endl;
	}

	void finishStatistics()
	{
		StatisticsVector<T>::finishStatistics();
		TimeAverageSurface /= nTimeAvg;
		TimeAverageSurface2 /= nTimeAvg;
		TimeAverageBase /= nTimeAvg;
		TimeAverageBase2 /= nTimeAvg;
		stringstream filename("");
		filename << getName() << ".height";
		stringstream text("");
		text << TimeAverageBase - 0.5 - sqrt(TimeAverageBase2-mathsFunc::square(TimeAverageBase))
			<< " " << TimeAverageSurface + 0.5 + sqrt(TimeAverageSurface2-mathsFunc::square(TimeAverageSurface));
		writeToFile(filename.str(),text.str());
	}

	bool writeToFile(string filename, string text){
		fstream file;
		file.open(filename.c_str(),fstream::out);
		file << text;
		file.close();
		return true;
	}
	
	bool appendToFile(string filename, string text){
		fstream file;
		file.open(filename.c_str(),fstream::out | fstream::app);
		file << text;
		file.close();
		return true;
	}
	
};

int main(int argc, char *argv[])
{
	//check for '-stattype' option (default Z)
	StatType T = Z;
	for (int i = 2; i<argc; i++) {
		if (!strcmp(argv[i],"-stattype")) {
			if (!strcmp(argv[i+1],"XYZ")) T = XYZ;
			else if (!strcmp(argv[i+1],"XY")) T = XY;
			else if (!strcmp(argv[i+1],"XZ")) T = XZ;
			else if (!strcmp(argv[i+1],"YZ")) T = YZ;
			else if (!strcmp(argv[i+1],"X")) T = X;
			else if (!strcmp(argv[i+1],"Y")) T = Y;
			else if (!strcmp(argv[i+1],"Z")) T = Z;
			else if (!strcmp(argv[i+1],"O")) T = O;
			else {cerr << "stattype unknown" << endl; exit(-1);}
		}
	}
	if (T==XY) { // averaging in z-direction
		cout << "averaging in z-direction" << endl;
		statistics_while_running<XY> stats; 
		stats.run(argc, argv);
	} else if (T==XZ) { // averaging in y-direction
		cout << "averaging in y-direction" << endl;
		statistics_while_running<XZ> stats; 
		stats.run(argc, argv);
	} else if (T==YZ) { // averaging in x-direction
		cout << "averaging in x-direction" << endl;
		statistics_while_running<YZ> stats; 
		stats.run(argc, argv);
	} else if (T==X) { // averaging in yz-direction
		cout << "averaging in yz-direction" << endl;
		statistics_while_running<X> stats; 
		stats.run(argc, argv);
	} else if (T==Y) { // averaging in yz-direction
		cout << "averaging in xz-direction" << endl;
		statistics_while_running<Y> stats; 
		stats.run(argc, argv);
	} else if (T==Z) { // averaging in yz-direction
		cout << "averaging in xy-direction" << endl;
		statistics_while_running<Z> stats; 
		stats.run(argc, argv);
	} else { //default: no averaging
		cout << "no averaging" << endl;
		statistics_while_running<XYZ> stats; 
		stats.run(argc, argv);
	}
}
