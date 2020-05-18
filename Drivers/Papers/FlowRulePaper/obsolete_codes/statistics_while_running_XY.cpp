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

template <StatType T> class statistics_while_running : public StatisticsVector<T>, public Chute
{
public:
	statistics_while_running<T>() : StatisticsVector<T>() , Chute() {}

	void actionsBeforeTimeStep(){};
		
	void actionsBeforeTimeLoop() {
		fix_hgrid();
		write(std::cout,false);
		cout<< endl
			<< StatisticsVector<T>::write() 
			<< endl;
		//~ for (unsigned int i=0; i<particleHandler.getNumberOfObjects(); i++) {
			//~ if (particleHandler.getObject(i)->isFixed() || i>10) {removeParticle(i); i--;}
		//~ }
	};
	
	void printTime() {
		static Mdouble tint = getTimeMax()-getTime();
		cout << "\r" 
		<< setprecision(7) << setw(10) << getTime() << ", dt" << tint << ", "
		<< setprecision(2) << setw(5) << (int)100.*(1-(getTimeMax()-getTime())/tint) 
		<< "%\r"
		;
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

	
};

int main(int argc, char *argv[])
{
	//select file name to restart from
	string name;
	if (argc>1) {name=argv[1]; argv++; argc--;} else name="../ini_statistics/H10A22L1M0.5B0.5";
	//select time interval
	Mdouble tint=.2; 
	if (argc>1) {tint=atof(argv[1]); argv++; argc--; cout << "tint=" << tint << endl;}
	argv--; argc++;
	//load restart data
	cout << "restart data: " << name << ".restart" << endl;
 	statistics_while_running<XY> problem;
 	problem.setName(name.c_str());
	problem.load_restart_data();
	//problem.setRestarted(false);
	problem.restartFile.setFileType(FileType::ONE_FILE);
	problem.writeRestartFile();
	//problem.auto_set_z();
	//keep file name but create files in the local directory, i.e. remove folder
	cout << "old name: " << problem.getName() << endl;
	size_t found=name.find_last_of("/\\");
	problem.setName(name.substr(found+1).c_str());
	cout << "new name: " << problem.getName() << endl;
	//set output to minimum
	problem.dataFile.setFileType(FileType::NO_FILE);
	problem.restartFile.setFileType(FileType::NO_FILE);
	problem.fStatFile.setFileType(FileType::NO_FILE);
	problem.eneFile.setFileType(FileType::ONE_FILE);
	//set statistical parameters
    problem.setDoPeriodicWalls(false);
 	//~ problem.setZMinStat(-1);
 	problem.setNX(50);
 	problem.setNY(10);
 	problem.setCGWidth(.1);
	problem.setSaveCount(25);
	//problem.set_infiniteStressForFixedParticles(true);
 	problem.setCGTimeMin(problem.getTime());
	problem.setTimeMax(problem.getTime()+tint);
	//solve and create live statistics
	problem.readStatArguments(argc, argv);
	problem.solve();
}
