/*
		FIRST ATTEMPT TO A PBM CODE FOR SCREW FEEDING
The initial aim is to correctly model the particle wear and breakage durinf the feedign process.
Mechanistic insights should be added in the breaking chance function and in the breaking pattern.

An ensemble of particles particelEnsemble composed of particels of class particle is created.
Whenever a particle of teh ensemble breaks, a new one is pushed back into the ensemble array, and teh original one modified.
Volume conservation is imposed in this transition to ensure mass conservation.
*/

#include <Mercury3D.h>
#include <math.h>
#include <fstream>

class particle : public Mercury3D
{
	public:
	// CONSTRUCTOR
	//particle(double r) : radius(r) {}
	particle(int i, double r)
	{
		particleLabel = i;
		radius = r;
		volume = 4.0*constants::pi*pow(radius,3.0)/3.0;
		lastInteractionTime = 0.0;
	}

	// SETTERS
	void setRadius(double r)
	{
		radius = r;
		volume = 4.0*constants::pi*pow(radius,3.0)/3.0;
	}
	void setVolume(double v)
	{
		volume = v;
		radius = pow(3.0*volume/(4.0*constants::pi),1.0/3.0);
	}
	void setLastInteractionTime(double t) {lastInteractionTime = t;}

	// GETTERS
	int getIndex() {return particleLabel;}
	double getRadius() {return radius;}
	double getVolume() {return volume;}
	double getLastInteractionTime() {return lastInteractionTime;}

	// PRINTERS
	void printRadius() {std::cout << radius << "\n";}
	void printStats() {std::cout << "# = " << particleLabel << ", r = " << radius << ", V = " << volume << ", t0 = " << lastInteractionTime << std::endl;}

	// VARIABLES
	private:
	int particleLabel;
	double radius;
	double volume;
	double lastInteractionTime;
};

class particleEnsemble : public Mercury3D
{
	public:
	// CONSTRUCTOR
	//particleEnsemble(int nP, double radius) : particles(nP, particle(nP, radius)) {std::cout << "Number of particles: " << particles.size() << std::endl;}
	particleEnsemble(int nP, double rMean, double rVar)
	{
		particles.reserve(nP);
		nParticles = nP;
		for (int i = 0; i < nP; i++) particles.push_back(particle(i, rMean*(1.0 + rVar*random.getRandomNumber(-1.0,1.0))));
	}

	// SETTERS
	void setMinBreakageVolume(double vMin) {minBreakageVolume = vMin;}

	void setCharacteristicBreakageTime(double t) {characteristicTime = t;}
	
	void setTimeEvolutionParameters(double tMax, double dt)
	{
		timeMax = tMax;
		integrationTime = dt;
	}

	// GETTERS
	double getNumberOfParticles() {return particles.size();}

	double getTotalParticleVolume()
	{
		double vTot = 0.0;
		for (int i=0; i<particles.size(); i++) vTot += particles[i].getVolume();
		return vTot;
	}

	double getAverageReactionTime()
	{
		double tAverage = 0.0;
		for (int i=0; i<particles.size(); i++) tAverage += time - particles[i].getLastInteractionTime();
		return tAverage/nParticles;
	}

	// DOERS
	void evolveEnsemble()
	{
		int numberOfParticlesBeforeEvolution;
		double randomVariable;

		saveInitialConfiguration();
		openOutputFile();

		do
		{
			numberOfParticlesBeforeEvolution = particles.size();
			std::cout << "t = " << time << ", nP = " << nParticles << ", vMean = " << getTotalParticleVolume()/nParticles << ", DtMean = " << getAverageReactionTime() << std::endl;
			writeToOutputFile();

			for (int i=0; i<numberOfParticlesBeforeEvolution; i++)
			{
				randomVariable = random.getRandomNumber(0.0,1.0);
				if (randomVariable > exp(-(time - particles[i].getLastInteractionTime())*probability/characteristicTime))
				{
					breakParticleI(i);
				}
			}

			nParticles = particles.size();
			time += integrationTime;
		} while (time <= timeMax);

		closeOutputFile();
		saveFinalConfiguration();
	}

	void breakParticleI(int i)
	{
		if (particles[i].getVolume() > minBreakageVolume)
		{
			int highestIndexParticlesAtBreakageEvent = particles.back().getIndex();
			double newVolumeI = particles[i].getVolume()*random.getRandomNumber(0.0,1.0);
			//double newRadiusI = pow(3.0*newVolumeI/(4.0*constants::pi),1.0/3.0);

			particles.push_back(particle(highestIndexParticlesAtBreakageEvent + 1, pow(3.0*(particles[i].getVolume() - newVolumeI)/(4.0*constants::pi),1.0/3.0)));
			particles[i].setVolume(newVolumeI);

			particles[i].setLastInteractionTime(time);
			particles.back().setLastInteractionTime(time);

			std::cout << time << ", " << i << " -> " << highestIndexParticlesAtBreakageEvent + 1 << std::endl;
		}
	}

	void saveDataOnce()
	{
		std::ofstream dataExpFile;
		dataExpFile.open("PBM.data", std::ios::out);
		for (int i=0; i<nParticles; i++) dataExpFile << particles[i].getIndex() << " " << particles[i].getRadius() << std::endl;
		dataExpFile.close();
	}

	// PRINTERS
	void printAllStats() {for (int i=0; i<particles.size(); i++) particles[i].printStats();}

	void printStats(int i)
	{
		if (i < particles.size()) particles[i].printStats();
		else std::cout << "ERROR: TRYING TO ACCESS VALUE OUT OF BOUNDS" << std::endl;
	}

	void openOutputFile()
	{
		std::cout << "Opening output file..." << std::endl;
		dataExpFile.open("PBM.data", std::ios::out);
		dataExpFile << "1) time   2) nP   3) vTot   4) vMean   5) DtMean" << std::endl;
	}

	void writeToOutputFile()
	{
		dataExpFile <<
		time << "   " <<
		getNumberOfParticles() << "   " <<
		getTotalParticleVolume() << "   " <<
		getTotalParticleVolume()/getNumberOfParticles() << "   " <<
		getAverageReactionTime() << std::endl;
	}

	void closeOutputFile()
	{
		std::cout << "Closing output file..." << std::endl;
		dataExpFile.close();
	}

	void saveInitialConfiguration()
	{
		std::ofstream initialConfigurationDataFile;
		initialConfigurationDataFile.open("PBM_I.data", std::ios::out);
		for (int i=0; i<nParticles; i++)
		{
			initialConfigurationDataFile << particles[i].getIndex() << "   "
			<< particles[i].getRadius() << "   "
			<< particles[i].getVolume() << "   "
			<< particles[i].getLastInteractionTime() << std::endl;
		}
		initialConfigurationDataFile.close();
	}

	void saveFinalConfiguration()
	{
		std::ofstream finalConfigurationDataFile;
		finalConfigurationDataFile.open("PBM_F.data", std::ios::out);
		for (int i=0; i<nParticles; i++)
		{
			finalConfigurationDataFile << particles[i].getIndex() << "   "
			<< particles[i].getRadius() << "   "
			<< particles[i].getVolume() << "   "
			<< particles[i].getLastInteractionTime() << std::endl;
		}
		finalConfigurationDataFile.close();
	}

	// VARIABLES
	private:
	int nParticles;
	std::vector<particle> particles;
	double minBreakageVolume;
	double characteristicTime;
	double probability = 0.001;
	double time = 0.0;
	double timeMax;
	double integrationTime;
	std::ofstream dataExpFile;
};

int main(int argc, char *argv[])
{
	int numberOfInitialParticles = 10000;
	double initialRadius = 1.0;
	double radiusVariance = 0.5;

	particleEnsemble ensemble(numberOfInitialParticles, initialRadius, radiusVariance);
	ensemble.setMinBreakageVolume(2.0);
	ensemble.setCharacteristicBreakageTime(0.1);
	ensemble.setTimeEvolutionParameters(20.0, 0.01);

	// save data once
	//ensemble.saveDataOnce();

	std::cout << "Ensemble before breakage" << std::endl;
	ensemble.printAllStats();
	ensemble.evolveEnsemble();
	//std::cout << "Ensemble after breakage" << std::endl;
	ensemble.printAllStats();



	//int nParticles = 10;
	//int nReactingParticles = nParticles;
	//double probability = 0.05;
	//double dt = 0.01;
	//double tMax = 10.0;
	//double tCharacteristic = 0.1;

	//double tLastReaction = 0.0;

	//double t = 0.0;

	//std::srand(std::time(nullptr));

	//do
	//{
		//for (int i = 0; i < nReactingParticles; i++)
		//{
			//if (std::rand() < (1.0 - exp(-(t-tLastReaction)*probability/tCharacteristic))*RAND_MAX) nParticles++;
		//}

		//std::cout << "time: " << t << ", nP(t): " << nReactingParticles << ", nP(t+dt): " << nParticles << ", dN: " << nParticles - nReactingParticles << std::endl;

		//if (nParticles != nReactingParticles) tLastReaction = t;
		//nReactingParticles = nParticles;
		//t += dt;
	//}
	//while (t < tMax);


	return 0;
}
