
//#include <Species/Species.h>
//#include <Species/LinearViscoelasticSpecies.h>
#include <Species/LinearViscoelasticFrictionSpecies.h>
#include <Mercury3D.h>
#include <Boundaries/PeriodicBoundary.h>
#include "Boundaries/LeesEdwardsBoundary.h"

class restart_simpleshear: public Mercury3D{
public:
    restart_simpleshear()
    {
        setName("p1");
        
        readRestartFile();
        setRestarted(false);
        particleSpecies = dynamic_cast<LinearViscoelasticFrictionSpecies*>(speciesHandler.getObject(0));

		std::cout << "input = " << getName() << std::endl;
		
    }
void setupInitialConditions()
    {
		
        Mdouble velocity =  getYMax();
        		
        boundaryHandler.removeObject(1);
		
        //! assign velocity to Lees Edward y-boundary
        LeesEdwardsBoundary leesEdwardsBoundary;
        leesEdwardsBoundary.set(
            [velocity] (double time) { return time*velocity; },
            [velocity] (double time UNUSED) { return velocity; },
            getXMin(),getXMax(),getYMin(),getYMax());
        boundaryHandler.copyAndAddObject(leesEdwardsBoundary);
		
        //! particleSpecies
        particleSpecies->setDensity(rhop);
        particleSpecies->setStiffnessAndRestitutionCoefficient(K,en,particleSpecies->getDensity()*constants::pi/6.0);
		particleSpecies->setSlidingStiffness(2.0/7.0*particleSpecies->getStiffness());
		particleSpecies->setRollingStiffness(2.0/5.0*particleSpecies->getStiffness());
		particleSpecies->setSlidingFrictionCoefficient(mu_slid);
		particleSpecies->setRollingFrictionCoefficient(mu_roll);
		
		//! particles properties and initial positions
		double N = particleHandler.getNumberOfObjects();
        for (int i=0; i < N; i++) {
			particleHandler.getObject(i)->setSpecies(particleSpecies);
		}
		
		double Rmin = particleHandler.getObject(0)->getRadius();
		double Rmax = particleHandler.getObject(0)->getRadius();
		
		for (int i=1; i < (N-1); i++) {
			Rmin = std::min(Rmin,particleHandler.getObject(i)->getRadius());
			Rmax = std::max(Rmax,particleHandler.getObject(i)->getRadius());
		}
			
		double particleDiameter = 2.0*Rmin;
        double mass = rhop*constants::pi*mathsFunc::cubic(particleDiameter)/6.0;
        double tc = std::sqrt( mass/2.0/K*( mathsFunc::square(constants::pi) + mathsFunc::square(log(en) ) ));

        setTimeStep(tc/50);
		setTimeMax(tmax);
		
		Mdouble Vp = 0;
		for (int i=0; i < N; i++) {
			Vp = Vp + particleHandler.getObject(i)->getVolume();
		}
		
		std::cout << "N = " << N << std::endl;
		
        std::cout << "Lx = " << getXMax() << ", Ly = " << getYMax() << ", Lz = " << getZMax() << std::endl;
        std::cout << "nu = " << Vp/(getXMax()*getYMax()*getZMax()) << std::endl;
        std::cout << "k = " << K << std::endl;
        std::cout << "output = " << getName() << std::endl;
        
        //std::cout << "delta t = " << getTimeStep() << std::endl;
        std::cout << "saving time = " <<  dataFile.getSaveCount()*getTimeStep() << std::endl;
		
    }
    
    LinearViscoelasticFrictionSpecies* particleSpecies;
    
    Mdouble tmax;
    Mdouble rhop, en, K, mu_slid, mu_roll;
};

int main(int argc UNUSED, char *argv[] UNUSED)
{
	restart_simpleshear problem;
//! never change! --------------------------------------------------
    problem.rhop = 1.0;
    problem.mu_slid = 0.0;
    problem.mu_roll = 0.0;
    problem.tmax = 10000.0;
    //! ----------------------------------------------------------------
    
//! assign particles properties: coef. of restitution, stiffness   

    problem.K = 10000;
    problem.en = 0.7;
    
    problem.setName("re_a1");
    
    problem.dataFile.setFileType(FileType::MULTIPLE_FILES_PADDED);
    problem.restartFile.setFileType(FileType::ONE_FILE);
    problem.fStatFile.setFileType(FileType::MULTIPLE_FILES_PADDED);
    problem.eneFile.setFileType(FileType::ONE_FILE);
 
    problem.solve();
}
