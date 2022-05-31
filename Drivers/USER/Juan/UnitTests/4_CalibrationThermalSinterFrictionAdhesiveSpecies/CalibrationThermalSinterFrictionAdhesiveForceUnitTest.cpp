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

#include "DPMBase.h"
#include <iostream>
#include <vector>
#include <Species/ThermalSinterLinFrictionReversibleAdhesiveSpecies.h>
#include <Logger.h>


///This code tests our plastic force model, as published in Luding 2008.
class calibrationThermalSinterFrictionAdhesiveForceUnitTest : public DPMBase{
public:

	void setupInitialConditions() override {
		setXMax(2.0*radius);
		setYMax(radius);
		setZMax(radius);
		setXMin(-getXMax());
		setYMin(-getYMax());
		setZMin(-getZMax());
		setGravity(Vec3D(0.0,0.0,0.0));
		
        Mdouble deltaStar = species->getPenetrationDepthMax() * radius;
        Mdouble dk = (species->getUnloadingStiffnessMax()/species->getLoadingStiffness() - 1.0)/deltaStar;
        Mdouble deltaMax = deltaStar + 1.0/dk;
		Mdouble mass=species->getDensity()*constants::pi/6.0*mathsFunc::cubic(2.0*radius);
        Mdouble relVelocity=sqrt(species->getLoadingStiffness()*mathsFunc::square(chi*deltaMax)*2.0/mass);
        
        particleHandler.clear();

        ThermalParticle P0,P1;
        P0.setSpecies(speciesHandler.getObject(0));
        P1.setSpecies(speciesHandler.getObject(0));
		P0.setPosition(Vec3D(-radius,0,0));
		P0.setVelocity(Vec3D(relVelocity/2.0,0.0,0.0));
		P0.setRadius(radius);
		P1.setPosition(Vec3D( radius,0,0));
		P1.setVelocity(Vec3D(-relVelocity/2.0,0.0,0.0));
		P1.setRadius(radius);
		
		particleHandler.copyAndAddObject(P0);
		particleHandler.copyAndAddObject(P1);
	}

    void set_chi(double new_){chi=new_;}
    double get_chi() {return chi;}

	double chi = 1.0;
    Mdouble radius = 0.5;
    ThermalSinterLinFrictionReversibleAdhesiveSpecies* species = speciesHandler.copyAndAddObject(ThermalSinterLinFrictionReversibleAdhesiveSpecies());
};

int main(int argc UNUSED, char *argv[] UNUSED)
{
    calibrationThermalSinterFrictionAdhesiveForceUnitTest sf;

    //define properties of particle-particle contacts
    sf.species->setDensity(1.04);
    sf.radius = 0.002;

    Mdouble mass=sf.species->getDensity()*((4.0/3.0)*constants::pi)*mathsFunc::cubic(sf.radius);

    double k2 = 0.03; //5000000.0

    //----------------------------------------------
    sf.species->setPlasticParameters(k2, 5.0*k2, 0.5*k2, 0.05);
    //----------------------------------------------
    sf.species->setDissipation(0.8e-5);

    //----------------------------------------------

    //----------------------------------------------
    Mdouble friction = 0.55;

    Mdouble k_s = 0.2*sf.species->getLoadingStiffness();
    Mdouble k_r = 0.1*sf.species->getLoadingStiffness();
    Mdouble k_o = 0.1*sf.species->getLoadingStiffness();

    Mdouble miu_s = 1.0*friction;
    Mdouble miu_r = 0.1*friction;
    Mdouble miu_o = 0.1*friction;

    Mdouble gamma_s = 0.2*sf.species->getDissipation();
    Mdouble gamma_r = 0.05*sf.species->getDissipation();
    Mdouble gamma_o = 0.05*sf.species->getDissipation();
    //-------->[Set Tangential Contribution]
    sf.species->setSlidingStiffness(k_s);
    sf.species->setRollingStiffness(k_r);
    sf.species->setTorsionStiffness(k_o);

    sf.species->setSlidingFrictionCoefficient(miu_s);
    sf.species->setRollingFrictionCoefficient(miu_r);
    sf.species->setTorsionFrictionCoefficient(miu_o);

    sf.species->setSlidingDissipation(gamma_s);
    sf.species->setRollingDissipation(gamma_r);
    sf.species->setTorsionDissipation(gamma_o);

    //----------------------------------------------
    //-------->[Set Adhesive Contribution]
    //Set Adhesive properties:
    Mdouble K_adh = k2;//0.5*k2;//1.0;
    Mdouble f_adh_max = 1e-5*k2;//5.0*(meanRadius*2.0*meanRadius*2.0)*5.0; ///set adhesive force max based on Bond number [Hao]

    sf.species->setAdhesionStiffness(K_adh);
    sf.species->setAdhesionForceMax(f_adh_max);

    //----------------------------------------------
//    sf.species->setHeatCapacity(500);
//    sf.species->setThermalConductivity(20);
//    // fa the sinterAdhesion (due to surface tension, which is dominant over van der Walls)
//    // sr the sinterRate (a temperature-dependent proportionality constant; here we assume that sr=0 below the glass temperature)
//    sf.species->setSinterRate(0.01); //I will include this parameter! Ask for Thoma
    sf.species->setSinterAdhesion(0.013*k2);
//
    sf.species->setSinterType(SINTER_APPROACH::FRENKEL);
//    sf.species->setSinterType(SINTERTYPE::TEMPERATURE_DEPENDENT_FRENKEL);

    /// initial and boundary value for temperature
    double backgroundTemperature = 400;
    /// temperature above which sintering starts
    double glassTemperature = 420;
    /// sintering ends if the temperature drops below this temperature
    double maxTemperatureAtEndOfSinter = backgroundTemperature + 0.5*(glassTemperature-backgroundTemperature);
    /// how much power (in Watt/l^2) to add to the top layer
    /// (note, the number of particles in the top layer, and thus the temperature added, strongly depends on dim)
    double sinterEnergy = 0.2e-3;

    //sf.species->setTemperatureDependentSinterRate([] (double temperature) { return 3.0*(temperature>glassTemperature); });

    //----------------------------------------------

    //----------------------------------------------
    sf.setParticleDimensions(3);
    sf.setSystemDimensions(3);

    sf.setTimeStep(sf.species->computeTimeStep(mass)*2.0);
    sf.setTimeMax(sf.getTimeStep()*400.0);
    sf.setFileType(FileType::ONE_FILE);
    sf.setSaveCount(1);

    logger(INFO,"Testing particle particles collision for elastic plastic forces. \n"
        "This will be done for several values of scaled relative velocity chi");
    
	//Set up constant data that will be used
    const std::vector<double> chi = {0.34, 0.69, 0.8, 0.9, 1.05, 1.25};

    //Loop over all test cases
	for (int i=0; i<4; i++)
    {
        logger(INFO, "Running for chi=%",chi[i]);
		sf.set_chi(chi[i]);
		std::stringstream ss("");
		ss << "CalibrationThermalSinterFrictionAdhesiveForceUnitTest" << sf.get_chi();
		sf.setName(ss.str().c_str());
		sf.solve();
        sf.writeRestartFile();
    }

    for (int i=0; i<4; i++)
    {
        logger(INFO, "Running for chi=%",chi[i]);
        sf.set_chi(chi[i]);
        std::stringstream ss2("");
        ss2 << "NeckGrowthTest" << sf.get_chi();
        sf.setName(ss2.str().c_str());
        sf.solve();
        sf.writeRestartFile();
    }
	//Neck
//    sf.set_chi(0.6);
//    sf.setName("NeckGrowthTest");
//    sf.setSaveCount(100);
//    sf.setTimeMax(5e-9);
//    sf.solve();
//    sf.writeRestartFile();
    
    std::cout << "Execute 'gnuplot CalibrationThermalSinterFrictionAdhesiveForceUnitTest.gnu' to view output" << std::endl;

    helpers::writeToFile("SinterForceUnitTest.gnu",
                         "set xlabel 'displacement'\n"
                         "set ylabel 'force'\n"
                         "plot 'CalibrationThermalSinterFrictionForceAdhesiveUnitTest1.05.fstat' u 7:9 w lp\n"
                         );

    std::cout << "Execute 'gnuplot NeckGrowthTest.gnu' to view output" << std::endl;

    helpers::writeToFile("NeckGrowthTest.gnu",
                         "set xlabel 'time [s]'\n"
                         "set ylabel 'x/a'\n"
                         "plot 'NeckGrowthTest.fstat' u ($1):(sqrt($7/0.002)) w lp lt rgb 'royalblue'\n"
    );
}
