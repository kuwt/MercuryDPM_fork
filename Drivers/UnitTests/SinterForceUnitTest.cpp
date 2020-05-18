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

#include "DPMBase.h"
#include <iostream>
#include <vector>
//! [T11:contactModel]
#include <Species/SinterSpecies.h>
//! [T11:contactModel]
#include <Logger.h>


///This code tests our plastic force model, as published in Luding 2008.
class SinterForceUnitTest : public DPMBase{
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
		
		SphericalParticle P0,P1;
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
    SinterSpecies* species = speciesHandler.copyAndAddObject(SinterSpecies());
};

int main(int argc UNUSED, char *argv[] UNUSED)
{	
	SinterForceUnitTest sf;
	double k1=5.0;
    sf.radius = 2e-6;
    sf.species->setDensity(6.0/constants::pi);
	sf.species->setPlasticParameters(k1, 5.0*k1, k1, 0.5);
    sf.setParticleDimensions(3);
    sf.setSystemDimensions(3);
    sf.species->setSinterRate(1000);

    Mdouble mass=sf.species->getDensity()*constants::pi/6.0*mathsFunc::cubic(2.0*sf.radius);
    sf.species->setDissipation(5e5*mass);
    sf.setTimeStep(sf.species->computeTimeStep(mass)*2.0);
    sf.setTimeMax(sf.getTimeStep()*400.0);
    sf.setFileType(FileType::ONE_FILE);
    sf.setSaveCount(1);

    logger(INFO,"Testing particle particles collision for elastic plastic forces. \n"
        "This will be done for several values of scaled relative velocity chi");
    
	//Set up constant data that will be used
    const std::vector<double> chi = {0.34, 0.69, 1.05, 1.25};

    //Loop over all test cases
	for (int i=0; i<4; i++)
    {
        logger(INFO, "Running for chi=%",chi[i]);
		sf.set_chi(chi[i]);
		std::stringstream ss("");
		ss << "SinterForceUnitTest" << sf.get_chi();
		sf.setName(ss.str().c_str());
		sf.solve();
        sf.writeRestartFile();
    }

    //A longer simulation, with sintering activated
    sf.set_chi(1.15);
    sf.setName("LongSinterForceUnitTest");
    sf.setSaveCount(2500);
    sf.setTimeMax(5e-4);
    sf.solve();
    sf.writeRestartFile();
    
    std::cout << "Execute 'gnuplot SinterForceUnitTest.gnu' to view output" << std::endl;
    helpers::writeToFile("SinterForceUnitTest.gnu",
                         "set xlabel 'displacement'\n"
                         "set ylabel 'force'\n"
                         "plot 'SinterForceUnitTest1.05.fstat' u 7:9 w lp\n"
                         );
    
    std::cout << "Execute 'gnuplot LongSinterForceUnitTest.gnu' to view output" << std::endl;
    helpers::writeToFile("LongSinterForceUnitTest.gnu",
                         "sinterRate = 1000\n"
                         "radius = 2e-6\n"
                         "set xlabel 'time'\n"
                         "set ylabel 'displacement'\n"
                         "plot 'LongSinterForceUnitTest.fstat' u 1:(sqrt(2.0*$7/radius)), sqrt((sinterRate*x)+0.5**2)\n"
                         //"plot 'LongSinterForceUnitTest.fstat' u 1:(sqrt(2.0*$7/radius)), sqrt((2.0/radius*sinterRate*x)+0.77**2)\n"
                         //"plot 'LongSinterForceUnitTest.fstat' u 1:(sqrt(2.0*$7/radius)), sqrt(sinterRate*x)\n"
                         );
}
