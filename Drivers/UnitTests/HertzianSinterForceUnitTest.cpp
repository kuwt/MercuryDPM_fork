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
#include <Species/HertzianSinterSpecies.h>
#include <Logger.h>


///This code tests our plastic force model, as published in Luding 2008.
class HertzianSinterForceUnitTest : public DPMBase
{
public:

    HertzianSinterForceUnitTest()
    {
        species = speciesHandler.copyAndAddObject(HertzianSinterSpecies());
    }

	void setupInitialConditions() override
	{
		setXMax(3.0);
		setYMax(2.0);
		setZMax(2.0);
		setGravity(Vec3D(0.0,0.0,0.0));
        species->setDensity(6.0/constants::pi);
		
		Mdouble depth = 0.05;
        Mdouble effectiveDiameter = 0.5;
		Mdouble maxFactor = 1.0 - mathsFunc::square(cbrt((species->getLoadingModulus()+species->getCohesionModulus())/species->getUnloadingModulusMax()));
        Mdouble deltaStar = species->getPenetrationDepthMax() * effectiveDiameter / maxFactor;    
        
        Mdouble delMax = chi*deltaStar;
        Mdouble contactRadius = sqrt(2.0*effectiveDiameter* delMax);
        Mdouble relVelocity = sqrt(4.0*0.4*4./3.*species->getLoadingModulus()*contactRadius*delMax*delMax);
        //v^2=2*2*int(4/3 E sqrt(d) del^3/2)_0^delMax = 4*(2/5*4/3 E sqrt(d) delMax^5/2)
        logger(INFO, " deltaStar: % \n contactRadius: % \n relVelocity: %",
               deltaStar, contactRadius, relVelocity);
        
        particleHandler.clear();
		
		SphericalParticle P0,P1;
        P0.setSpecies(speciesHandler.getObject(0));
        P1.setSpecies(speciesHandler.getObject(0));
		P0.setPosition(Vec3D( 1.0,1.0,1.0));
		P0.setVelocity(Vec3D(0.5*relVelocity,0.0,0.0));
		P0.setRadius(0.5);
		P1.setPosition(Vec3D( 2.0,1.0,1.0));
		P1.setVelocity(Vec3D(-0.5*relVelocity,0.0,0.0));
		P1.setRadius(0.5);
		
		particleHandler.copyAndAddObject(P0);
		particleHandler.copyAndAddObject(P1);
				
		double mass=1.0;
        setTimeStep(species->computeTimeStep(mass)*2.0);
		setTimeMax(getTimeStep()*500.0);
        setFileType(FileType::ONE_FILE);
        restartFile.setFileType(FileType::NO_FILE);
        dataFile.setFileType(FileType::NO_FILE);
        setSaveCount(1);
	}

    void set_chi(double new_){chi=new_;}
    double get_chi() {return chi;}

	double chi;
    HertzianSinterSpecies* species;
};

class LongHertzianSinterForceUnitTest : public HertzianSinterForceUnitTest
{
public:

    LongHertzianSinterForceUnitTest() : HertzianSinterForceUnitTest() {}
    
    void printTime()	const override
    {
        static unsigned int counter = 0;
        if (counter++ > 100)
        {
            DPMBase::printTime();
            counter=0;
        }
    }
    
	void setupInitialConditions() override
    {
        HertzianSinterForceUnitTest::setupInitialConditions();
        setSaveCount(100);
		setTimeMax(1000);
    }
};

int main(int argc UNUSED, char *argv[] UNUSED)
{	
	HertzianSinterForceUnitTest sf;
	double k1=100.0;
	sf.species->setPlasticParameters(k1, 10.0*k1, k1, 1.0);
    sf.species->setDissipation(0);
    sf.setParticleDimensions(3);
    sf.setSystemDimensions(3);

    logger(INFO,"Testing particle particles collision for elastic plastic forces. \n"
        "This will be done for several values of scaled relative velocity chi");
    
    //sf.restartFile.getFstream().precision(20);
    
	//Set up constant data that will be used
    const std::vector<double> chi = {0.34, 0.69, 1.1, 1.37};
    const std::vector<Vec3D> leftFinalVecloity = {
        Vec3D(-0.032721738352012,0.0,0.0),
        Vec3D(-0.0138683231953154,0.0,0.0),
        Vec3D(-0.0204655358555405,0.0,0.0),
        Vec3D(-0.163049415300304,0.0,0.0)};
    const std::vector<Vec3D> leftFinalPosition = {
        Vec3D(0.995546292935715,1.0,1.0),
        Vec3D(1.00695193269955,1.0,1.0),
        Vec3D(1.00840467123501,1.0,1.0),
        Vec3D(0.969386085767181,1.0,1.0)};

    //Loop over all test cases
	for (int i=0; i<4; i++)
    {
        logger(INFO, "Running for chi=%", chi[i]);
		sf.set_chi(chi[i]);
		std::stringstream ss("");
		ss << "HertzianSinterForceUnitTest" << sf.get_chi();
		sf.setName(ss.str().c_str());
		sf.solve();
        sf.writeRestartFile();

//        //Now check the particles are in the right place for each of the 4 cases
//        auto pIt = sf.particleHandler.begin();
//        if (!(*pIt)->getPosition().isEqualTo(leftFinalPosition[i], 1e-10))
//                logger(FATAL,"Left particle is in the wrong position. It is at % and should be %",(*pIt)->getPosition(),leftFinalPosition[i]);
//        if (!(*pIt)->getVelocity().isEqualTo(leftFinalVecloity[i]  , 1e-10))
//                logger(FATAL,"Left particle has the wrong velocity. It is at % and should be %",(*pIt)->getVelocity(),leftFinalVecloity[i]);
    }

    //A longer simulation, with sintering activated
    LongHertzianSinterForceUnitTest lsf;
    lsf.species->setPlasticParameters(k1, 10.0*k1, k1, 1.0);
    lsf.species->setSinterRate(1e-3);
    lsf.species->setDissipation(2e-1);
    lsf.setParticleDimensions(3);
    lsf.setSystemDimensions(3);
	logger(INFO, "Running longer");
    lsf.set_chi(0.47);
    lsf.setName("LongHertzianSinterForceUnitTest");
    lsf.solve();
    lsf.writeRestartFile();
    
    
    
    //(1-((k1+kc)/k2)^1.5) / (1-(k1/k2)^1.5)   
    
    logger(INFO,"Execute 'gnuplot HertzianSinterForceUnitTest.gnu' to view "
                       "output");
    helpers::writeToFile("HertzianSinterForceUnitTest.gnu",
                         "set xlabel 'displacement'\n"
                         "set ylabel 'force'\n"
                         "plot 'HertzianSinterForceUnitTest1.1.fstat' u 7:9, 4./3.*100*x**1.5, 4./3.*1000*(x-0.1)**1.5-4./3.*100*x**1.5, -4./3.*100*x**1.5/(x<0.1)\n"
                         );
    
    logger(INFO,"Execute 'gnuplot LongHertzianSinterForceUnitTest.gnu' to view "
                  "output");
    helpers::writeToFile("LongHertzianSinterForceUnitTest.gnu",
                         "set xlabel 'time'\n"
                         "set ylabel 'displacement'\n"
                         //"plot 'LongHertzianSinterForceUnitTest1.1.fstat' u 7:9 w lp\n"
                         "plot 'LongHertzianSinterForceUnitTest.fstat' u 1:7, sqrt(0.001*(x+120))\n"
                         );
}
