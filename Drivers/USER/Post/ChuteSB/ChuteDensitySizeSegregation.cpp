//Copyright (c) 2013-2014, The MercuryDPM Developers Team. All rights reserved.
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

#include <sstream>
#include <iostream>
#include <iomanip>
#include <cmath>

//This code is based of chute
#include "Chute.h"
#include "Boundaries/PeriodicBoundary.h"
#include "Walls/InfiniteWall.h"

#include <sys/types.h>
#include <sys/stat.h>
#include <Species/LinearViscoelasticFrictionSpecies.h>

using namespace std;

/**
 * \brief This class does segregation problems in a periodic chute
 * It uses species to create two type of partices.
 * It the sets contact properties of the collisions such that coefficient of resitution and contact time are the same for all collisions
 **/
class SegregationPeriodic : public Chute{
public:

///Allows the user to set what is written into the info column in the data file
double getInfo(const BaseParticle &P) const {
	return P.getIndSpecies();
	}


///This code requires you do not nothing special after each time step
void actionsBeforeTimeStep(){};

///This is the info call
void write(std::ostream& os, bool print_all = false) const
{
//	os << "This is a density size-segregation chute code problem" << endl;
//	os << "\n \n \n"<< endl;
	 MercuryBase::write(os, print_all);
	 os << "particle specie-0 size     : " << radius_0 << endl;
	 os << "particle specie-1 size     : " << radius_1 << endl;
	 os << "particle specie-2 size     : " << radius_2 << endl;
	 os << "particle specie-0 rho      : " << rho_0 << endl;
	 os << "particle specie-1 rho      : " << rho_1 << endl;
	 os << "particle specie-2 rho      : " << rho_2 << endl;
	 os << "particle size-ratio        : " << radius_2/radius_1 << endl;
	 os << "particle density-ratio     : " << rho_2/rho_1 << endl;
         os << "No. of species-1 particles : " << N_1 << endl;
         os << "No. of species-2 particles : " << N_2 << endl;
}

/// This setup the intial conditions, generates volume fraction of particle 1.
/// Sets the program to be periodic in x.
/// \bug This code is not non-dimensionalised at the moment, should do this shortly, but at the moment. Should swap this to Silbert particles shortly
void setupInitialConditions()
{

	//Check if the run has been done before. If yes, skip and start next run
    //if (helpers::fileExists(getDataFile().getName()))
	//{
	//	//If it has move on to teh next run immedently
	//	cout << "This run has been done " << endl;
	//	//launch_new("density-size-segregation",true);
	//	exit(0);
	//}

	//Set up a 12 by 12 study

	//vector<int> study_num=get2DParametersFromRunNumber(3,1);
    //std::cout<<"study_num = "<<study_num[0]<<study_num[1]<<study_num[2]<<std::endl;

	//If study 0 is complete quit
	//if (study_num[0] > 0)
	//{
	//	cout << "Study is complete " << endl;
	//	exit(0);
	//}
	//else
	//If the study is not complete save the data to disk and move on
	//{

    //	writeRestartFile();
    //      launchNewRun("segregation");
	//}

	//  PARTICLE PROPERTIES//
	/////////////////////////

	//********************ratio of radii / diameters (not volumes!)****************************************
	double shat = sizeRatio;//(0.2 + 0.1*study_num[1]); // d2/d1 // 1.7; //
	double rhat = densityRatio;//study_num[2];//(0.2 + 0.1*study_num[2]) // 1.0; // // rho2/rho1


	//
	// All the parameters are non-dimensional
	// mean particle dia = 1 , mean particle mass = 1 and gravity g=1
	// implying the mean particle density = 6/pi

	// Box dimensions, Volume of the box = L*W*H
	double Vbox = 15.0*10.0*40.0;

	//******************** volume fraction of species type-1*********************************************
	double phi = 1.0;

	// total volume occupied by the particles
	//double Vp = Vbox*(constants::pi/6.);

	// particle diameters
	double dm = 1.0; // mean particle
	double d1 = 1./(phi + (1-phi)*shat); // species type-1
	double d2 = d1*shat; // species type-2
	double tolo = 1.e-12;
        if (abs(d1 - d2) <= tolo)
        {
                d2 = d1 - 0.0001;
        }

	// particle species radii
	radius_0 = 0.5*dm; // mean particle
	radius_1 = 0.5*d1;
	radius_2 = 0.5*d2;

	// particle densities
	rho_0 = 6./constants::pi;
	rho_1 = (6./constants::pi)/(phi + (1-phi)*rhat);
	rho_2 = rhat*rho_1;

	// no of particles
	int N1 = (phi*Vbox)/pow(d1,3);
	int N2 = 1.+((1-phi)*Vbox)/pow(d2,3);
	N_1 = N1;
    N_2 = N2;
	//
	double tc  = 2.5e-4;//0.005*pow(d,0.5);//1.e-2;//2.e-1;//1e-2;//5e-2;//5e-3;
	double r   = 0.88249690258;

	//
        double mass_0 = 1.0;// mean particle mass
	double mass_1 = 4./3.*constants::pi*pow(radius_1,3.0)*rho_1;
	double mass_2 = 4./3.*constants::pi*pow(radius_2,3.0)*rho_2;
	//

	double mass_00 = (mass_0*mass_0)/(mass_0 + mass_0);
	double mass_10 = (mass_1*mass_0)/(mass_1 + mass_0);
        double mass_20 = (mass_2*mass_0)/(mass_2 + mass_0);
	double mass_11 = (mass_1*mass_1)/(mass_1 + mass_1);
	double mass_22 = (mass_2*mass_2)/(mass_2 + mass_2);
	double mass_12 = (mass_1*mass_2)/(mass_1 + mass_2);
	//

	double gamma_n00 = -2.0*mass_00*log(r)/tc;
	double gamma_n11 = -2.0*mass_11*log(r)/tc;
	double gamma_n22 = -2.0*mass_22*log(r)/tc;
        double gamma_n12 = -2.0*mass_12*log(r)/tc;
        double gamma_n10 = -2.0*mass_10*log(r)/tc;
        double gamma_n20 = -2.0*mass_20*log(r)/tc;
	//
        double const1 = pow(tc/constants::pi,2.0);
        double k_n00 = mass_00*(1./const1 + pow(gamma_n00/(2*mass_00),2.0));
        double k_n11 = mass_11*(1./const1 + pow(gamma_n11/(2*mass_11),2.0));
        double k_n22 = mass_22*(1./const1 + pow(gamma_n22/(2*mass_22),2.0));
        double k_n12 = mass_12*(1./const1 + pow(gamma_n12/(2*mass_12),2.0));
        double k_n10 = mass_10*(1./const1 + pow(gamma_n10/(2*mass_10),2.0));
        double k_n20 = mass_20*(1./const1 + pow(gamma_n20/(2*mass_20),2.0));
	//

        auto S0 = speciesHandler.copyAndAddObject(LinearViscoelasticFrictionSpecies());
        auto S1 = speciesHandler.copyAndAddObject(S0);
        auto S2 = speciesHandler.copyAndAddObject(S1);
        auto S01 = speciesHandler.getMixedObject(S0, S1);
        auto S02 = speciesHandler.getMixedObject(S0, S2);
        auto S12 = speciesHandler.getMixedObject(S1, S2);
        //
        S0->setDensity(rho_0);
        S0->setStiffness(k_n00);
        S0->setSlidingStiffness((2.0/7.0)*k_n00);
        S0->setDissipation(gamma_n00);
        S0->setSlidingDissipation(S0->getDissipation());//  Set the tangential dissipation equal to the normal disipation for 1-1 collsions
	setInflowParticleRadius(std::min(radius_1,radius_2),std::max(radius_1,radius_2));
	//set_HGRID_cell_to_cell_ratio(1.00001*max(radius_1,radius_2)/min(radius_1,radius_2));
	//set_HGRID_num_buckets_to_power(particleHandler.getNumberOfObjects());
        S0->setSlidingFrictionCoefficient(0.5);
	//
	S1->setDensity(rho_1);
        S1->setStiffness(k_n11);
	double k_t=2.0/7.0*k_n11;
	S1->setSlidingStiffness(k_t);
        S1->setDissipation(gamma_n11);
	S1->setSlidingFrictionCoefficient(0.5);
	S1->setSlidingDissipation(S1->getDissipation());//  Set the tangential dissipation equal to the normal disipation for 2-2 c
	//
	S2->setDensity(rho_2);
        S2->setStiffness(k_n22);
	k_t=2.0/7.0*k_n22;
	S2->setSlidingStiffness(k_t);
        S2->setDissipation(gamma_n22);
	S2->setSlidingFrictionCoefficient(0.5);
	S2->setSlidingDissipation(S2->getDissipation());//  Set the tangential dissipation equal to the normal disipation for 2-2 collision
	//
        S01->setStiffness(k_n10);
	k_t=2.0/7.0*k_n10;
	S01->setSlidingStiffness(k_t);
        S01->setDissipation(gamma_n10);
	S01->setSlidingDissipation(S01->getDissipation());
	S01->setSlidingFrictionCoefficient(0.5);
	//
        S12->setStiffness(k_n12);
	k_t=2.0/7.0*k_n12;
	S12->setSlidingStiffness(k_t);
        S12->setDissipation(gamma_n12);
	S12->setSlidingDissipation(S12->getDissipation());
	S12->setSlidingFrictionCoefficient(0.5);
        //
        S02->setStiffness(k_n20);
	k_t=2.0/7.0*k_n20;
	S02->setSlidingStiffness(k_t);
        S02->setDissipation(gamma_n20);
	S02->setSlidingDissipation(S02->getDissipation());
	S02->setSlidingFrictionCoefficient(0.5);
        //
	//Setup the base i.e. the chute particles - This has to be done after the particle properties are set, but the inflow partilces are created.
	Chute::createBottom();

	InfiniteWall w0;
	w0.set(Vec3D(0.0,0.0,-1.0), -(getZMin()-0.5));
        w0.setSpecies(S0);
        wallHandler.copyAndAddObject(w0);

        //w0.set(Vec3D(0.0,0.0,1.0), getZMax());
        //w0.setSpecies(S0);
        //wallHandler.copyAndAddObject(w0);

	PeriodicBoundary b0;
        b0.set(Vec3D(1.0,0.0,0.0), getXMin(), getXMax());
        boundaryHandler.copyAndAddObject(b0);

        b0.set(Vec3D(0.0,1.0,0.0), getYMin(), getYMax());
        boundaryHandler.copyAndAddObject(b0);

        // Create Particles
        std::cout << "no. of N1 particles" << N1 << std::endl;
        std::cout << "no. of N2 particles" << N2 << std::endl;
        //
        inflowParticle_.setSpecies(S2);
        Mdouble radius; // set the radius for the first particle to be inserted
        radius = radius_2;
        inflowParticle_.setRadius(radius);
        //
        int numberOfN1ParticlesInserted=0;
        int numberOfN2ParticlesInserted=0;
        while( (N1>0) || (N2>0) )
        {

          //specify random components separately
          Vec3D pos;
          pos.X = random.getRandomNumber(getXMin()+radius, getXMax()-radius);
          pos.Y = random.getRandomNumber(getYMin()+radius, getYMax()-radius);
          pos.Z = random.getRandomNumber(getZMin()+radius, getZMax()-radius);
          inflowParticle_.setPosition(pos);

          if (checkParticleForInteraction(inflowParticle_))
          {
              particleHandler.copyAndAddObject(inflowParticle_);
              if (radius == radius_1)
              {
                  N1--;
                  ++numberOfN1ParticlesInserted;
                  //std::cout << "N1 particles inserted" << numberOfN1ParticlesInserted << std::endl;
              }
              else
              {
                  N2--;
                  ++numberOfN2ParticlesInserted;
                  //std::cout << "N2 particles inserted" << numberOfN2ParticlesInserted << std::endl;
              }

              //only modify the radius if a particle is successfully inserted
              int grn = random.getRandomNumber(1,N1+N2);
              //if( random.getRandomNumber(1,N1+N2) > N2)
              //std::cout << "grn=" << grn << std::endl;
              if( grn > N2)
              {
                  radius = radius_1;
                  inflowParticle_.setSpecies(S1);
                  inflowParticle_.setRadius(radius);
              }
              else
              {
                  radius = radius_2;
                  inflowParticle_.setSpecies(S2);
                  inflowParticle_.setRadius(radius);
              }
          }
          else
              setZMax(getZMax()+.0001);
        }
//
	// CREATE THE PARTICLES
/*
        int numberOfN1ParticlesInserted=0;
        int numberOfN2ParticlesInserted=0;
        while( (N1>0) || (N2>0) )
        {
        if( random.getRandomNumber(1,N1+N2) > N1)
        {
                inflowParticle_.setRadius(radius_2);
                inflowParticle_.setSpecies(S2);
                inflowParticle_.setPosition(Vec3D(random.getRandomNumber(getXMin()+inflowParticle_.getRadius(),getXMax()-inflowParticle_.getRadius()),
                                      random.getRandomNumber(getYMin()+inflowParticle_.getRadius(),getYMax()-inflowParticle_.getRadius()),
                                      random.getRandomNumber(getZMin(),getZMax())));
                inflowParticle_.computeMass();
                if (checkParticleForInteraction(inflowParticle_))
                {
                    particleHandler.copyAndAddObject(inflowParticle_);
                    N2--;
                    ++numberOfN2ParticlesInserted;
                    std::cout << "Created N2 particle " << numberOfN2ParticlesInserted << std::endl;
                }
                else
                {
                    setZMax(getZMax()+.001);
                }
        }

        else
        {
            inflowParticle_.setRadius(radius_1);
            inflowParticle_.setSpecies(S1);
            inflowParticle_.setPosition(Vec3D(random.getRandomNumber(getXMin()+inflowParticle_.getRadius(),getXMax()-inflowParticle_.getRadius()),
                                      random.getRandomNumber(getYMin()+inflowParticle_.getRadius(),getYMax()-inflowParticle_.getRadius()),
                                      random.getRandomNumber(getZMin(),getZMax())));
            inflowParticle_.computeMass();
            if (checkParticleForInteraction(inflowParticle_))
            {
                particleHandler.copyAndAddObject(inflowParticle_);
                N1--;
                ++numberOfN1ParticlesInserted;
                std::cout << "Created a N1 particle " << numberOfN1ParticlesInserted << std::endl;
            }
            else
            {
               setZMax(getZMax()+.001);
            }
        }
        }
*/
	//Write the info to the screen and save a copy to the disk

	std::cout << "Finished creating N1 particles=" << numberOfN1ParticlesInserted << "Finished creating N2 particles=" << numberOfN2ParticlesInserted << std::endl;

        write(std::cout,false);
        writeRestartFile();

}

void actionsAfterTimeStep()
{
    if(getTime() >= 10.0)
    {
        setChuteAngleAndMagnitudeOfGravity(23.0,1.0);
        setSaveCount(2000);
    }
}

void setSizeAndDensityRatio (double sr, double dr)
{
    sizeRatio = sr;
    densityRatio = dr;
}

private:
double rho_0;
double rho_1;
double rho_2;
double radius_0;
double radius_1;
double radius_2;
int N_1;
int N_2;
SphericalParticle inflowParticle_;
double sizeRatio;
double densityRatio;

};


int main(int argc UNUSED, char *argv[] UNUSED)
{
	SegregationPeriodic problem;
	//std::vector<int> studyNum=problem.get2DParametersFromRunNumber(3,1);
    //std::cout<<"studyNum = "<<studyNum[0]<<studyNum[1]<<studyNum[2]<<std::endl;

	//std::stringstream nameStream;
	//std::string nameBase = "Chute";
    //nameStream << nameBase <<"_study_"<<studyNum[0]<<"_"<<studyNum[1]<<"_"<<studyNum[2];
    //std::string fileName = nameStream.str();
    //problem.setName(fileName);
    //problem.incrementRunNumberInFile();

	double dRatio;
	double rhoRatio;

	if (argc < 2)
    {
        double aCounter = problem.readRunNumberFromFile();
        std::cout<<"aCounter = "<<aCounter<<std::endl;

        dRatio = 1.0+0.25*aCounter;
        rhoRatio = 1.0;

        std::stringstream nameStream;
        std::string nameBase = "Chute_";
        nameStream << nameBase << "Sr" << dRatio << "Dr" << rhoRatio << "Box_15x10x40";
        std::string fullName = nameStream.str();
        problem.setName(fullName);
        problem.incrementRunNumberInFile();

    }
    if (argc>1)
    {
        //fileName = argv[1];
        //std::cout << "restarting from " << fileName << std::endl;
        //problem.readRestartFile(fileName);
        //problem.setAppend(true);
    }

	//This should be set to 100 for full problems.
	problem.setTimeMax(3000);
        //problem.speciesHandler.copyAndAddObject(S0);
        //problem.speciesHandler.copyAndAddObject(S1);
        //problem.speciesHandler.copyAndAddObject(S2);

        //problem.autoNumber();
        //std::cout << "studyNum2=" << studyNum[2] << std::endl;
        Mdouble gravity= 1.0;//2.0 + (studyNum[2]-1)*2.0;
	//std::cout << "Gravity" << gravity << std::endl;
        // Chute properties
	problem.setFixedParticleRadius(0.5);
	problem.setRoughBottomType(MULTILAYER);
    problem.setChuteAngleAndMagnitudeOfGravity(0.0,gravity);// Is changed to 26 degrees after t=10 DPM units
	problem.setChuteLength(15);
	problem.setChuteWidth(10);
	problem.setZMax(30);
	problem.setMaxFailed(6);
	problem.makeChutePeriodic();
	problem.setSizeAndDensityRatio(dRatio,rhoRatio);//size ratio and density ratio

	//solve
	problem.setTimeStep(2.5e-4 / 50.0);//initially it was auto-set
	problem.setSaveCount(4000);
	std::cout << "dt=" << problem.getTimeStep() << std::endl;

	//problem.random.randomise();
        //
        problem.dataFile.setFileType(FileType::MULTIPLE_FILES);
        problem.restartFile.setFileType(FileType::ONE_FILE);
        problem.fStatFile.setFileType(FileType::NO_FILE);
        problem.eneFile.setFileType(FileType::NO_FILE);

	//This set to colouring based of size and small vectors
	problem.setXBallsColourMode(7);
	problem.setXBallsVectorScale(1);
	problem.setXBallsAdditionalArguments("-v0 -solidf");

	//solves the problems
	problem.solve(argc,argv);
}
