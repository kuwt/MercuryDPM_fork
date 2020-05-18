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

#include <sstream>
#include <iostream>
#include <iomanip>
#include <cmath>

//This code is based on the chute
#include "Chute.h"
#include "Boundaries/PeriodicBoundary.h"

#include <sys/types.h>
#include <sys/stat.h>
#include <Species/LinearViscoelasticFrictionSpecies.h>

using namespace std;

/** \brief This class does segregation problems in a periodic chute
 * It uses species to create two type of particles. One for the large and one for the small
 * It the sets contact properties of the collisions such that coefficient of restitution and contact time are the same for all collisions
 */
class SegregationPeriodic : public Chute
{
public:

    ///This code requires you do not nothing special after each time step
    void actionsBeforeTimeStep()
    {
    }

    ///This is the info call
    //void write(  	std::ostream &   	 os, bool  	print_all = false)
    //{
    //	os << "This is a segregation chute code problem " << endl;
    //	os << "\n \n \n"<< endl;
    //
    //
    //	 HGRID_base::write(os, print_all);
    //
    //	 os << "Large particle size : " << radius_l << endl;
    //	 os << "Small particle size : " << radius_s << endl;
    //
    //
    //
    //}


    /// This setup the initial conditions, generates small volume fraction of particles.
    /// Sets the program to be periodic in x.
    /// \bug This code is not non-dimensionalised at the moment, should do this shortly, but at the moment. Should swap this to Silbert particles shortly
    void setupInitialConditions()
    {

        //Check if the run has been done before. If yes, skip and start next run
        if (helpers::fileExists(dataFile.getName()))
        {
            //If it has move on to the next run immediately
            cout << "This run has been done " << endl;
            launchNewRun("./segregation", true);
            exit(0);
        }

        //Set up a 10 by 10 study
        vector<int> study_num = get2DParametersFromRunNumber(10, 1);


        /*
         This part was in setupInitialConditions, but then creates an infinite loop of setting up initial conditions.
         It might be better to put it in tasksAfterSolve, or something like that.
        //If study 0 is complete quit
        if (study_num[0] > 0)
        {
            cout << "Study is complete " << endl;
            exit(0);
        }
        else //If the study is not complete save the data to disk and move on
        {
            writeRestartFile();
            launchNewRun("./segregation");
        }*/

        //CREATE THE WALLS//
        ////////////////////
        createWalls();

        //  PARTICLE PROPERTIES//
        /////////////////////////

        //Number of small particles
        int numberOfSmallParticles = 10;
        //Small particle radius
        radius_s = 0.5;
        //Radius of large particles, changes from study to study.
        radius_l = radius_s * (1.0 + study_num[1] / 10.0);
        //Number of large particles, fixed to the keep the volume fraction of large and small particles equal.
        int numberOfLargeParticles = pow(radius_s / radius_l, 3) * numberOfSmallParticles;

        setSpeciesProperties();

        //Setup the base i.e. the chute particles - This has to be done after the particle properties are set, but before the inflow particles are created.
        Chute::setupInitialConditions();

        // CREATE THE PARTICLES
        createParticles(numberOfSmallParticles, numberOfLargeParticles);

        //Write the info to the screen and save a copy to the disk
        cout << "Finished creating particles" << endl;
        write(std::cout, false);
        writeRestartFile();

        setChuteProperties();

    }

    void setSpeciesProperties()
    {
        //Set the contact time (tc), restitution coefficient (r) and density (rho) for small for all particles
        double tc = 1e-5;
        double r = 0.88;
        double rho = 6 / constants::pi;

        double mass_small = 4 / 3 * constants::pi * pow(radius_s, 3.0) * rho;
        double mass_large = 4 / 3 * constants::pi * pow(radius_l, 3.0) * rho;

        auto S0 = speciesHandler.copyAndAddObject(LinearViscoelasticFrictionSpecies());
        auto S1 = speciesHandler.copyAndAddObject(S0);
        auto S01 = speciesHandler.getMixedObject(S0, S1);

        S0->setDensity(rho);
        S0->setCollisionTimeAndRestitutionCoefficient(tc, r, mass_small);
        S0->setSlidingDissipation(S0->getDissipation()); //  Set the tangential dissipation equal to the normal dissipation for small-small collisions
        setInflowParticleRadius(0.5, 1.0);
        S0->setSlidingFrictionCoefficient(0.5);
        S1->setCollisionTimeAndRestitutionCoefficient(tc, r, mass_large);

        S1->setSlidingDissipation(S1->getDissipation()); //  Set the tangential dissipation equal to the normal dissipation for large-large collision
        S01->setCollisionTimeAndRestitutionCoefficient(tc, r, mass_small, mass_large);
        S01->setSlidingDissipation(S01->getDissipation()); //  Set the tangential dissipation equal to the normal dissipation for mixed collision

    }

    void createWalls()
    {
	    PeriodicBoundary B0;
	    B0.set(Vec3D(1, 0, 0), getXMin(), getXMax());
	    boundaryHandler.copyAndAddObject(B0);
    }

    void createParticles(int numberOfSmallParticles, int numberOfLargeParticles)
    {
        //Generate a large particle: set radius to large radius subtract one of the list of large particles to be generated
        inflowParticle_.setRadius(radius_l);
        inflowParticle_.setSpecies(speciesHandler.getObject(1));
        numberOfLargeParticles--;

        //randomize particle position, zero initial velocity
        inflowParticle_.setPosition(Vec3D(random.getRandomNumber(getXMin(), getXMax()),
        random.getRandomNumber(getYMin(), getYMax()),
        random.getRandomNumber(getZMin(), getZMax())));
        inflowParticle_.setVelocity(Vec3D(0.0, 0.0, 0.0));

        //Add the new particle to the list of current particles
        particleHandler.copyAndAddObject(inflowParticle_);
        hGridRebuild();

        while ((numberOfSmallParticles > 0) && (numberOfLargeParticles > 0))
        {
            //random to see if want to generate a large or small particles, helps makes the initial conditions homogeneous
            if (random.getRandomNumber(1.0, numberOfLargeParticles + numberOfSmallParticles) > numberOfLargeParticles)
            {
                //Generate a small particle: set radius to small radius subtract one off the list of small particles to be generated
                inflowParticle_.setRadius(radius_s);
                inflowParticle_.setSpecies(speciesHandler.getObject(0));
                numberOfSmallParticles--;
            }
            else
            {
                //Generate a large particle: set radius to large radius subtract one of the list of large particles to be generated
                inflowParticle_.setRadius(radius_l);
                inflowParticle_.setSpecies(speciesHandler.getObject(1));
                numberOfLargeParticles--;
            }

            //randomize particle position, zero initial velocity
            inflowParticle_.setPosition(Vec3D(random.getRandomNumber(getXMin(), getXMax()),
                                              random.getRandomNumber(getYMin(), getYMax()),
                                              random.getRandomNumber(getZMin(), getZMax())));
            inflowParticle_.setVelocity(Vec3D(0.0, 0.0, 0.0));

            //Add the new particle to the list of current particles
            particleHandler.copyAndAddObject(inflowParticle_);
        }
    }

    void setChuteProperties()
    {
        // Chute properties
        setFixedParticleRadius(0.5);
        setRoughBottomType(MONOLAYER_DISORDERED);
        setChuteAngleAndMagnitudeOfGravity(25.0, 1.0);
        setChuteLength(20.0);
        setChuteWidth(10.0);
        setZMax(10.0);
        setMaxFailed(6);
        makeChutePeriodic();
    }

private:
	double radius_s;
	double radius_l;
	SphericalParticle inflowParticle_;
};


int main(int argc UNUSED, char *argv[] UNUSED)
{
    SegregationPeriodic problem;

    // Problem parameters, name tmax and two types of particles
    problem.setName("segregation");
    //This should be set to 100 for full problems.
    problem.setTimeMax(.1);
    problem.setTimeStep(1e-4);
    problem.setupInitialConditions();

    //solve
    ///\todo TW we need a replacement for BaseParticle::calculateMaximumVelocity
    //std::cout << "Maximum allowed speed of particles: " << problem.particleHandler.getSmallestParticle()->calculateMaximumVelocity() << std::endl; // speed allowed before particles move through each other!
    //problem.setTimeStepByParticle();
    //This is based on the fact in general you get too much data, so prob at worst you want to turn it into a 20 at 60fps (which is its self overkill)
    problem.setSaveCount(helpers::getSaveCountFromNumberOfSavesAndTimeMaxAndTimeStep(20 * 60, problem.getTimeMax(), problem.getTimeStep()));
    //problem.setSaveCount(helpers::getSaveCountFromNumberOfSavesAndTimeMaxAndTimeStep(20*60,getTimeMax(),getTimeStep()));
    //problem.setSaveCount(1);

    problem.autoNumber();


    //This set to colouring based of size and small vectors
    problem.setXBallsColourMode(7);
    problem.setXBallsVectorScale(1);
    problem.setXBallsAdditionalArguments("-v0 -solidf");


    //solves the problems
    problem.solve();

    //Make sure the restart data is upto date at the end
    //problem.writeRestartFile();
}


