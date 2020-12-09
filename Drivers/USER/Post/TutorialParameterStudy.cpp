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

#include "Mercury3D.h"
#include "Math/Helpers.h"
#include "Species/LinearViscoelasticFrictionSpecies.h"

#include <iostream>
#include <vector>

class TutorialParameterStudy : public Mercury3D
{
public:
    void setupInitialConditions()
    {
        //logger(INFO,"Call to setupInitialConditions()");

        std::vector<int> study_num = DPMBase::get2DParametersFromRunNumber(2,2);

        SphericalParticle P0,P1;

        P0.setSpecies(speciesHandler.getObject(0));
        P1.setSpecies(speciesHandler.getObject(0));

        particleHandler.copyAndAddObject(P0);
        particleHandler.copyAndAddObject(P1);

        //Set up the stuff that is the same for all runs
        particleHandler.getObject(0)->setPosition(Vec3D(0.006,0.005,0.0));
        particleHandler.getObject(1)->setPosition(Vec3D(0.004,0.005,0.0));

        particleHandler.getObject(0)->setVelocity(Vec3D(-0.1,0.0,0.0));
        particleHandler.getObject(1)->setVelocity(Vec3D( 0.1,0.0,0.0));

        //We have three different size particles for both left and right particles
        particleHandler.getObject(0)->setRadius(0.0001*study_num[1]);
        particleHandler.getObject(1)->setRadius(0.0001*study_num[2]);
    }

    void actionsBeforeTimeLoop() override
    {
        std::string fileName = this->dataFile.getName();
        //Set up a 2 by 2 study
        std::vector<int> study_num = DPMBase::get2DParametersFromRunNumber(2,2);

        if (study_num[0] > 0)
        {
            std::cout << "Whole study is complete 1" << std::endl;
            std::exit(0);
        }
        else if (helpers::fileExists(fileName))
        {
            //If it has move on to the next run immedently
            logger(INFO,"This run has been done: %", this->dataFile.getName());
            launchNewRun("./TutorialParameterStudy");
            std::exit(0);
        }
    }

    void actionsAfterSolve() override
    {
        std::vector<int> study_num = DPMBase::get2DParametersFromRunNumber(2,2);

        //If study 0 is complete quit
        if (study_num[0] > 0)
        {
            std::cout << "Whole study is complete 2" << std::endl;
            std::exit(0);
            logger(INFO,"Hanging here");

        }
        else
            //If the study is not complete save the data to disk and move on
        {
            writeRestartFile();
            launchNewRun("./TutorialParameterStudy");
            std::exit(0);
        }

    }


};

int main(int argc UNUSED, char *argv[] UNUSED)
{


    ///Start off by solving the default problem
    TutorialParameterStudy problem;
    problem.autoNumber();

    //make the species of the particle and wall
    LinearViscoelasticFrictionSpecies species;
    species.setDensity(2000);
    //specify contact properties
    //normal forces
    species.setStiffness(1e5);
    species.setDissipation(0.63);
    //tangential (sliding) forces
    species.setSlidingFrictionCoefficient(0.5);
    species.setSlidingStiffness(1.2e4);
    species.setSlidingDissipation(0.16);
    //tangential (rolling) torques
    species.setRollingFrictionCoefficient(0.2);
    species.setRollingStiffness(1.2e4);
    species.setRollingDissipation(6.3e-2);
    //normal (torsion/spin) torques
    species.setTorsionFrictionCoefficient(0.1);
    species.setTorsionStiffness(1.2e4);
    species.setSlidingDissipation(6.3e-2);
    problem.speciesHandler.copyAndAddObject(species);

    ///Autonumber turns on file numbering
    //problem.autoNumber();
    problem.setName("parm_demo"); //elastic_collision
    problem.setTimeStep(1e-5);
    problem.setTimeMax(0.1);
    problem.setSaveCount(500);

    problem.setXMin(0.00);
    problem.setXMax(0.01);
    problem.setYMin(0.00);
    problem.setYMax(0.01);
    problem.setZMin(0.00);
    problem.setZMax(0.01);

    problem.solve();
    return 0;
    exit(0);
}
