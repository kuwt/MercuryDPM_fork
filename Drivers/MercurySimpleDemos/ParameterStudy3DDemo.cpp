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
//! [PAR_SIM3D:headers]
#include "Mercury3D.h"
#include "Math/Helpers.h"
#include "Species/LinearViscoelasticFrictionSpecies.h"
#include <iostream>
#include <vector>
//! [PAR_SIM3D:headers]

//! [PAR_SIM3D:class]
class ParameterStudy3DDemo : public Mercury3D
{
public:
    ParameterStudy3DDemo() = default;
    ~ParameterStudy3DDemo() = default;

    //! [PAR_SIM3D:setupInit]
    void setupInitialConditions() override
    {
        studyNum = getCurrentStudyNum();
        createSpecies();

        SphericalParticle P;
        P.setSpecies(speciesHandler.getObject(0));
        // Create 2 particles
        particleHandler.copyAndAddObject(P);
        particleHandler.copyAndAddObject(P);
        //Set up the stuff that is the same for all runs
        particleHandler.getObject(0)->setPosition(Vec3D(0.006,0.005,0.0));
        particleHandler.getObject(1)->setPosition(Vec3D(0.004,0.005,0.0));
        particleHandler.getObject(0)->setVelocity(Vec3D(-0.1,0.0,0.0));
        particleHandler.getObject(1)->setVelocity(Vec3D( 0.1,0.0,0.0));
        // We have two different size particles for both left and right particles
        // The value gets updated based on what studyNum is, these come from the automatic counter
        particleHandler.getObject(0)->setRadius(0.0001);
        particleHandler.getObject(1)->setRadius(0.0001);
    }
    //! [PAR_SIM3D:setupInit]

    //! [PAR_SIM3D:actBefTLoop]
    void actionsBeforeTimeLoop() override
    {
        std::string fileName = this->dataFile.getName();
        if (studyNum[0] > 0)
        {
            std::cout << "Study is complete." << std::endl;
            std::exit(0);
        }
        else if (helpers::fileExists(fileName))
        {
            //If it has move on to the next run immediately
            logger(INFO,"This run has been done: %", this->dataFile.getName());
            logger(INFO,"Launching next run \n\n");
            launchNewRun("./ParameterStudy3DDemo");
            std::exit(0);
        }
    }
    //! [PAR_SIM3D:actBefTLoop]

    //! [PAR_SIM3D:actAftSolve]
    void actionsAfterSolve() override
    {
        if (studyNum[0] > 0)
        {
            std::cout << "Study was already completed, no new simulations ran." << std::endl;
            std::exit(0);
        }
        else
        {
            //If the study is not complete save the data to disk and move on
            writeRestartFile();
            logger(INFO,"Launching next run \n\n");
            launchNewRun("./ParameterStudy3DDemo");
            std::exit(0);
        }
    }
    //! [PAR_SIM3D:actAftSolve]

    //! [PAR_SIM3D:createSpecies]
    void createSpecies()
    {
        //Lets say we want to update the particle stiffness, dissipation and density based on studyNum.
        double k_ref = 1.0e5;
        double diss_ref = 0.63;
        double rho_ref = 2000.0;

        // Arbitrary change to the values based on the three study numbers stored in studyNum
        double k_run = k_ref * studyNum[1];
        double diss_run = diss_ref + 0.05*studyNum[2];
        double rho_run = rho_ref + (studyNum[3] - 1)*500;

        // Setup the contact model for the species
        LinearViscoelasticFrictionSpecies species;
        species.setDensity(rho_run);
        //specify contact properties
        //normal forces
        species.setStiffness(k_run);
        species.setDissipation(diss_run);
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
        speciesHandler.copyAndAddObject(species);
    }
    //! [PAR_SIM3D:createSpecies]

    //! [PAR_SIM3D:setgetFunc]
    void setStudyDimensions(int dim1, int dim2, int dim3)
    {
        dim1_ = dim1;
        dim2_ = dim2;
        dim3_ = dim3;
    }

    std::vector<int> getCurrentStudyNum()
    {
        return DPMBase::get3DParametersFromRunNumber(dim1_, dim2_, dim3_);
    }
    //! [PAR_SIM3D:setgetFunc]

private:
    //! [PAR_SIM3D:privMem]
    int dim1_ = 0;
    int dim2_ = 0;
    int dim3_ = 0;
    std::vector<int> studyNum;
    //! [PAR_SIM3D:privMem]
};
//! [PAR_SIM3D:class]

//! [PAR_SIM3D:main]
int main(int argc UNUSED, char *argv[] UNUSED)
{
    //Problem constructor, using class definition above
    ParameterStudy3DDemo problem;
    // Setup the study dimensions that you want to perform
    problem.setStudyDimensions(2,2,3);
    // Setup the automatic numbering system (create/read COUNTER_DO_NOT_DEL file)
    problem.autoNumber();
    // Setup output file name
    problem.setName("ParameterStudy3DDemo");
    // General simulation settings
    problem.setTimeStep(1e-5);
    problem.setTimeMax(0.1);
    problem.setSaveCount(5000);
    problem.setXMin(0.0);
    problem.setXMax(1.0);
    problem.setYMin(0.0);
    problem.setYMax(1.0);
    problem.setZMin(0.0);
    problem.setZMax(1.0);

    // Solve the problem
    problem.solve();
    // Give a successful exit code to the terminal
    return 0;
}
//! [PAR_SIM3D:main]
