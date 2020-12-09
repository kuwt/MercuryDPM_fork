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
//#include <Species/Species.h>
//#include <Species/LinearViscoelasticSpecies.h>
#include <Species/LinearPlasticViscoelasticFrictionSpecies.h>
#include <Mercury3D.h>
#include <Boundaries/PeriodicBoundary.h>
#include <Boundaries/StressStrainControlBoundary.h>
#include "Boundaries/LeesEdwardsBoundary.h"
#include <CMakeDefinitions.h>

///In this file, we create the selftest for the StressStrainControlBoundary
///Note that here we load the restart file and only tested simple shear xy with strainrate control
///There are multiple ways of defining the control modes, please check documentation for more details.

class StressStrainControl : public Mercury3D
{
public:
    
    StressStrainControl(const Matrix3D& stressGoal, const Matrix3D& strainRate, const Matrix3D& gainFactor,
                        bool isStrainRateControlled)
            : stressGoal_(stressGoal), strainRate_(strainRate), gainFactor_(gainFactor),
              isStrainRateControlled_(isStrainRateControlled)
    {
        setName("StressStrainControlBoundarySelfTestInput");
        readRestartFile(getMercurySourceDir()+"/Drivers/SelfTests/Boundaries/InputData/StressStrainControlBoundarySelfTestInput");
        setRestarted(false);
        particleSpecies = dynamic_cast<LinearPlasticViscoelasticFrictionSpecies*>(speciesHandler.getLastObject());
        for (auto& p : particleHandler)
        {
            p->setSpecies(particleSpecies);
        }
    }


private:
    
    
    void setupInitialConditions() override
    {
        
        //  set particle properties for the species--------------------------------------------------
        const Mdouble mean_particleDiameter = 2.0;        //set particle diameter
        const Mdouble rhop = 2000.0;                //set particle density
        const Mdouble en = 0.804;                    //set restitution coefficient
        const Mdouble K1 = 100000;                //set loading stiffness
        const Mdouble K2 = 100000;                //set unloading stiffness
        const Mdouble Kc = 0.0;                    //set cohesive stiffness
        
        const Mdouble mu_slid = 0.0;            //set sliding friction coefficient
        const Mdouble mu_roll = 0.0;                //set rolling friction coefficient
        const Mdouble mu_tor = 0.0;                //set torsional friction coefficient
        const Mdouble Phic = 0.5;                // penetration DepthMax, the maximum depth of linear plastic-viscoelastic normal force
        const Mdouble poly = 3;                  //set polydispersity
        
        
        
        //calculate the contact duration and safe timestep
        Mdouble rmin = particleHandler.getObject(0)->getRadius();
        Mdouble rmax = particleHandler.getObject(0)->getRadius();
        
        for (auto& p : particleHandler)
        {
            rmin = std::min(rmin, p->getRadius());
            rmax = std::max(rmax, p->getRadius());
            p->setVelocity(Vec3D(0.0, 0.0, 0.0));
            p->setSpecies(particleSpecies);
        }
        
        double mass =
                rhop * constants::pi * mathsFunc::cubic(2.0 * rmin) / 6.0; //calculate the mass from smallest particle
        double tc = std::sqrt(mass / 2.0 / K1 * (mathsFunc::square(constants::pi) +
                                                 mathsFunc::square(log(en)))); //calculate the contact duration
        setTimeStep(tc / 50); //set the timestep based on the shortest contact duration
        setGravity(Vec3D(0.0, 0.0, 0.0)); //set gravity in the system
        
        
        //! particleSpecies    set the species parameters
        particleSpecies->setDensity(rhop);
        particleSpecies->setCollisionTimeAndRestitutionCoefficient(tc, en, mass);
        particleSpecies->setPlasticParameters(K1, K2, Kc, Phic);
        if (mu_slid == 0)
        {
            particleSpecies->setSlidingStiffness(0.0);
        }
        else
        {
            particleSpecies->setSlidingStiffness(2.0 / 10.0 * particleSpecies->getLoadingStiffness());
        }
        if (mu_roll == 0)
        {
            particleSpecies->setRollingStiffness(0.0);
        }
        else
        {
            particleSpecies->setRollingStiffness(2.0 / 10.0 * particleSpecies->getLoadingStiffness());
        }
        if (mu_tor == 0)
        {
            particleSpecies->setTorsionStiffness(0.0);
        }
        else
        {
            particleSpecies->setTorsionStiffness(2.0 / 10.0 * particleSpecies->getLoadingStiffness());
        }
        particleSpecies->setSlidingFrictionCoefficient(mu_slid);
        particleSpecies->setSlidingFrictionCoefficientStatic(mu_slid);
        particleSpecies->setRollingFrictionCoefficient(mu_roll);
        particleSpecies->setRollingFrictionCoefficientStatic(mu_roll);
        particleSpecies->setTorsionFrictionCoefficient(mu_tor);
        particleSpecies->setTorsionFrictionCoefficientStatic(mu_tor);
        particleSpecies->setSlidingDissipation(2.0 / 10.0 * particleSpecies->getDissipation());
        particleSpecies->setRollingDissipation(2.0 / 10.0 * particleSpecies->getDissipation());
        particleSpecies->setTorsionDissipation(2.0 / 10.0 * particleSpecies->getDissipation());
        
        //Set up the new boundaries based on user input
        boundaryHandler.clear();                                //! Delete all exist boundaries
        StressStrainControlBoundary boundary;
        boundary.setHandler(&boundaryHandler);
        boundary.set(stressGoal_, strainRate_, gainFactor_, isStrainRateControlled_);
        boundaryHandler.copyAndAddObject(boundary);

        
    }
    
    LinearPlasticViscoelasticFrictionSpecies* particleSpecies;
    Matrix3D stressGoal_;
    Matrix3D strainRate_;
    Matrix3D gainFactor_;
    bool isStrainRateControlled_;
    
};

int main(int argc UNUSED, char* argv[] UNUSED)
{
    //helpers::writeToFile("xyz.restart","abc");
    
    Matrix3D stressGoal;
    stressGoal.XX = 0.0;
    stressGoal.YY = 0.0;
    stressGoal.ZZ = 0.0;
    stressGoal.XY = -2.0;
    
    Matrix3D strainRate;
    strainRate.XY = 0.0;
    Matrix3D gainFactor;
    gainFactor.XY = 0.0001;
    gainFactor.XX = 0.0001;
    gainFactor.YY = 0.0001;
    gainFactor.ZZ = 0.0001;
    bool isStrainRateControlled = true;
    
    //instantiate the class
    StressStrainControl dpm(stressGoal, strainRate, gainFactor, isStrainRateControlled);
    
    //set name
    dpm.setName("StressStrainControlBoundarySelfTest");
    //set output and time stepping properties
    dpm.setTimeMax(0.1);
    dpm.eneFile.setSaveCount(1000);
    //dpm.dataFile.setFileType(FileType::NO_FILE);
    //dpm.restartFile.setFileType(FileType::ONE_FILE);
    //dpm.fStatFile.setFileType(FileType::NO_FILE);
    dpm.eneFile.setFileType(FileType::ONE_FILE);
    //solve
    dpm.solve();
    //dpm.writeRestartFile();
}
