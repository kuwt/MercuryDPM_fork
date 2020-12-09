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
//! [REV_ISO:headers]
#include <Mercury3D.h>
#include <Species/LinearViscoelasticSpecies.h>
#include <Boundaries/StressStrainControlBoundary.h>
#include <Boundaries/CubeInsertionBoundary.h>
//! [REV_ISO:headers]

/**
 * This code simulates the isotropic compression in a 3D cuboid box using the StressStrainControlBoundary
 * Particles are inserted into the whole domain (0,0,0) to (1,1,1)
 * and compressed from x-y-z directions with user input constant strain rate tensor.
 */
//! [REV_ISO:class]
class StressStrainControl : public Mercury3D
{
    //! [REV_ISO:construct]
public:
    //Create the class and passing the user inputs using the constructors
    StressStrainControl(const Matrix3D& stressGoal, const Matrix3D& strainRate, const Matrix3D& gainFactor,
                        bool isStrainRateControlled)
            : stressGoal_(stressGoal), strainRate_(strainRate), gainFactor_(gainFactor),
              isStrainRateControlled_(isStrainRateControlled)
    {
        //Define the domain size
        setMin(Vec3D(0, 0, 0));
        setMax(Vec3D(1, 1, 1));
        //Calculate the mass from smallest particle
        double mass = 2000 * constants::pi * mathsFunc::cubic(2.0 * 0.08) / 6.0;
        //Define the species for the particles
        LinearViscoelasticSpecies species;
        species.setDensity(2000);
        species.setStiffnessAndRestitutionCoefficient(10000, 0.4, mass);
        speciesHandler.copyAndAddObject(species);
    }
    //! [REV_ISO:construct]

    //! [REV_ISO:setIni]
private:
    //Here we set up the system parameters, adding particles and define boundaries
    void setupInitialConditions() override
    {
        //Set up the micro parameters.
        double rhop = 2000;     //particle density
        double K1 = 10000;      //particle stiffness
        double en = 0.4;        //restitution coefficient
        double radius = 0.08;   //particle radius
        //Calculate the mass from smallest particle
        double mass = rhop * constants::pi * mathsFunc::cubic(2.0 * radius) / 6.0;
        //Calculate the contact duration
        double tc = std::sqrt(mass / 2.0 / K1 * (mathsFunc::square(constants::pi) +
                                                 mathsFunc::square(log(en))));
        setSystemDimensions(3); //set the 3d simulation
        setTimeStep(tc / 50); //set the timestep based on the shortest contact duration
        setGravity(Vec3D(0.0, 0.0, 0.0)); //set gravity in the system

        //Add particles
        CubeInsertionBoundary insertion;
        SphericalParticle particle(speciesHandler.getObject(0));
        //Insert 96 particles in the subvolume = 0.8*0.8*0.8 = 0.512, if failed by 100 times, it stops
        insertion.set(&particle, 100, 0.8*getMin(), 0.8*getMax(), Vec3D(0, 0, 0), Vec3D(0, 0, 0), radius, radius);
        insertion.setInitialVolume(1.0);
        insertion.checkBoundaryBeforeTimeStep(this);
        logger(INFO,"Inserted % particles",particleHandler.getSize());

        //Delete all existing boundaries
        boundaryHandler.clear();
        //Set up the deformation mode using the boundaries that are based on user's input
        StressStrainControlBoundary boundary;
        boundary.setHandler(&boundaryHandler);
        boundary.set(stressGoal_, strainRate_, gainFactor_, isStrainRateControlled_);
        boundaryHandler.copyAndAddObject(boundary);

        
    }
    //! [REV_ISO:setIni]

    //! [REV_ISO:datamembers]
    //initialize the private variables for passing the user's inputs in the constructor
    Matrix3D stressGoal_;
    Matrix3D strainRate_;
    Matrix3D gainFactor_;
    bool isStrainRateControlled_;
    //! [REV_ISO:datamembers]
};
//! [REV_ISO:class]

//! [REV_ISO:main]
//Here we define all the control parameters and solve the problem
int main(int argc UNUSED, char* argv[] UNUSED)
{
    //We want strainrate control, so here we set the target Stress to free, all to zero
    Matrix3D stressGoal;
    stressGoal.XX = 0.0;
    stressGoal.YY = 0.0;
    stressGoal.ZZ = 0.0;
    stressGoal.XY = 0.0;

    //Here we set the isotropic strainrate tensor, negative sign means compression
    Matrix3D strainRate;
    strainRate.XX = -0.02;
    strainRate.YY = -0.02;
    strainRate.ZZ = -0.02;
    strainRate.XY = 0.0;

    //This is default gainFactor, if you choose to use stress control, you might want to tune this one
    Matrix3D gainFactor;
    gainFactor.XY = 0.0001;
    gainFactor.XX = 0.0001;
    gainFactor.YY = 0.0001;
    gainFactor.ZZ = 0.0001;

    //This is by default set to true, so particles are controlled by affine movement.
    //If set it to false, then the particles will only follow the boundary movement.
    bool isStrainRateControlled = true;
    
    //Define the object and problem to solve
    StressStrainControl dpm(stressGoal, strainRate, gainFactor, isStrainRateControlled);
    
    //Set name
    dpm.setName("Tutorial_IsotropicCompression");
    //Set output and time stepping properties
    dpm.setTimeMax(10);
    //Set the SaveCount, i.e. 100 timesteps saving one snapshot
    dpm.setSaveCount(100);
    //Currently all the normal file outputs are switched off, simply switch it on by replacing NO_FILE to ONE_FILE
//    dpm.dataFile.setFileType(FileType::NO_FILE);
//    dpm.restartFile.setFileType(FileType::NO_FILE);
//    dpm.fStatFile.setFileType(FileType::NO_FILE);
//    dpm.eneFile.setFileType(FileType::NO_FILE);
    //Set output the particle information in VTK for ParaView Visualisation
    dpm.setParticlesWriteVTK(true);
    //Because of periodic boundary, out put wall files is not necessary in this case
    //dpm.setWallsWriteVTK(true);
    //Solve the problem
    dpm.solve();
}
//! [REV_ISO:main]
