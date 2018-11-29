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
//#include <Species/Species.h>
//#include <Species/LinearViscoelasticSpecies.h>
#include <Species/LinearPlasticViscoelasticFrictionSpecies.h>
#include <Mercury3D.h>
#include <Boundaries/PeriodicBoundary.h>
#include <Boundaries/StressStrainControlBoundary.h>
#include <Boundaries/CubeInsertionBoundary.h>
#include <Walls/AxisymmetricIntersectionOfWalls.h>
#include <Boundaries/DeletionBoundary.h>
#include <Species/LinearPlasticViscoelasticFrictionReversibleAdhesiveSpecies.h>
#include "Boundaries/LeesEdwardsBoundary.h"
#include "addSpecies.h"

/**
 * This simulation consists of two parts:
 * A)
 * Particles, with given volume fracion, are
 * placed in a cylindral column with a bottom, no lid.
 * Then all particles are relaxed.
 * B)
 * Remove the cylindrical wall and wait until the particles
 * are not move anymore, then we extract the angle of repose.
 * Note:
 * user can set baseRadius, fillHeight, and give
 * input material name like APAPP with possiblity of changing
 * restitutionCoefficient, relativeCohesionStiffness,
 * slidingFriction and rollingFriction.
 **/
class GranuHeap : public Mercury3D
{
public:
    
    /**
     *
     * @param stressGoal the stress tensor that the material should be subjected to
     * @param strainRate the applied strain rate
     * @param gainFactor the sensitivity by which volume is modified to achieve the stress, or strain, goal
     * @param isStrainRateControlled if this is switched off, then stress is applied and strain is controlled
     */
    GranuHeap(SpeciesType type) : AngleOfRepose_(0)
    {
        const Mdouble initialVolumeFraction = 0.3;
        
        //add species
        auto psd = getPSD(type);
        logger.assert(hasPSD(type),"Material type % not allowed",type);
        auto species = addSingleSpecies(type, getMedianParticleRadius(psd), getMinParticleRadius(psd), *this, true);
        logger(INFO, "%", *species);
        
        baseRadius = 30.0*getMedianParticleRadius(psd);
        fillHeight = 3.0*baseRadius;
        
    
        setName("GranuHeap");
        setGravity({0,0,-9.8}); //set gravity in the system
        
        AxisymmetricIntersectionOfWalls base;
        base.setSpecies(species);
        base.setAxis({0,0,1});
        base.addObject({-1,0,0},{baseRadius,0,0});
        base.addObject({0,0,-1},{0,0,0});
        wallHandler.copyAndAddObject(base);

        AxisymmetricIntersectionOfWalls casing;
        casing.setSpecies(species);
        casing.setAxis({0,0,1});
        casing.addObject({1,0,0},{baseRadius,0,0});
        wallHandler.copyAndAddObject(casing);
   
        setMax({baseRadius,baseRadius,fillHeight});
        setMin({-baseRadius,-baseRadius,0});
        setTimeMax(10);
        setSaveCount(1000);
        
        
        //add particles
        CubeInsertionBoundary insertion;
        BaseParticle particle(speciesHandler.getObject(0));
        insertion.set(&particle, 10000, getMin(), getMax(), Vec3D(0, 0, 0), Vec3D(0, 0, 0), 0, 0);
        insertion.setPSD(getPSD(type));
        insertion.setInitialVolume(100);
        insertion.checkBoundaryBeforeTimeStep(this);
        logger(INFO,"Inserted % particles",particleHandler.getSize());
        
        //write(std::cout, false);
    }

    
//    void computeExternalForces (BaseParticle * CI) override
//    {
//        DPMBase::computeExternalForces (CI);
//        if (fillStage && !CI->isFixed())
//        {
//            // Background damping
//            CI->addForce(-dampingCoefficient * CI->getVelocity());
//            //logger(INFO,"damping coefficient %",dampingCoefficient);
//        }
//    }
    
    
    void printTime() const override {
        logger(INFO, "t % COM % Ene % N %", getTime(), getCentreOfMass().Z, getKineticEnergy()/getElasticEnergy(), particleHandler.getNumberOfObjects());
    }
    
    void actionsAfterTimeStep() override
    {
        //check only every 1000th time
        static unsigned count = 0;
        if (++count>999) count=0; else return;
        
        if (fillStage) {
            if (getKineticEnergy()<1e-4*getElasticEnergy())
            {
                fillStage = false;
                logger(INFO,"remove outer cylinder");
    
                wallHandler.removeLastObject();
                DeletionBoundary boundary;
                boundary.set({0,0,-1},0);
                boundaryHandler.copyAndAddObject(boundary);
            }
        } else {
            if (getKineticEnergy()<1e-5*getElasticEnergy() && timeSwitch) {
                Mdouble COM = getCentreOfMass().Z;
                Mdouble height = 4.0* COM;
                Mdouble angleOfRepose = atan(height/baseRadius)*180.0/pi;
                AngleOfRepose_ = angleOfRepose;
                logger(INFO,"AoR=%, based on conical heap shape",angleOfRepose);
                setTimeMax(getTime()+0.5);
                timeSwitch = false;
            }else if (getKineticEnergy()<1e-5*getElasticEnergy() && !timeSwitch) {
                Mdouble COM = getCentreOfMass().Z;
                Mdouble height = 4.0* COM;
                Mdouble angleOfRepose = atan(height/baseRadius)*180.0/pi;
                AngleOfRepose_ = angleOfRepose;
                logger(INFO,"AoR=%, based on conical heap shape",angleOfRepose);
            }
        }
    }
    
    Mdouble getAngleOfRepose() {return AngleOfRepose_;}
    
    Mdouble baseRadius;
    Mdouble fillHeight;
    bool fillStage = true;
    bool timeSwitch = true;
    Mdouble dampingCoefficient;
private:
    Mdouble AngleOfRepose_;
};

int main(int argc, char* argv[])
{
    //Reading in the five parameters
    SpeciesType speciesType;
    const std::string speciesTypeString = helpers::readFromCommandLine(argc,argv,"-speciesType",std::string("SD100"));
    if (speciesTypeString=="APAPP") speciesType = SpeciesType::APAPP;
    else if (speciesTypeString=="MPT") speciesType = SpeciesType::MPT;
    else if (speciesTypeString=="PH101") speciesType = SpeciesType::PH101;
    else if (speciesTypeString=="SD100") speciesType = SpeciesType::SD100;
    else logger(ERROR, "SpeciesType % not found",speciesType);
    const Mdouble restitutionCoefficient = helpers::readFromCommandLine(argc,argv,"-restitutionCoefficient",0.5);
    const Mdouble relativeCohesionStiffness = helpers::readFromCommandLine(argc,argv,"-relativeCohesionStiffness",0.1);
    const Mdouble slidingFriction = helpers::readFromCommandLine(argc,argv,"-slidingFriction",0.5);
    const Mdouble rollingFriction = helpers::readFromCommandLine(argc,argv,"-rollingFriction",0.1);
    
    
    std::stringstream restartName;
    //restartName << std::fixed << std::setprecision(5) << "AoRFE" << "_e" << restitutionCoefficient << "_rC" << relativeCohesionStiffness << "_mus" << slidingFriction << "_mur" << rollingFriction;
    restartName << std::setprecision(6) << std::scientific << "AoRFE" << "_e" << restitutionCoefficient << "_rC" << relativeCohesionStiffness << "_mus" << slidingFriction << "_mur" << rollingFriction;
    const auto newName = restartName.str();
    std::cout << newName << std::endl;
    
    GranuHeap dpm(speciesType);
    modifySpecies(dpm.speciesHandler.getLastObject(), restitutionCoefficient, relativeCohesionStiffness, slidingFriction, rollingFriction, speciesType);
//    auto species = dynamic_cast<LinearPlasticViscoelasticFrictionSpecies*>(dpm.speciesHandler.getLastObject());
//    dpm.dampingCoefficient = 0.1*species->getDissipation();
    dpm.setName(newName);
    logger(INFO,"Species: %", *dpm.speciesHandler.getLastObject());
    dpm.removeOldFiles();
    dpm.setParticlesWriteVTK(false);
    dpm.setWallsWriteVTK(false);
    dpm.dataFile.setFileType(FileType::NO_FILE);
    dpm.eneFile.setFileType(FileType::NO_FILE);
    dpm.fStatFile.setFileType(FileType::NO_FILE);
    dpm.restartFile.setFileType(FileType::ONE_FILE);
    dpm.restartFile.setSaveCount(100000000);
    dpm.solve();
    
    std::stringstream ss;
    ss << dpm.getAngleOfRepose();
    
    std::cout << "Results:\n" << ss.str();
    helpers::writeToFile(newName+".txt", ss.str());
    
    return 0;
}
