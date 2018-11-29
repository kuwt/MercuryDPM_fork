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
#include "Boundaries/LeesEdwardsBoundary.h"
#include "addSpecies.h"

/**
 * This simulation consists of two parts:
 * A)
 * Particles, with given volume fracion, are
 * placed in a cubic box of dimensions Xmax, Ymax, Zmax
 * and period boundaries are applied in all directions.
 * With a given strain-rate, all 6 period boundaries are moved inwards
 * isotropically to compress the particle gas to target normal stresses.
 * B)
 * After the system is relaxed, it is sheared without stress control.
 * (doesn't work yet)
 **/
class ConstantRestitutionTest : public Mercury3D
{
public:
    
    /**
     *
     * @param stressGoal the stress tensor that the material should be subjected to
     * @param strainRate the applied strain rate
     * @param gainFactor the sensitivity by which volume is modified to achieve the stress, or strain, goal
     * @param isStrainRateControlled if this is switched off, then stress is applied and strain is controlled
     */
    ConstantRestitutionTest(SpeciesType type)
    {
        const Mdouble initialVolumeFraction = 0.3;
        
        //add species
        auto psd = getPSD(type);
        logger.assert(hasPSD(type),"Material type % not allowed",type);
        auto species = addSingleSpecies(type, getMedianParticleRadius(psd), getMinParticleRadius(psd), *this, true);
        logger(INFO, "%", *species);
    
        setName("ConstantRestitutionAPAPSelfTest");
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
        casing.setVTKVisibility(false);
        wallHandler.copyAndAddObject(casing);
   
        setMax({baseRadius,baseRadius,fillHeight});
        setMin({-baseRadius,-baseRadius,0});
        setTimeMax(0.2);
        setSaveCount(200);
        
        //add particles
        CubeInsertionBoundary insertion;
        BaseParticle particle(speciesHandler.getObject(0));
        insertion.set(&particle, 10000, Vec3D(getXMin(),getYMin(),getZMax()), getMax(), Vec3D(0, 0, 0), Vec3D(0, 0, 0), 0, 0);
        insertion.setPSD(getPSD(type));
        insertion.setInitialVolume(100);
        insertion.checkBoundaryBeforeTimeStep(this);
        logger(INFO,"Inserted % particles",particleHandler.getSize());
        
        //write(std::cout, false);
        setParticlesWriteVTK(true);
        setWallsWriteVTK(true);
    }
    
    void printTime() const override {
        logger(INFO, "t % COM % Ene % N %", getTime(), getCentreOfMass().Z, getKineticEnergy()/getElasticEnergy(), particleHandler.getNumberOfObjects());
    }
    
    void actionsAfterTimeStep() override
    {
        //check only every 1000th time
        static unsigned count = 0;
        if (++count>999) count=0; else return;
        
        if (fillStage) {
            if (getKineticEnergy()<1e-5*getElasticEnergy())
            {
                fillStage = false;
                logger(INFO,"remove outer cylinder");
    
                wallHandler.removeLastObject();
                DeletionBoundary boundary;
                boundary.set({0,0,-1},0);
                boundaryHandler.copyAndAddObject(boundary);
            }
        } else {
            if (getKineticEnergy()<1e-5*getElasticEnergy()) {
                logger(INFO,"system is relaxed");
                Mdouble COM = getCentreOfMass().Z;
                Mdouble height = 4.0* COM;
                Mdouble angleOfRepose = atan(height/baseRadius)*180.0/pi;
                logger(INFO,"AoR=%, based on conical heap shape",angleOfRepose);
                setTimeMax(getTime());
            }
        }
    }
    
    const Mdouble baseRadius = 20e-3;
    const Mdouble fillHeight = 5.0*baseRadius;
    bool fillStage = true;
};

int main(int argc UNUSED, char* argv[] UNUSED)
{
    ConstantRestitutionTest dpm(SpeciesType::APAPP);
    auto species = dynamic_cast<LinearPlasticViscoelasticFrictionSpecies*>(dpm.speciesHandler.getLastObject());
    species->setCohesionStiffness(0.1*species->getCohesionStiffness());
    //species->setDissipation(0);
    //species->setUnloadingStiffnessMax(species->getLoadingStiffness());
    logger(INFO,"Bo %",species->computeBondNumberMax(2.5e-3,9.8));
    //dpm.removeOldFiles();
    //dpm.setParticlesWriteVTK(true);
    //dpm.setWallsWriteVTK(true);
    dpm.solve();
    return 0;
}
