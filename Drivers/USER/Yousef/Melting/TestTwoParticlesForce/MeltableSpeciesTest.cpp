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
#include "Species/MeltableSpecies.h"
#include "Particles/ThermalParticle.h"
#include <Walls/InfiniteWall.h>

/**
 * Simulates the interaction of two particles using the MeltableSpecies
 */
class DPM : public Mercury3D
{
public:
    void setupInitialConditions() override {
        // initial parameters set by the user
        Mdouble timeMax = 0.056; //0.01; //5.0; //0.1;//0.001;
        Mdouble timeStep = 1e-5; //1e-3; //1e-8;
        unsigned int saveCount = 1;
        //TW I set fstatSaveCount equal to savecount, as it produced more output than necessary
        unsigned int fstatSaveCount =  1;
        //TW I set writeVTK to false by default
        bool writeVTK = false;
        Mdouble radius1 = 100e-6; //100e-6;//20e-6;
        Mdouble radius2 = 100e-6; //100e-6;//20e-6;
        Mdouble initialOverlap = 1.0*radius1; //1.0*radius1; //0.03*radius1; //1.6*radius1;
        Mdouble moltenLaterThickness1 = 0.8*radius1; //0.8*radius1; //0.01*radius1;//0.5*radius1;
        Mdouble moltenLaterThickness2 = 0.8*radius1; //0.8*radius2; //0.01*radius2;//0.5*radius2;
        Mdouble initialBondingOverlap = 0.0; //1.6*radius2; //0.0;
        Mdouble ambientTemperature = 155 + 273.15;
        Mdouble meltingTemperature = 178 + 273.15;
        Mdouble initialParticleTemperature = meltingTemperature;//ambientTemperature; //ambientTemperature; //meltingTemperature;
        // mass / density scaling factor
        Mdouble scaleMass = 1e8; //1e8;
        Mdouble density = 1050.0 * scaleMass;
        Mdouble elasticModulus = 1e7; //0.8*2.95e9;
        Mdouble poissonRatio = 0.4;//0.2;
        Mdouble effectiveElasticModulus = 0.5*elasticModulus/((1.0-poissonRatio*poissonRatio));
        Mdouble dissipation = 0.1; //0.0; //0.1; // * sqrt(scaleMass); //0.1 / sqrt(scaleMass); //0.95 / sqrt(scaleMass);
        Mdouble latentHeat = 56.4e3;
        Mdouble heatCapacity = 1200;
        Mdouble thermalConductivityCoefficient = 0.0;
        Mdouble thermalConvectionCoefficient = 0.0;
        Mdouble materialEmissivity = 0.0;
        Mdouble surfaceTension = 0.035 * sqrt(scaleMass); //0.035 * 1e-3 * sqrt(scaleMass);// * 2e-3;
        Mdouble refViscosity = 1e-2 * sqrt(scaleMass); //2e-4 * sqrt(scaleMass);
        Mdouble deltaT = 20.0;
        Mdouble gravity = -9.81/scaleMass;
        //
        Mdouble radiusHarmonicMean = 2.0 * ((radius1*radius2)/(radius1+radius2));
        Mdouble a0 = sqrt(radiusHarmonicMean*initialOverlap);
        
        // now use these parameters to setup a simulation
        setName("MeltableSpeciesTest");
        // delete files from the previous run
        removeOldFiles();
        // set gravitational acceleration
        setGravity(Vec3D(0.0, 0.0,gravity));
        // set domain
        //TW changed domain size, so the contact point is at the origin and the particles are fully visible, now you don't have to correct xballs arguments
        Mdouble R = std::max(radius1, radius2);
        setDomain(Vec3D(-R, -2*R, -R), Vec3D(R, 2*R, R));
        // set simulation time
        setTimeMax(timeMax);
        
        //set species properties
        auto species = speciesHandler.copyAndAddObject(MeltableSpecies());
        species->setDensity(density);
        species->setElasticModulus(elasticModulus);
        species->setPoissonRatio(poissonRatio);
        species->setDissipation(dissipation);
        species->setAmbientTemperature(ambientTemperature);
        species->setMeltingTemperature(meltingTemperature);
        species->setLatentHeat(latentHeat);
        species->setSolidHeatCapacity(heatCapacity);
        species->setThermalConductivityCoefficient(thermalConductivityCoefficient);
        species->setThermalConvectionCoefficient(thermalConvectionCoefficient);
        species->setMaterialEmissivity(materialEmissivity);
        species->setSurfaceTension(surfaceTension);
        species->setRefViscosity(refViscosity);
        species->setDeltaT(deltaT);
        species->setThermalExpansionOnOff(0);
        
        //set time step
        //TW has to be based on the minimum radius
        //species->setMaxTimeStep(std::min(radius1, radius2),density,effectiveElasticModulus);
        setTimeStep(timeStep);//0.3*species->getMaxTimeStep()); //8e-10 , 1e-7, 8e-9
        
        //set output file properties
        setSaveCount(saveCount);
        //fStatFile.setSaveCount(fstatSaveCount);
        //TW now it only saves restart files at the first and last time step, reducing simulation time.
        restartFile.writeFirstAndLastTimeStep();
        setFileType(FileType::ONE_FILE);
        setParticlesWriteVTK(writeVTK);
        setWallsWriteVTK(writeVTK);
        setXBallsAdditionalArguments("-solidf -v0");

        InfiniteWall w0;
        w0.setSpecies(speciesHandler.getObject(0));
        w0.set({0.0, 0.0, -1.0}, {0.0, 0.0, 0.0});
        wallHandler.copyAndAddObject(w0);
        //set properties of 0 particle
/*        ThermalParticle baseP;
        baseP.setSpecies(speciesHandler.getObject(0));
        baseP.setTemperature(initialParticleTemperature);
        baseP.setRadius(radius1);
        //TW I centered the particles, so the contact point is at the origin
        baseP.setPosition(Vec3D(0.0, radius1-0.5*initialOverlap, radius1));
        //TW why add, not set?
        baseP.addMoltenLayerThickness(moltenLaterThickness1);
        auto baseP0 = particleHandler.copyAndAddObject(baseP);*/

        //set properties of 1st particle
        ThermalParticle p;
        p.setSpecies(speciesHandler.getObject(0));
        p.setTemperature(initialParticleTemperature);
        p.setRadius(radius1);
        //TW I centered the particles, so the contact point is at the origin
        p.setPosition(Vec3D(0.0, radius1-0.5*initialOverlap, radius1)); //+2.0*radius1));
        //TW why add, not set?
        p.addMoltenLayerThickness(moltenLaterThickness1);
        auto p0 = particleHandler.copyAndAddObject(p);
    
        //set properties of 2nd particle
        ThermalParticle I;
        I.setSpecies(speciesHandler.getObject(0));
        I.setTemperature(initialParticleTemperature);
        I.setRadius(radius2);
        // Next to eachOther:
        I.setPosition(Vec3D(0.0, -radius2+0.5*initialOverlap, radius2));
        // Top:
        //I.setPosition(Vec3D(0.0, radius2-0.5*initialOverlap, (radius1+2.0*radius2)-initialOverlap)); //4.0*radius2)-initialOverlap));
        I.setVelocity(Vec3D(0.0, 0.0, 0.0));
        I.addMoltenLayerThickness(moltenLaterThickness2);
        auto p1 = particleHandler.copyAndAddObject(I);

        // set initial bonding overlap
        // get pointer to the interaction
        auto i = interactionHandler.getInteraction(p0,p1,1);
        auto mi = dynamic_cast<MeltableInteraction*>(i);
        logger.assert_always(mi!= nullptr,"mi not set");
        //TW why do you need to compute the force here?
        //mi->computeNormalForce();
        mi->addBondingOverlap(initialBondingOverlap);
        //
        //Mdouble a0 = sqrt(radiusHarmonicMean*mi->getOverlap());
        mi->setContactRadiusMeltable(a0);
    }
    
    // Write custom output to the ene file
    void writeEneTimeStep(std::ostream &os) const override
    {
        //write header line
        if (eneFile.getCounter() == 1) os << "time relativeVelocity tensionForce visForce bondingForce dampingForce elasticForce "
                                             "overlap solidOverlap bondingOverlap meltRateP meltRateI "
                                             "contactRadius solidContactRadius bondingContactRadius visCoeff visTimeStep timeStep contactPoint NormalForce\n";
        //write interaction properties at every time step
        auto i = interactionHandler.getObject(0); //getLastObject() -> P-W
        logger.assert_debug(i,"No interaction exists");
        auto mi = dynamic_cast<const MeltableInteraction *>(i);
        logger.assert_debug(mi,"Interaction is not of type MeltableInteraction");
        os <<  getTime() << ' '
            << mi->getNormalRelativeVelocity() << ' '
            << mi->getTensionForce() << ' '
            << mi->getViscousForce() << ' '
            << mi->getBondingForce() << ' '
            << mi->getDampingForce() << ' '
            << mi->getElasticForce() << ' '
            << mi->getOverlap() << ' '
            << mi->getSolidOverlap() << ' '
            << mi->getBondingOverlap() << ' '
            << mi->getMoltenLayerThicknessRateP() << ' '
            << mi->getMoltenLayerThicknessRateI() << ' '
            << mi->getContactRadiusMeltable() << ' '
            << mi->getSolidContactRadius() << ' '
            << mi->getBondingContactRadius() << ' '
            << mi->getVisCoeff() << ' '
            << mi->getHandler()->getDPMBase()->getTimeStep() << ' '
            << mi->getContactPoint() << ' '
            << mi->getAbsoluteNormalForce() << std::endl;
        //alternative to plot all:
        //os << getTime() << ' ' << *mi << std::endl;
        //You can plot a force balance in gnuplot using e.g.
        //p 'MeltableSpeciesTest.ene' u 1:3 t 'tension force', '' u 1:4 t 'viscous force', '' u 1:6 t 'elastic force'
    }
    
    //Write to the console
    void printTime() const override {
        //TW printTime outputted a lot to the screen, which is costly as well.
        //TW Instead of screen output I now write a custom ene file, see writeEneTimeStep.
        logger(INFO,"Time % TimeMax %", getTime(), getTimeMax());
    }
};

int main()
{
    DPM dpm;
    dpm.solve();
    return 0;
}
