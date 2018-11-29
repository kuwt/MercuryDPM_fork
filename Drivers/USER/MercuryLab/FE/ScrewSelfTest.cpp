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
#include <Walls/Screw.h>
#include "Boundaries/LeesEdwardsBoundary.h"
#include "addSpecies.h"

/**
 **/
class ScrewSelfTest : public Mercury3D
{
public:

    ScrewSelfTest()
    {
        //basic properties
        setName("ScrewSelfTest");
        setGravity({0,0,-9.8}); //set gravity in the system

        //add species
        SpeciesType type = SpeciesType::APAPP;
        auto psd = getPSD(type);
        logger.assert(hasPSD(type),"Material type % not allowed",type);
        auto species = addSingleSpecies(type, getMedianParticleRadius(psd), getMinParticleRadius(psd), *this, true);

        // define screw
        Mdouble casingRadius = 35e-3; //< inner radius of the casing
        Mdouble shaftRadius = 15e-3; //< outer radius of the shaft
        Mdouble bladeRadius = casingRadius - 1.5e-3; //< outer radius of the blades
        Mdouble bladeDistance = 30e-3; //< distance between the blades
        Mdouble bladeInclination = 20.0; //< angle by which the blade normal is inclined with respect to the x-axis
        Mdouble bladeOpening = 60.0; //< arc angle of the slice
        Mdouble bladeThickness = 2.5e-3; //< width of the blade
        Mdouble screwPitch = 20e-3; //offset between the screw blades at the in- and outflow
        Mdouble inflowScrewLength = 2*screwPitch; //length of the inflow screw
        Mdouble rpm = 0;
        Screw inflowScrew({0,0,0}, inflowScrewLength, bladeRadius, -inflowScrewLength/screwPitch, rpm/60.0, bladeThickness/2.0, ScrewType::singleHelix);
        inflowScrew.setSpecies(species); //make the helix rubber-like to reduce backsplatter
        inflowScrew.setPosition({0,0,0});
        inflowScrew.setOrientationViaNormal({1,0,0});
        wallHandler.copyAndAddObject(inflowScrew);

        CubeInsertionBoundary insertionBoundary;
        BaseParticle particle(speciesHandler.getObject(0));
        Vec3D pos = Vec3D(0.625*inflowScrewLength,0,2.0*casingRadius);
        insertionBoundary.set(&particle, 100, pos, pos, Vec3D(0, 0, 0), Vec3D(0, 0, 0), 1, 1);
        insertionBoundary.setPSD(getPSD(type));
        boundaryHandler.copyAndAddObject(insertionBoundary);

        setDomain({0,-casingRadius,-casingRadius},{inflowScrewLength,casingRadius,casingRadius});
        setTimeMax(0.2);
        setSaveCount(200);
        setParticlesWriteVTK(true);
        setWallsWriteVTK(true);
    }
    
};

int main()
{
    ScrewSelfTest dpm;
    dpm.solve();
    return 0;
}
