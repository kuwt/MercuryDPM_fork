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

//based on /storage2/usr/people/sluding/MDCC/C3DshearXL30/MU0_LONG2
#include <Species/LinearViscoelasticFrictionLiquidMigrationWilletSpecies.h>
#include <Species/LinearViscoelasticSpecies.h>
#include "Particles/LiquidFilmParticle.h"
#include "Mercury3D.h"
#include "Boundaries/AngledPeriodicBoundary.h"
#include "Walls/AxisymmetricIntersectionOfWalls.h"
#include <iomanip>

class ShearCell3D : public Mercury3D {
public:

    ShearCell3D (Mdouble splitRadius, Mdouble angularVelocity)
    {
        setName("ShearCell3DInitialConditions");
        readRestartFile();
        setRestarted(false);
        setName("ShearCell3D_vol_20nl_surf1_30degree");
        setSaveCount(2000);

        // get pointer to species loaded from initial conditions
        auto originalSpecies = dynamic_cast<LinearViscoelasticSpecies*>(speciesHandler.getObject(0));

        // create new species that includes frictional properties
        species = new LinearViscoelasticFrictionLiquidMigrationWilletSpecies();
        species->setDensity(originalSpecies->getDensity());
        species->setStiffness(originalSpecies->getStiffness());
        species->setDissipation(originalSpecies->getDissipation());

		//default species parameters
        species->setSlidingDissipation(0.5e-3);
        species->setSlidingStiffness(120);
	    species->setSlidingFrictionCoefficient(0.01);
        species->setLiquidBridgeVolumeMax(40e-12);
		species->setSurfaceTension(0.020);
		species->setContactAngle(20.0*constants::pi/180.0);

        // replace original species with the new one
        speciesHandler.clear();
        speciesHandler.addObject(species);
        for (BaseParticle* p : particleHandler)
            p->setSpecies(species);
        for (BaseWall* w : wallHandler)
            w->setSpecies(species);

        unsigned int i = 0;
        for (BaseParticle* p : particleHandler)
        {
            if (p->isFixed() && mathsFunc::square(p->getPosition().X)+mathsFunc::square(p->getPosition().Y)>mathsFunc::square(splitRadius))
            {
                MPIndex.push_back(p->getIndex());
                Vec3D pos = p->getPosition();
                MPRadius.push_back(sqrt(pos.X * pos.X + pos.Y * pos.Y));
                MPAngle.push_back(atan(pos.Y / pos.X));
                MPHeight.push_back(atan(pos.Z));
                i++;
            }
        }
        AngularVelocity = angularVelocity;
        AngledPeriodicBoundary* b = dynamic_cast<AngledPeriodicBoundary*>(boundaryHandler.getObject(0));
        openingAngle = b->getOpeningAngle();
        //write(std::cout,false);
        std::cout << "openingAngle " << openingAngle << std::endl;
        
        boundaryHandler.clear();
    }

private:
    /////////////////////////////////////////////////////////////////////////////////////////////////////
    ///This is were the walls are implemented
    /////////////////////////////////////////////////////////////////////////////////////////////////////
    void setupInitialConditions()
    {
    }

    void printTime() const
    {
		Mdouble volTot = 0.0;
		for (BaseParticle* p : particleHandler)
		{
			LiquidFilmParticle* l = dynamic_cast<LiquidFilmParticle*>(p);
			volTot += l->getLiquidVolume();
			
		}
		for (BaseInteraction* p : interactionHandler)
		{
			LiquidMigrationWilletInteraction* l = dynamic_cast<LiquidMigrationWilletInteraction*>(p);
			volTot += l->getLiquidBridgeVolume();
			
		}
		
		
		
		
		std::cout 
		<< "t=" << std::setprecision(3) << std::left << std::setw(6) << getTime()
        << ", tmax=" << std::setprecision(3) << std::left << std::setw(6) << getTimeMax()
        << ", VolTot=" << std::setprecision(3) << std::left << std::setw(6) << volTot
        << std::endl;
		std::cout.flush();
        //std::cout
            //<< "t=" << getTime()
            //<< " rotationAngleDeg " << getTime()*AngularVelocity*180.0/2.0/ constants::pi
            //<< " Ene " << getKineticEnergy()/getElasticEnergy()
            //<< std::endl;
            Mdouble volTotP = 0.0;
        unsigned int nLB = 0;
        for (BaseParticle * const p : particleHandler)
        {
            LiquidFilmParticle* l = dynamic_cast<LiquidFilmParticle*>(p);
            //             if  (l->getLiquidVolume())
            //                 std::cout << " P" << l->getIndex() << " " << l->getLiquidVolume() <<std::endl;
            if (l == nullptr)
                std::cerr << "Error" << std::endl;
            volTotP += l->getLiquidVolume();
        }
        Mdouble volTotI = 0.0;
        for (BaseInteraction* i : interactionHandler)
        {
            LiquidMigrationWilletInteraction* l = dynamic_cast<LiquidMigrationWilletInteraction*>(i);
            if (l == nullptr)
                std::cerr << "Error" << std::endl;
            //            if  (l->getLiquidBridgeVolume())
            //                std::cout << " I" << l->getIndex() << " " << l->getLiquidBridgeVolume();
            volTotI += l->getLiquidBridgeVolume();
            if (l->getLiquidBridgeVolume() != 0.0)
                nLB++;
        }
        std::cout
            << "t=" << std::setprecision(6) << std::left << std::setw(6) << getTime()
            //<< ", tmax="    << std::setprecision(3) << std::left << std::setw(6) << getTimeMax()
            << ", Vol=" << std::setprecision(8) << std::left << std::setw(12) << volTotP + volTotI
            << "(P=" << std::setprecision(6) << std::left << std::setw(12) << volTotP
            << "), nLB=" << std::setprecision(6) << std::left << std::setw(12) << nLB
            << std::endl;
        std::cout.flush();
    }

    void actionsAfterTimeStep () {
        Mdouble dAlpha = getTime()*AngularVelocity;
        for (unsigned int i=0; i<MPIndex.size(); i++) {
            Mdouble a = MPAngle[i]+dAlpha;
            if (a> openingAngle) MPAngle[i] -= openingAngle;
            //if (a<=0.0) MPAngle[i] += openingAngle;
            Mdouble s = MPRadius[i]*sin(a);
            Mdouble c = MPRadius[i]*cos(a);
            BaseParticle* p = particleHandler.getObject(MPIndex[i]);
            p->setPosition
                (Vec3D(c, s, MPHeight[i] ));
            p->setVelocity
                (Vec3D(-AngularVelocity*s, AngularVelocity*c, 0.0));
            p->setAngularVelocity
                (Vec3D(0.0,0.0,AngularVelocity));
        }
    }

    Mdouble openingAngle;
    Mdouble AngularVelocity;
    std::vector<unsigned int> MPIndex;
    std::vector<Mdouble> MPRadius;
    std::vector<Mdouble> MPAngle;
    std::vector<Mdouble> MPHeight;
public:
    LinearViscoelasticFrictionLiquidMigrationWilletSpecies* species;
    
};


int main() {

    Mdouble splitRadius = 85e-3;
    Mdouble angularVelocity = 0.1*2.0*constants::pi;
    ShearCell3D SC(splitRadius, angularVelocity);
	
    for (BaseParticle* p : SC.particleHandler)
    {
		LiquidFilmParticle* l = dynamic_cast<LiquidFilmParticle*>(p);
		l->setLiquidVolume(20e-12);
	}

	//SC.species->setSlidingFrictionCoefficient(0.01); //Change back to 0.01; for one simulation to 0.5
	SC.species->setSlidingDissipation(0.5e-3);
	SC.species->setSlidingStiffness(120);
	SC.species->setSlidingFrictionCoefficient(0.01);
    SC.species->setLiquidBridgeVolumeMax(40e-12);
	SC.species->setSurfaceTension(0.020);
	SC.species->setContactAngle(20.0*constants::pi/180.0);
    SC.setXBallsAdditionalArguments("-v0 -solidf");
    SC.setHGridMaxLevels(2);
    SC.setTimeMax(0.7);
    SC.dataFile.setFileType(FileType::MULTIPLE_FILES_PADDED);
    SC.fStatFile.setFileType(FileType::MULTIPLE_FILES_PADDED);
    SC.restartFile.setFileType(FileType::MULTIPLE_FILES_PADDED);
    SC.setSaveCount(500);
    SC.solve();
    return 0;
}
