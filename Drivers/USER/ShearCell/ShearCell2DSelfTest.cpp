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

#include "Species/LinearViscoelasticFrictionLiquidMigrationWilletSpecies.h"
#include "Mercury2D.h"
#include "Boundaries/AngledPeriodicBoundary.h"
#include "Particles/LiquidFilmParticle.h"
#include <iomanip>
#include <CG/Fields/LiquidMigrationFields.h>
#include <CG/TimeAveragedCG.h>
#include "CG/CG.h"

class ShearCell2D : public Mercury2D
{
public:

    ShearCell2D(Mdouble OuterRadius_, Mdouble InnerRadius_, Mdouble AngularVelocity_)
        : OuterRadius(OuterRadius_), InnerRadius(InnerRadius_), AngularVelocity(AngularVelocity_)
    {
        species = speciesHandler.copyAndAddObject(LinearViscoelasticFrictionLiquidMigrationWilletSpecies());

        //check input variables
        if (InnerRadius < 0)
        {
            std::cerr << "error in constructor of ShearCell2D; InnerRadius has to be positive" << std::endl;
            exit(-1);
        }
        else if (InnerRadius > OuterRadius)
        {
            std::cerr << "error in constructor of ShearCell2D; InnerRadius has to be smaller than OuterRadius" << std::endl;
            exit(-1);
        }

        //set default name
        setName("ShearCell2DSelfTest");
        //unset gravity
        setGravity(Vec3D(0.0, 0.0, 0.0));

        //use particles of unit mass and diameter 
        species->setDensity(6.0 / constants::pi);

        //begin: set geometry

        //set domain size
        setXMin(0.0);
        setYMin(0.0);
        setZMin(-0.5);
        setXMax(OuterRadius);
        setYMax(getXMax());
        setZMax(0.5);

        //set boundaries
        AngledPeriodicBoundary B0;
        Vec3D normal_left(1.0, 0.0, 0.0);
        Vec3D normal_right(0.0, -1.0, 0.0);
        Vec3D origin(0.0, 0.0, 0.0);
        B0.set(normal_left, normal_right, origin);
        boundaryHandler.copyAndAddObject(B0);

        //set walls
        LiquidFilmParticle F0;
        //F0.setLiquidVolume(2e-11);
        F0.setSpecies(speciesHandler.getObject(0));
        F0.setVelocity(Vec3D(0.0, 0.0, 0.0));
        //F0.setOrientation(Vec3D(0.0, 0.0, 0.0));
        F0.setAngularVelocity(Vec3D(0.0, 0.0, 0.0));
        F0.setRadius(0.5);
        F0.fixParticle();

        Mdouble Circumference = OuterRadius * 0.5 * constants::pi;
        Mdouble N = std::floor(Circumference / (2.0 * F0.getRadius()));
        Mdouble dAlpha = 0.5 * constants::pi / (double)N;
        for (Mdouble Alpha = dAlpha / 2.0; Alpha < 0.5 * constants::pi; Alpha += dAlpha)
        {
            F0.setPosition(Vec3D(OuterRadius * mathsFunc::sin(Alpha), OuterRadius * mathsFunc::cos(Alpha), 0.0));
            particleHandler.copyAndAddObject(F0);
        }

        Circumference = InnerRadius * 0.5 * constants::pi;
        N = std::floor(Circumference / (2.0 * F0.getRadius()));
        dAlpha = 0.5 * constants::pi / (double)N;
        for (Mdouble Alpha = dAlpha / 2.0; Alpha < 0.5 * constants::pi; Alpha += dAlpha)
        {
            F0.setPosition(Vec3D(InnerRadius * mathsFunc::sin(Alpha), InnerRadius * mathsFunc::cos(Alpha), 0.0));
            particleHandler.copyAndAddObject(F0);
        }

        //end: set geometry
        particleHandler.setStorageCapacity(200);
        writeRestartFile();
    }


private:
    /////////////////////////////////////////////////////////////////////////////////////////////////////
    ///This is were the walls are implemented
    /////////////////////////////////////////////////////////////////////////////////////////////////////

    void setupInitialConditions()
    {
        //introduce particles randomly
        LiquidFilmParticle P0;
        P0.setLiquidVolume(2e-11);
        P0.setSpecies(speciesHandler.getObject(0));
        P0.setVelocity(Vec3D(0.0, 0.0, 0.0));
        //P0.setOrientation(Vec3D(0.0, 0.0, 0.0));
        P0.setAngularVelocity(Vec3D(0.0, 0.0, 0.0));
        P0.setRadius(0.5);

        Mdouble dRadius = 1.8 * P0.getRadius();
        //minimum distance between particles
        Mdouble NRadial = std::floor((OuterRadius - InnerRadius) / dRadius);
        Mdouble dRadiusRadial = (OuterRadius - InnerRadius) / NRadial;
        //actual radial distance between particles
        for (int i = 1; i < NRadial; i++)
        {
            Mdouble Radius = OuterRadius - i*dRadiusRadial;
            Mdouble Circumference = Radius * 0.5 * constants::pi;
            Mdouble N = std::floor(Circumference / dRadius);
            Mdouble dAlpha = 0.5 * constants::pi / (double)N;
            //actual angular distance between particles
            for (Mdouble Alpha = dAlpha / 4.0 + dAlpha / 2.0 * (double)((int)N % 2); Alpha < 0.5 * constants::pi; Alpha += dAlpha)
            {
                P0.setPosition(Vec3D(Radius * mathsFunc::sin(Alpha), Radius * mathsFunc::cos(Alpha), 0.0));
                particleHandler.copyAndAddObject(P0);
            }
        }
        setHGridMaxLevels(1);
        //setHGridNumberOfBucketsToPower(particleHandler.getNumberOfObjects());
        write(std::cout, false);
    }

    void actionsAfterTimeStep()
    {
        Mdouble Circumference = OuterRadius * 0.5 * constants::pi;
        Mdouble N = std::floor(Circumference / (2.0 * particleHandler.getObject(0)->getRadius()));
        Mdouble dAlpha = 0.5 * constants::pi / (double)N;
        for (unsigned int i = 0; i < N; i++)
        {
            Mdouble Alpha = dAlpha * (0.5 + (double)i) + getTime() * AngularVelocity;
            particleHandler.getObject(i)->setPosition(Vec3D(OuterRadius * mathsFunc::sin(Alpha), OuterRadius * mathsFunc::cos(Alpha), 0.0));
        }
        //if (getTime()>10.6)
        //    setSaveCount(1);

    }

    void outputXBallsData(std::ostream& os) const
    {
        os << particleHandler.getNumberOfObjects()*2.0 + interactionHandler.getNumberOfObjects()
            << " " << getTime()
            << " " << getXMin()
            << " " << getYMin()
            << " " << getZMin()
            << " " << getXMax()
            << " " << getYMax()
            << " " << getZMax()
            << " " << std::endl;

        // This outputs the particle data
        for (unsigned int i = 0; i < particleHandler.getNumberOfObjects(); i++)
        {
            outputXBallsDataParticle(i, 14, os);
        }
        // This outputs the particle data
        for (auto i : particleHandler)
        {
            LiquidFilmParticle* j = dynamic_cast<LiquidFilmParticle*>(i);
            os
                << j->getPosition().X << " "
                << j->getPosition().Y << " "
                << j->getPosition().Z + 1 << " "
                << 100*(j->getLiquidVolume() > 0) << " 0 0 "
                << sqrt(j->getLiquidVolume() / 2e-11)*0.5 / 5
                << " 0 0 0 0 0 0 0" << std::endl;
        }
        // This outputs the particle data
        for (auto i : interactionHandler)
        {
            LiquidMigrationWilletInteraction* j = dynamic_cast<LiquidMigrationWilletInteraction*>(i);
            os
                << j->getContactPoint().X << " "
                << j->getContactPoint().Y << " "
                << j->getContactPoint().Z + 1 << " "
                << 100*(j->getLiquidBridgeVolume() > 0) << " 0 0 "
                << sqrt(j->getLiquidBridgeVolume() / 2e-11)*0.5 / 5
                << " 0 0 0 0 0 0 0" << std::endl;
        }
#ifdef DEBUG_OUTPUT
        std::cerr << "Have output the properties of the problem to disk " << std::endl;
#endif
    }

    void printTime() const
    {
        Mdouble volTotP = 0.0;
        unsigned int nLB = 0;
        for (BaseParticle * const p : particleHandler)
        {
            LiquidFilmParticle* l = dynamic_cast<LiquidFilmParticle*>(p);
            logger.assert(l,"Error in outputXBallsData");
            volTotP += l->getLiquidVolume();
        }
        Mdouble volTotI = 0.0;
        for (BaseInteraction* i : interactionHandler)
        {
            LiquidMigrationWilletInteraction* l = dynamic_cast<LiquidMigrationWilletInteraction*>(i);
            logger.assert(l,"Error in outputXBallsData");
            volTotI += l->getLiquidBridgeVolume();
            if (l->getLiquidBridgeVolume() != 0.0)
                nLB++;
        }
        logger(DEFAULT,"t=%, Vol=% (bridges: % film: %), #LB=%",getTime(),volTotP+volTotI, volTotI, volTotP, nLB);
    }

    Mdouble OuterRadius;
    Mdouble InnerRadius;
    Mdouble AngularVelocity;
public:
    LinearViscoelasticFrictionLiquidMigrationWilletSpecies* species;
};

int main()
{
    ShearCell2D SC(10.0, 5.5, 0.01 * 2.0 * constants::pi);
    SC.setSystemDimensions(3);
    SC.setXBallsAdditionalArguments("-v0 -solid -3dturn 1");
    
    SC.setFileType(FileType::ONE_FILE);
    SC.setTime(0.0);
    SC.species->setStiffness(1000.0);
    SC.species->setDissipation(25.0); //dissipation_^2<<4km
    SC.species->setLiquidBridgeVolumeMax(2e-11);
    SC.species->setDistributionCoefficient(0.8);
    SC.species->setSurfaceTension(0.020);
    SC.species->setContactAngle(20.0 * constants::pi / 180.0);

    SC.setTimeStep(1e-3);
    SC.setTimeMax(2.0);
    //SC.setTimeMax(0.0);
    //comment the next line to get output
    SC.setSaveCount(10);
    //uncomment the next lines to get VTK output
    //SC.setParticlesWriteVTK(true);
    //SC.interactionHandler.setWriteVTK(FileType::MULTIPLE_FILES);

    //define object for global coarse-graining
    BaseCG* cg0 = SC.cgHandler.copyAndAddObject(CG<CGCoordinates::O,CGFunctions::Lucy,CGFields::LiquidMigrationFields>());
    cg0->statFile.setSaveCount(10);
    //gnu file to visualise the cg data
    helpers::writeToFile("ShearCell2DSelfTest.1D.gnu","p 'ShearCell2DSelfTest.0.stat' u 1:2, '' u 1:3, '' u 1:($2+$3)");

    //define object for local coarse-graining of the liquid properties
    BaseCG* cg1 = SC.cgHandler.copyAndAddObject(TimeAveragedCG<CGCoordinates::XY,CGFunctions::Heaviside,CGFields::LiquidMigrationFields>());
    cg1->statFile.setSaveCount(10);
    cg1->setN(30);
    cg1->setWidth(0.25);
    //matlab file to plot the cg output
    helpers::writeToFile("ShearCell2DSelfTest.m","raw  = importdata('ShearCell2DSelfTest.1.stat',' ',2);\n"
     "x = raw.data(:,2);\n"
     "nx = length(unique(x))\n"
     "ny = length(x)/nx\n"
     "x = reshape(raw.data(:,2),[ny,nx]);\n"
     "y = reshape(raw.data(:,3),[ny,nx]);\n"
     "b = reshape(raw.data(:,4),[ny,nx]);\n"
     "v = reshape(raw.data(:,5),[ny,nx]);\n"
     "figure(1)\n"
     "contourf(x,y,b)\n"
     "xlabel('x')\n"
     "ylabel('y')\n"
     "title('Density of liquid bridge volume')");

    //define object to test radial coarse-graining
    BaseCG* cg2 = SC.cgHandler.copyAndAddObject(CG<CGCoordinates::R,CGFunctions::Gauss>());
    cg2->statFile.setSaveCount(10);
    cg2->setN(100);
    cg2->setWidth(0.25);
    cg2->setTimeMin(SC.getTimeMax());
    //gnu file to visualise the cg data
    helpers::writeToFile("ShearCell2DSelfTest.2D.gnu","p 'ShearCell2DSelfTest.2.stat' u 2:3 w l t 'bridge volume', '' u 2:4 w l t 'film volume ', '' u 2:($3+$4) w l t 'total volume'");

    SC.solve();

    return 0;
}
