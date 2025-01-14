//Copyright (c) 2013-2023, The MercuryDPM Developers Team. All rights reserved.
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

#include <Mercury3D.h>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <cmath>
#include "Species/LinearViscoelasticFrictionSpecies.h"
#include "Walls/AxisymmetricIntersectionOfWalls.h"
#include "Walls/IntersectionOfWalls.h"
#include "Walls/InfiniteWall.h"
#include "Boundaries/PeriodicBoundary.h"
#include <chrono>
/*
* This is our problem description. Everything is set up here.
* We inherit from Mercury3D, since this gives full flexibility.
* For more predefined problems (for instance, chutes), please see the
* documentation.
*/
class RotatingDrum : public Mercury3D
{
public:
    /* We define our own 'setupInitialConditions' function here,
     * which defines all the specifics about our simulation here.
     */
    void setupInitialConditions() override
    {
        //The first step: set any properties which are always true for your system.
        // (for instance, if gravity remains constant, set it here)
        fractionalPolydispersity = 0.0;
        radiusS1 = 0.0015; // 3mm diameter
        radiusS2 = radiusS1*sizeRatio;

        rhoS1 = 2500.0;
        rhoS2 = densityRatio*rhoS1;

        massS1 = 4 / 3 * constants::pi * pow(radiusS1, 3.0) * rhoS1;
        massS2 = 4 / 3 * constants::pi * pow(radiusS2, 3.0) * rhoS2;

        double fillVolume = drumFillFraction*4*drumRadius*0.024;
        //numS1 = volumeFraction*fillVolume/(4./3. * constants::pi*pow(radiusS1,3.0));
        //numS2 = (1. - volumeFraction)*fillVolume/(4./3. * constants::pi*pow(radiusS2,3.0));
        numS1 = 10; // for aquick test
        numS2 = 10 ; // for a quick test

        numS1ToBeInserted = numS1;
        numS2ToBeInserted = numS2;

        tc = 1 / 2000.0;
        //Now, decide what Species you need for your system.
        speciesHandler.clear();

        auto speciesDrum = speciesHandler.copyAndAddObject(LinearViscoelasticFrictionSpecies());
        auto speciesS1 =  speciesHandler.copyAndAddObject(LinearViscoelasticFrictionSpecies());
        auto speciesS2 =  speciesHandler.copyAndAddObject(LinearViscoelasticFrictionSpecies());

        double RPSInitial = 0.0;

        speciesDrum->setDensity(rhoS1);
        speciesDrum->setCollisionTimeAndRestitutionCoefficient(tc, CORDrum, massS1);

        speciesDrum->setSlidingDissipation(speciesDrum->getDissipation()*2./7.);
        speciesDrum->setSlidingStiffness(speciesDrum->getStiffness()*2./7.);
        speciesDrum->setSlidingFrictionCoefficient(slidingFrictionDrum);

        speciesDrum->setRollingStiffness(speciesDrum->getStiffness()*2.0/7.0);
        speciesDrum->setRollingFrictionCoefficient(rollingFrictionDrum);
        speciesDrum->setRollingDissipation(speciesDrum->getDissipation()*2./7.);

        speciesDrum->setTorsionStiffness(speciesDrum->getStiffness()*2.0/7.0);
        speciesDrum->setTorsionFrictionCoefficient(torsionFrictionDrum);
        speciesDrum->setTorsionDissipation(speciesDrum->getDissipation()*2./7.);
        //

        //
        speciesS1->setDensity(rhoS1);
        speciesS1->setCollisionTimeAndRestitutionCoefficient(tc, CORS1, massS1);

        speciesS1->setSlidingDissipation(speciesS1->getDissipation()*2./7.);
        speciesS1->setSlidingStiffness(speciesS1->getStiffness()*2./7.);
        speciesS1->setSlidingFrictionCoefficient(slidingFriction1);

        speciesS1->setRollingStiffness(speciesS1->getStiffness()*2.0/7.0);
        speciesS1->setRollingFrictionCoefficient(rollingFriction1);
        speciesS1->setRollingDissipation(speciesS1->getDissipation()*2./7.);

        speciesS1->setTorsionStiffness(speciesS1->getStiffness()*2.0/7.0);
        speciesS1->setTorsionFrictionCoefficient(torsionFriction1);
        speciesS1->setTorsionDissipation(speciesS1->getDissipation()*2./7.);
        //

        speciesS2->setDensity(rhoS2);
        speciesS2->setCollisionTimeAndRestitutionCoefficient(tc, CORS2, massS2);

        speciesS2->setSlidingDissipation(speciesS2->getDissipation()*2./7.);
        speciesS2->setSlidingStiffness(speciesS2->getStiffness()*2./7.);
        speciesS2->setSlidingFrictionCoefficient(slidingFriction2);

        speciesS2->setRollingStiffness(speciesS2->getStiffness()*2.0/7.0);
        speciesS2->setRollingFrictionCoefficient(rollingFriction2);
        speciesS2->setRollingDissipation(speciesS2->getDissipation()*2./7.);

        speciesS2->setTorsionStiffness(speciesS2->getStiffness()*2.0/7.0);
        speciesS2->setTorsionFrictionCoefficient(torsionFriction2);
        speciesS2->setTorsionDissipation(speciesS2->getDissipation()*2./7.);

        auto speciesDrumAndS1 = speciesHandler.getMixedObject(speciesDrum,speciesS1);

        speciesDrumAndS1->setCollisionTimeAndRestitutionCoefficient(tc, ((CORS1 + CORDrum) / 2) , massS1, massS1);

        speciesDrumAndS1->setSlidingDissipation(speciesDrumAndS1->getDissipation()*2./7.);
        speciesDrumAndS1->setSlidingFrictionCoefficient( ((slidingFrictionDrum + slidingFriction1)/2));
        speciesDrumAndS1->setSlidingStiffness(speciesDrumAndS1->getStiffness()*2.0/7.0);

        speciesDrumAndS1->setRollingStiffness(speciesDrumAndS1->getStiffness()*2.0/7.0);
        speciesDrumAndS1->setRollingFrictionCoefficient(((rollingFrictionDrum + rollingFriction1)/2));
        speciesDrumAndS1->setRollingDissipation(speciesDrumAndS1->getDissipation()*2./7.);

        speciesDrumAndS1->setTorsionStiffness(speciesDrumAndS1->getStiffness()*2.0/7.0);
        speciesDrumAndS1->setTorsionFrictionCoefficient(((torsionFrictionDrum + torsionFriction1)/2));
        speciesDrumAndS1->setTorsionDissipation(speciesDrumAndS1->getDissipation()*2./7.);
        //
        auto speciesDrumAndS2 = speciesHandler.getMixedObject(speciesDrum,speciesS2);

        speciesDrumAndS2->setCollisionTimeAndRestitutionCoefficient(tc, ((CORDrum + CORS2) / 2), massS1, massS2);

        speciesDrumAndS2->setSlidingDissipation(speciesDrumAndS2->getDissipation()*2./7.);
        speciesDrumAndS2->setSlidingFrictionCoefficient(((slidingFrictionDrum + slidingFriction2)/2));
        speciesDrumAndS2->setSlidingStiffness(speciesDrumAndS2->getStiffness()*2.0/7.0);

        speciesDrumAndS2->setRollingStiffness(speciesDrumAndS2->getStiffness()*2.0/7.0);
        speciesDrumAndS2->setRollingFrictionCoefficient(((rollingFrictionDrum + rollingFriction2)/2));
        speciesDrumAndS2->setRollingDissipation(speciesDrumAndS2->getDissipation()*2./7.);

        speciesDrumAndS2->setTorsionStiffness(speciesDrumAndS2->getStiffness()*2.0/7.0);
        speciesDrumAndS2->setTorsionFrictionCoefficient(((torsionFrictionDrum + torsionFriction2)/2));
        speciesDrumAndS2->setTorsionDissipation(speciesDrumAndS2->getDissipation()*2./7.);
        //
        auto speciesS1AndS2 = speciesHandler.getMixedObject(speciesS1,speciesS2);

        speciesS1AndS2->setCollisionTimeAndRestitutionCoefficient(tc, ((CORS1 + CORS2) / 2), massS1, massS2);

        speciesS1AndS2->setSlidingDissipation(speciesS1AndS2->getDissipation()*2./7.);
        speciesS1AndS2->setSlidingFrictionCoefficient(((rollingFriction1 + rollingFriction2)/2));
        speciesS1AndS2->setSlidingStiffness(speciesS1AndS2->getStiffness()*2.0/7.0);

        speciesS1AndS2->setRollingStiffness(speciesS1AndS2->getStiffness()*2.0/7.0);
        speciesS1AndS2->setRollingFrictionCoefficient(((rollingFriction1 + rollingFriction2)/2));
        speciesS1AndS2->setRollingDissipation(speciesS1AndS2->getDissipation()*2./7.);

        speciesS1AndS2->setTorsionStiffness(speciesS1AndS2->getStiffness()*2.0/7.0);
        speciesS1AndS2->setTorsionFrictionCoefficient(((torsionFriction1 + torsionFriction2)/2));
        speciesS1AndS2->setTorsionDissipation(speciesS1AndS2->getDissipation()*2./7.);

        //Add your walls below, and don't forget to set the species!
        Vec3D drumCenter = {0.0 , 0.0 , 0.0};
        wallHandler.clear();
        IntersectionOfWalls drumWall1;
        drumWall1.setAngularVelocity(Vec3D(0.0,RPSInitial*constants::pi*2.0 ,0.0));
        drumWall1.setSpecies(speciesDrum);
        drumWall1.addObject(Vec3D(1,0,0),  Vec3D(drumRadius, 0.0 ,0.0));
        wallHandler.copyAndAddObject(drumWall1);

        IntersectionOfWalls drumWall2;
        drumWall2.addObject(Vec3D(-1,0,0), Vec3D(-drumRadius,0.0,0.0));
        drumWall1.setAngularVelocity(Vec3D(0.0,RPSInitial*constants::pi*2.0 ,0.0));
        drumWall1.setSpecies(speciesDrum);
        wallHandler.copyAndAddObject(drumWall2);


        IntersectionOfWalls drumWall3;
        drumWall2.addObject(Vec3D(0,0,1), Vec3D(0.0,0.0,drumRadius));
        drumWall1.setAngularVelocity(Vec3D(0.0,RPSInitial*constants::pi*2.0 ,0.0));
        drumWall1.setSpecies(speciesDrum);
        wallHandler.copyAndAddObject(drumWall3);


        IntersectionOfWalls drumWall4;
        drumWall4.addObject(Vec3D(0,0,-1), Vec3D(0.0,0.0,-drumRadius));
        drumWall1.setAngularVelocity(Vec3D(0.0,RPSInitial*constants::pi*2.0 ,0.0));
        drumWall1.setSpecies(speciesDrum);
        wallHandler.copyAndAddObject(drumWall4);



        InfiniteWall w0;
        w0.setSpecies(speciesDrum);
        w0.set(Vec3D(0.,-1,0.),Vec3D(0,getYMin(),0));
        wallHandler.copyAndAddObject(w0);
        w0.set(Vec3D(0.,1,0.),Vec3D(0,getYMax(),0));
        wallHandler.copyAndAddObject(w0);

        //Now, either add particles or Insertion/Deletion boundaries

        SphericalParticle P0;
        double radius = 0.0;
        int numS1Inserted=0;
        int numS2Inserted=0;
        Vec3D pos;
        double x, z, y;
        int failCounter = 0;
        while( (numS1Inserted < numS1) || (numS2Inserted < numS2) )
        {
            int grn = random.getRandomNumber(1,numS1ToBeInserted+numS2ToBeInserted);

            if( grn > numS2ToBeInserted)
            {
                P0.setSpecies(speciesS1);
                P0.setRadius(radiusS1);

                failCounter = 0;
                do
                {
                    x = random.getRandomNumber(-drumRadius+radiusS1 , drumRadius-radiusS1);
                    y = random.getRandomNumber(getYMin()+radiusS1 , getYMax()-radiusS1);
                    z = random.getRandomNumber(-drumRadius+radiusS1 , drumRadius-radiusS1);
                    pos.X = x;
                    pos.Y = y;
                    pos.Z = z;

                    P0.setPosition(pos);
                    P0.setVelocity(Vec3D(0.0,0.0,0.0));

                    failCounter++;
                    if (failCounter==1000) break;

                } while (checkParticleForInteraction(P0));

                numS1ToBeInserted--;
                numS1Inserted++;
            }
            else
            {
                P0.setSpecies(speciesS2);
                P0.setRadius(radiusS2);

                failCounter = 0;
                do
                {
                    x = random.getRandomNumber(-drumRadius+radiusS2 , drumRadius-radiusS2);
                    y = random.getRandomNumber(getYMin()+radiusS2 , getYMax()-radiusS2);
                    z = random.getRandomNumber(-drumRadius+radiusS2 , drumRadius-radiusS2);
                    pos.X = x;
                    pos.Y = y;
                    pos.Z = z;

                    P0.setPosition(pos);
                    P0.setVelocity(Vec3D(0.0,0.0,0.0));

                    failCounter++;
                    if (failCounter==1000) break;
    
                } while (checkParticleForInteraction(P0));
    
    
                numS2ToBeInserted--;
                numS2Inserted++;
            }
            particleHandler.copyAndAddObject(P0);
    
            hGridRebuild();
        }
        logger(INFO, "Finished creating particles\n"
                     "Number of S1 particles inserted %\n"
                     "Number of S2 particles inserted %", numS1Inserted, numS2Inserted);
    
        //hGridRebuild();
    
        if ((numS1ToBeInserted == 0) && (numS2ToBeInserted == 0))
        {
            step = 2;
            logger(INFO, "\n \n \n"
                         "Particles settling down"
                         "--------------------------"
                         "\n\n\n");
            checkTime = getTime() + 1.0;
        }
    }
    void actionsBeforeTimeStep() override {

        if (step==2)
        {
            if (getTime() > checkTime)
            {
                logger(INFO, "Current KE %", getKineticEnergy());
                if (getKineticEnergy() < (10.0))
                {
                    step = 3;
                    double drumStartTime = getTime();
                    logger(INFO, "\n \n \n"
                                 "Start the drum rotation"
                                 "--------------------------"
                                 "\n\n\n");
                    // rotate the drum
                    wallHandler.getObject(0)->setAngularVelocity(
                            Vec3D(0.0, revolutionsPerSecond * constants::pi * 2.0, 0.0));
                    wallHandler.getObject(1)->setAngularVelocity(
                            Vec3D(0.0, revolutionsPerSecond * constants::pi * 2.0, 0.0));
                    wallHandler.getObject(2)->setAngularVelocity(
                            Vec3D(0.0, revolutionsPerSecond * constants::pi * 2.0, 0.0));
    
    
                }
                else
                {
                    checkTime = getTime() + 0.1;
                }
            }
        }
    }

    void setDrumRadius (double radius)
    {
        drumRadius = radius;
    }

    void setRevolutionSpeed (double rpm)
    {
        revolutionsPerSecond = rpm/60.0;
    }

    void setSizeAndDensityRatio (double sr, double dr)
    {
        sizeRatio = sr;
        densityRatio = dr;
    }

    void setFractionalPolydispersity(double fpd)
    {
        fractionalPolydispersity = fpd;
    }

    void setDrumFillFraction (double dff)
    {
        drumFillFraction = dff;
    }

    void setSpeciesVolumeFraction(double vf)
    {
        volumeFraction = vf;
    }

    void setFrictionCoeff(double pwf, double ppf)
    {
        particleWallFriction = pwf;
        particleParticleFriction = ppf;
    }

    void setCOR (double drumCOR, double COR1, double COR2)
    {
        CORDrum = drumCOR;
        CORS1 = COR1;
        CORS2 = COR2;
    }
    //a series of functions by which to easily set particles'
    //various frictional coefficients
    void setSlidingFriction (double drum, double f1, double f2)
    {
        slidingFrictionDrum = drum;
        slidingFriction1 = f1;
        slidingFriction2 = f2;
    }

    void setRollingFriction (double drum, double f1, double f2)
    {
        rollingFrictionDrum = drum;
        rollingFriction1 = f1;
        rollingFriction2 = f2;
    }

    void setTorsionFriction (double drum, double f1, double f2)
    {
        torsionFrictionDrum = drum;
        torsionFriction1 = f1;
        torsionFriction2 = f2;
    }

    double getDrumRadius()
    {
        return drumRadius;
    }

    //new functions to set the frequency and amplitude with which the drum is
    //vibrated
    void setVibrationFrequency (double f) {
        vibrationFreq = f;
    }
    void setVibrationAmplitude (double A) {
        vibrationAmp = A;
    }
private:

    double radiusS1,radiusS2;
    double rhoS1, rhoS2;
    double massS1, massS2;

    double CORDrum,CORS1,CORS2,tc;

    double sizeRatio;
    double densityRatio;
    double drumFillFraction;
    double volumeFraction;

    double particleWallFriction,particleParticleFriction;

    int numS1,numS2;
    int numS1ToBeInserted,numS2ToBeInserted;

    double drumRadius;
    double revolutionsPerSecond;
    double fractionalPolydispersity;

    double slidingFrictionDrum;
    double slidingFriction1;
    double slidingFriction2;

    double rollingFrictionDrum;
    double rollingFriction1;
    double rollingFriction2;

    double torsionFrictionDrum;
    double torsionFriction1;
    double torsionFriction2;

    int step;
    double checkTime;

    //New parameters to give the strength of applied vibrations
    double vibrationAmp;
    double vibrationFreq;
};

int main(int argc, char *argv[])
{
    RotatingDrum problem;
    /* Make sure to add an (unique) name here */
    problem.setName("SquareDrumTest");

    //setting locally-available variables to define key
    //parameters such that they can be automatically included in
    //the name of the file
    //*******************Excitation Properties*****************************************
    //Drum rotation rate (rpm)
    double rotRateBasal = 10.0;

    //vibration amplitude and frequency
    double fVib = 0.0;
    double aVib = 0.0;

    //*******************Frictional Properties*****************************************
    //sliding friction for species 1, 2 and wall
    double muSWall = 0.6;
    double muS1 = 0.16;
    double muS2 = 0.16;

    //rolling friction for species 1, 2 and wall
    double muRWall = 0.06;
    double muR1 = 0.016;
    double muR2 = 0.016;

    //torsion friction for species 1, 2 and wall
    double muTWall = 0.0;
    double muT1 = 0.000;
    double muT2 = 0.000;

    //the density ratio of particles
    double rhoRatio = 5000.0/2500.0;
    //the size ratio of particles
    double dRatio = 1;

    // computes rpm based on froude number and drumRad
    //the fraction of species 1 particles
    double specFrac = 0.5;

    double drumRad = 0.150/2; //diameter is 0.0150 m so drumRad is halve of the width
    double drumLength = 0.024; // depth of the drum is 0.024 m

    problem.setTimeMax(5.0);
    problem.setTimeStep(1.0/(2000.0 * 50.0));
    double froudeNumber = 0.22; //was 0.22
    double rotRate = 20; //pow(froudeNumber*9.81/drumRad,0.5);
    logger(INFO, "rotRate = %", rotRate);

    //*******************Setting Up Filename******************************************

    problem.setDrumRadius(drumRad);// in meters

    problem.setXMin(-(pow(2*pow(drumRad,2.0),0.5)));
    problem.setYMin(-drumLength/2);
    problem.setZMin(-(pow(2*pow(drumRad,2.0),0.5)));

    problem.setXMax((pow(2*pow(drumRad,2.0),0.5)));
    problem.setYMax(drumLength/2);// in meters
    problem.setZMax((pow(2*pow(drumRad,2.0),0.5)));

    problem.setGravity(Vec3D(0.,0.,-9.81));
    problem.setCOR(0.97,0.97,0.97);//drumWall, species1, species2 was all 0.97
    problem.setSizeAndDensityRatio(dRatio,rhoRatio);//size ratio and density ratio
    problem.setFractionalPolydispersity(0.00);
    problem.setDrumFillFraction(0.4);// At 0.5 the drum is 3/4 filled.
    problem.setSpeciesVolumeFraction(specFrac);//Species1 volume fraction

    problem.setSlidingFriction(muSWall,muS1,muS2); //wall, species1, species2
    problem.setRollingFriction(muRWall,muR1,muR2); //wall, species1, species2
    problem.setTorsionFriction(muTWall,muT1,muT2); //wall, species1, species2

    problem.setRevolutionSpeed(rotRate);//rpm

    problem.setSaveCount(2000);
    problem.readArguments(argc,argv);

    problem.dataFile.setFileType(FileType::ONE_FILE);
    problem.restartFile.setFileType(FileType::ONE_FILE);
    problem.fStatFile.setFileType(FileType::ONE_FILE);
    problem.eneFile.setFileType(FileType::ONE_FILE);
    problem.wallHandler.setWriteVTK(true);
    problem.setParticlesWriteVTK(true);
    problem.solve();
    return 0;
}
