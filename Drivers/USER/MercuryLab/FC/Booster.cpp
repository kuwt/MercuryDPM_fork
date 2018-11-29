//Copyright (c) 2013-2017, The MercuryDPM Developers Team. All rights reserved.
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

#include "Booster.h"

class BoosterWithStraightBlade : public Booster
{
public:

    BoosterWithStraightBlade(Mdouble bladeThickness, Mdouble bladeHeight)
    : Booster(), bladeThickness_(bladeThickness), bladeHeight_(bladeHeight)
    {   }

    void addBlades()
    {
        //blade
        std::cout << "creating blade0" << std::endl;
        blade0 = wallHandler.copyAndAddObject(IntersectionOfWalls());
        blade1 = wallHandler.copyAndAddObject(IntersectionOfWalls());
        blade2 = wallHandler.copyAndAddObject(IntersectionOfWalls());
        blade3 = wallHandler.copyAndAddObject(IntersectionOfWalls());
        blade0->setSpecies(wallSpecies);
        blade1->setSpecies(wallSpecies);
        blade2->setSpecies(wallSpecies);
        blade3->setSpecies(wallSpecies);
        setBladeAngle(blade0,0.0);
//        blade0->addObject(Vec3D(1,0,0), 1);
//        blade0->addObject(Vec3D(-1,0,0), 1);
//        blade0->addObject(Vec3D(0,1,0), 1);
//        blade0->addObject(Vec3D(0,-1,0), 1);
//        blade0->addObject(Vec3D(0,0,1), 1);
//        blade0->addObject(Vec3D(0,0,-1), 1);
        setBladeAngle(blade1,constants::pi/2.0);
        setBladeAngle(blade2,2.0*constants::pi/2.0);
        setBladeAngle(blade3,3.0*constants::pi/2.0);
    }

    void setBladeAngle(IntersectionOfWalls* blade, Mdouble angle)
    {
        blade->clear();
        Vec3D longAxis = Vec3D(0.0, std::sin(angle), -std::cos(angle));
        Vec3D shortAxis = Vec3D(0.0, -longAxis.Z, longAxis.Y);
        blade->addObject(  longAxis, (getZMax()-bladeHeight_)*longAxis);
        blade->addObject( shortAxis, -0.5*bladeThickness_*shortAxis);
        blade->addObject(-shortAxis,  0.5*bladeThickness_*shortAxis);
    }

    void actionsAfterTimeStep()
    {
        Mdouble angle = revolutionsPerSecond * getTime() * 2.0 * constants::pi;
        setBladeAngle(blade0,angle);
        setBladeAngle(blade1,angle+constants::pi/2.0);
        setBladeAngle(blade2,angle+2.0*constants::pi/2.0);
        setBladeAngle(blade3,angle+3.0*constants::pi/2.0);
    }

protected:
    IntersectionOfWalls* blade0, * blade1, * blade2, * blade3;
    Mdouble bladeThickness_;
    Mdouble bladeHeight_;
};

int main(int argc UNUSED, char *argv[] UNUSED)
{
    Mdouble bladeInclination = 30.0*constants::pi/180.0; // default angle of 5.71059314 deg such that tan(angle)=0.1 (it's easy to check the normals then)
    Mdouble bladeThickness = 7.2e-2;
    Mdouble bladeWidth = 45.0/180.0*constants::pi*64.25e-2; //set in degrees
    Mdouble bladeHeight = 0.5*64.25e-2;
    unsigned int bladeNumber = 3;
    BoosterWithStraightBlade booster(bladeThickness,bladeHeight);

    booster.drumInclination = 5.0*constants::pi/180.0;
    booster.revolutionsPerSecond = 8.0/60.0;
    booster.particleDiameter = 7.2e-2/1.2599; //cuberoot(2)=1.2599
    booster.particleSpecies->setDensity(2000.0);
    booster.collisionTime = 0.01; //soft particles, but relative overlap due to gravity, o/d ~ (tc/tg)^2 = 0.0136, is still small
    booster.restitutionCoefficient = 0.001; // very dissipative particles
    booster.setXBallsAdditionalArguments(" -v0 -solidf  -w 1200 -s 0.76 -noborder 3");
    booster.setTimeMax(5.0*60.0/8.0); //run 5 revolutions
    //booster.setTimeMax(0.1); //for testing
    booster.autoNumber();
     // Remove this line to go back to old style simulations
    booster.setFillingFraction(0.5);
    booster.makePolydispersed();
    booster.readArguments(argc, argv);
    booster.solve();
   
   
    return true;
}
