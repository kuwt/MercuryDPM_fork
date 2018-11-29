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

class BoosterWithWedges : public Booster {
public:

    BoosterWithWedges(Mdouble wedgeHeight, unsigned int wedgeNumber)
        : Booster(), wedgeHeight_(wedgeHeight), wedgeNumber_(wedgeNumber)
    {
        setName("BoosterWithWedges");
    }

    void addBlades()
    {
        //wedge
        wedgesDown.reserve(wedgeNumber_);
        wedgesRight.reserve(wedgeNumber_);
        wedgesUp.reserve(wedgeNumber_);
        wedgesLeft.reserve(wedgeNumber_);

        IntersectionOfWalls defaultWedge;
        defaultWedge.setSpecies(wallSpecies);

        Mdouble lidWidth = 216.73e-3;
        Mdouble upperWidth = 976.57e-3;
        Mdouble centreWidth = 2000.99e-3;

        double x = lidWidth+upperWidth;
        double dx = 0.5*centreWidth/static_cast<double>(wedgeNumber_-0.5);
        for (unsigned int i=0; i<wedgeNumber_; ++i)
        {
            //double x = (static_cast<double>(i)+0.5)*getXMax()/ static_cast<double>(wedgeNumber_);
            wedgesDown.push_back(wallHandler.copyAndAddObject(defaultWedge));
            wedgesDown.back()->setPosition(Vec3D(x,0.0,0.0));
            wedgesUp.push_back(wallHandler.copyAndAddObject(defaultWedge));
            wedgesUp.back()->setPosition(Vec3D(x,0.0,0.0));
            x += dx;
            wedgesRight.push_back(wallHandler.copyAndAddObject(defaultWedge));
            wedgesRight.back()->setPosition(Vec3D(x,0.0,0.0));
            wedgesLeft.push_back(wallHandler.copyAndAddObject(defaultWedge));
            wedgesLeft.back()->setPosition(Vec3D(x,0.0,0.0));
            x += dx;
        }

        actionsAfterTimeStep();
    }

    virtual void actionsAfterTimeStep()
    {
        Mdouble angle = revolutionsPerSecond * getTime() * 2.0 * constants::pi;
        for (IntersectionOfWalls* wedge : wedgesDown)
            setWedgeAngle(wedge,angle);
        for (IntersectionOfWalls* wedge : wedgesRight)
            setWedgeAngle(wedge,angle+constants::pi/2.0);
        for (IntersectionOfWalls* wedge : wedgesUp)
            setWedgeAngle(wedge,angle+2.0*constants::pi/2.0);
        for (IntersectionOfWalls* wedge : wedgesLeft)
            setWedgeAngle(wedge,angle+3.0*constants::pi/2.0);
    }

    virtual void setWedgeAngle(IntersectionOfWalls* wedge, Mdouble angle)
    {
        wedge->clear();
        Mdouble s = std::sin(angle), c = std::cos(angle);
        Mdouble s2 = std::sqrt(2), s3=std::sqrt(3);
        //Corners of the wedge
        static Vec3D D = Vec3D(0,0,s2/s3);
        static Vec3D A = Vec3D(-0.5,-0.5/s3,0.0);
        static Vec3D B = Vec3D(0.5,-0.5/s3,0.0);
        static Vec3D C = Vec3D(0,1/s3,0.0);
        static Vec3D DA = A-D;
        static Vec3D DB = B-D;
        static Vec3D DC = C-D;
        static Vec3D normalA = Vec3D::cross(DC,DB);
        static Vec3D normalB = Vec3D::cross(DA,DC);
        static Vec3D normalC = Vec3D::cross(DB,DA);
        Vec3D DRotated = Vec3D(D.X, c*D.Y + s*D.Z,  -s*D.Y+c*D.Z);
        Vec3D normalARotated = Vec3D(normalA.X, c*normalA.Y - s*normalA.Z,  s*normalA.Y+c*normalA.Z);
        Vec3D normalBRotated = Vec3D(normalB.X, c*normalB.Y - s*normalB.Z,  s*normalB.Y+c*normalB.Z);
        Vec3D normalCRotated = Vec3D(normalC.X, c*normalC.Y - s*normalC.Z,  s*normalC.Y+c*normalC.Z);
        Vec3D anchorPosition = Vec3D(wedge->getPosition().X, (getZMax()-wedgeHeight_)*s, -(getZMax()-wedgeHeight_)*c);
        wedge->addObject(normalARotated, anchorPosition);
        wedge->addObject(normalBRotated, anchorPosition);
        wedge->addObject(normalCRotated, anchorPosition);
        std::cout << "anchorPosition" << anchorPosition << std::endl;
        std::cout << normalA << std::endl;
        std::cout << normalARotated << std::endl;
        std::cout << normalB << std::endl;
        std::cout << normalBRotated << std::endl;
    }

    void write(std::ostream& os, bool x) const
    {
        Mercury3D::write(os, x);
        os  << " wedgeNumber " << wedgeNumber_ << std::endl
            << " wedgeHeight " << wedgeHeight_ << std::endl;
    }
protected:
    // pointers to the wedges, split into groups denoting the four initial positions of the anchor points
    std::vector<IntersectionOfWalls*> wedgesDown, wedgesRight, wedgesUp, wedgesLeft;
public:
    Mdouble wedgeHeight_;
    unsigned int wedgeNumber_;
};

int main(int argc UNUSED, char *argv[] UNUSED)
{
    Mdouble wedgeHeight = 0.5*64.25e-2;
    unsigned int wedgeNumber = 2;
    BoosterWithWedges booster(wedgeHeight, wedgeNumber);
    
    //booster.autoNumber();
    //std::vector<int> studyNum=booster.get2DParametersFromRunNumber(3,3);
    //booster.wedgeInclination_ = (10.0 + 10.0*(studyNum[1]-1.))*constants::pi/180.0;
    //booster.wedgeHeight_ = (0.2 + 0.2*(studyNum[2]-1.))*64.25e-2;
    
    booster.drumInclination = 5.0*constants::pi/180.0;
    booster.revolutionsPerSecond = 8.0/60.0;
    booster.particleDiameter = 7.2e-2/1.2599; //cuberoot(2)=1.2599
    booster.particleSpecies->setDensity(2000.0);
    booster.collisionTime = 0.01; //soft particles, but relative overlap due to gravity, o/d ~ (tc/tg)^2 = 0.0136, is still small
    booster.restitutionCoefficient = 0.001; // very dissipative particles
    booster.setXBallsAdditionalArguments(" -v0 -solidf  -w 1200 -s 0.76 -noborder 3");
    booster.setTimeMax(5.0*60.0/8.0); //run 5 revolutions
    //booster.setTimeMax(1e-20); //run 5 revolutions
    booster.solve();
    return 0;
}
