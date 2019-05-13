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

//based on /storage2/usr/people/sluding/MDCC/C3DshearXL30/MU0_LONG2
#include <Species/LinearViscoelasticFrictionSpecies.h>
#include <Particles/BaseParticle.h>
#include "Mercury3D.h"
#include "Boundaries/AngledPeriodicBoundary.h"
#include "Walls/AxisymmetricIntersectionOfWalls.h"

class ShearCell3DSelfTest : public Mercury3D {
public:
    ShearCell3DSelfTest (Mdouble outerRadius, Mdouble splitRadius, Mdouble innerRadius, Mdouble height, Mdouble angularVelocity)
        : outerRadius_(outerRadius), splitRadius_(splitRadius), innerRadius_(innerRadius), height_(height), angularVelocity_(angularVelocity)
    {
        //check input variables
        if (innerRadius_ <0) {
            std::cerr << "error in constructor of ShearCell3DSelfTest; InnerRadius has to be positive" << std::endl;
            exit(-1);
        } else if (innerRadius_ > splitRadius_) {
            std::cerr << "error in constructor of ShearCell3DSelfTest; InnerRadius has to be smaller than SplitRadius" << std::endl;
            exit(-1);
        } else if (splitRadius_ > outerRadius_) {
            std::cerr << "error in constructor of ShearCell3DSelfTest; SplitRadius has to be smaller than OuterRadius" << std::endl;
            exit(-1);
        } else if (height_ <0) {
            std::cerr << "error in constructor of ShearCell3DSelfTest; DomainHeight has to be positive" << std::endl;
            exit(-1);
        }

        //set default name
        setName("ShearCell3DSelfTest");
        
        //set gravity (extra strong since we have few particle layers)
        setGravity(Vec3D(0.0,0.0,-10.0));

        //use particles of unit mass and diameter 
        auto particleSpecies =speciesHandler.copyAndAddObject(LinearViscoelasticFrictionSpecies());
        particleSpecies->setDensity(6.0/constants::pi);

        //set inter-particle contact properties
        particleSpecies->setStiffness(1e3); //tc \sim pi/sqrt(k/m) = 1e-1
        particleSpecies->setDissipation(25); //dissipation_<<sqrt(4km)=50
        particleSpecies->setSlidingStiffness(0.4* particleSpecies->getStiffness());
        particleSpecies->setSlidingDissipation(0.4* particleSpecies->getDissipation());
        particleSpecies->setSlidingFrictionCoefficient(0.1);
        
        // wall-particle contacts
        auto wallSpecies =speciesHandler.copyAndAddObject(particleSpecies);
        wallSpecies->setSlidingStiffness(0);
        wallSpecies->setSlidingDissipation(0);
        wallSpecies->setSlidingFrictionCoefficient(0);
        wallSpecies->setSlidingFrictionCoefficientStatic(0);
        auto particleWallSpecies =speciesHandler.getMixedObject(particleSpecies, wallSpecies);
        particleWallSpecies->setStiffness(particleSpecies->getStiffness()); //tc \sim pi/sqrt(k/m) = 1e-1
        particleWallSpecies->setDissipation(particleSpecies->getDissipation()); //dissipation_<<sqrt(4km)=50
        particleWallSpecies->setSlidingStiffness(0.4* particleSpecies->getStiffness());
        particleWallSpecies->setSlidingDissipation(0.4* particleSpecies->getDissipation());
        particleWallSpecies->setSlidingFrictionCoefficient(1.0);
        particleWallSpecies->setRollingStiffness(0.4* particleSpecies->getStiffness());
        particleWallSpecies->setRollingDissipation(0.4* particleSpecies->getDissipation());
        particleWallSpecies->setRollingFrictionCoefficient(1.0);

        //set time step
        setTimeMax(10.0);
        setTimeStep(5e-3);
        setSaveCount(200);

        //set domain
        setXMin(0.0);
        setYMin(0.0);
        setZMin(0.0);
        setXMax(outerRadius_);
        setYMax(outerRadius_);
        setZMax(height_);

        //set boundaries
        //quarter cell
        AngledPeriodicBoundary b;
        Vec3D normal_left(1.0,0.0,0.0);
        Vec3D normal_right(0.0,-1.0,0.0);
        Vec3D origin(0.0,0.0,0.0);
        b.set(normal_left,normal_right,origin);
        boundaryHandler.copyAndAddObject(b);

        //set walls
        AxisymmetricIntersectionOfWalls w;
        Vec3D AxisDirection(0.,0.,1.);
        Vec3D PointOnAxis(0.,0.,0.);
        Vec3D PrismAxis(0.,-1.,0.);
        w.setPosition(PointOnAxis);
        w.setOrientation(AxisDirection);
        w.setSpecies(wallSpecies);

        //add points in anti-clockwise direction around the prism axis
        std::vector<Vec3D> Points(2);
        Points[0] = Vec3D(innerRadius_,0.0,getZMin());
        Points[1] = Vec3D(innerRadius_,0.0,getZMax());
        w.createOpenPrism(Points,PrismAxis);
        wallHandler.copyAndAddObject(w);

        Points.resize(3);
        Points[0] = Vec3D(splitRadius_,0.0,getZMin()-1.0);
        Points[1] = Vec3D(splitRadius_,0.0,getZMin());
        Points[2] = Vec3D(innerRadius_,0.0,getZMin());
        w.createOpenPrism(Points,PrismAxis);
        wallHandler.copyAndAddObject(w);

        w.setAngularVelocity(Vec3D(0.0, 0.0, angularVelocity_));

        Points.resize(2);
        Points[0] = Vec3D(outerRadius_,0.0,getZMax());
        Points[1] = Vec3D(outerRadius_,0.0,getZMin());
        w.createOpenPrism(Points,PrismAxis);
        wallHandler.copyAndAddObject(w);

        Points.resize(3);
        Points[0] = Vec3D(outerRadius_,0.0,getZMin());
        Points[1] = Vec3D(splitRadius_,0.0,getZMin());
        Points[2] = Vec3D(splitRadius_,0.0,getZMin()-1.0);
        w.createOpenPrism(Points,PrismAxis);
        wallHandler.copyAndAddObject(w);

    }

private:
    /////////////////////////////////////////////////////////////////////////////////////////////////////
    ///This is were the walls are implemented
    /////////////////////////////////////////////////////////////////////////////////////////////////////
    void setupInitialConditions() 
    {
        SphericalParticle p;
        p.setSpecies(speciesHandler.getObject(0));
        p.setVelocity(Vec3D(0.0,0.0,0.0));
        p.setAngularVelocity(Vec3D(0.0,0.0,0.0));
        p.setRadius(0.5);

        Mdouble dRadius = 2.0* p.getRadius(); //minimum distance between particles
        Mdouble nRadial = std::floor((outerRadius_ - innerRadius_)/dRadius)-1;
        Mdouble dRadiusRadial = (outerRadius_ - innerRadius_)/ (nRadial+1); //actual radial distance between particles
        for (Mdouble z=0.5*dRadius; z< height_; z+=dRadius) {
            for (int i=1; i<= nRadial; i++) {
                Mdouble Radius = outerRadius_ -i*dRadiusRadial;
                Mdouble circumference = Radius * 0.5 * constants::pi;
                Mdouble n = std::floor(circumference / dRadius);
                Mdouble dAlpha = 0.5 * constants::pi / (double) n; //actual angular distance between particles
                for (Mdouble alpha=dAlpha/4.0+dAlpha/2.0*fmod(i,2); alpha<0.5*constants::pi; alpha+=dAlpha) {
                    p.setPosition(Vec3D(Radius*sin(alpha+z),Radius*cos(alpha+z),z));
                    particleHandler.copyAndAddObject(p);
                }
            }
        }
        setHGridMaxLevels(1);
    }

    void actionsBeforeTimeStep ()
    {
        //this is because the angular velocity of the wall does not calculate the orientation correctly
        wallHandler.getObject(2)->setOrientation(Vec3D(0.0, 0.0, 1));
        wallHandler.getObject(3)->setOrientation(Vec3D(0.0, 0.0, 1));
    }

    Mdouble outerRadius_;
    Mdouble splitRadius_;
    Mdouble innerRadius_;
    Mdouble height_;
    Mdouble angularVelocity_;
};


int main() {
    ShearCell3DSelfTest SC(10.0,7.75,5.0,3.0,0.01);
    SC.setXBallsAdditionalArguments("-v0 -solidf -3dturn 1");
    //comment the next line to get output
    SC.setSaveCount(1e8);
    SC.solve();
    return 0;
}
