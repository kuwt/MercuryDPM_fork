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

#include <Species/LinearViscoelasticSlidingFrictionSpecies.h>
#include <Species/LinearViscoelasticSpecies.h>
#include "DPMBase.h"
#include "Boundaries/AngledPeriodicBoundary.h"
class AngledPeriodicBoundaryTestA : public DPMBase {
public:
    /////////////////////////////////////////////////////////////////////////////////////////////////////
    ///This is were the walls are implemented
    /////////////////////////////////////////////////////////////////////////////////////////////////////
    void setupInitialConditions() override
    {
        setXMax(10.0);
        setYMax(10.0);

        auto S=speciesHandler.copyAndAddObject(LinearViscoelasticSpecies());
        S->setDissipation(1);
            
        AngledPeriodicBoundary B0;
        Vec3D normal_left(1.0,0.0,0.0);
        Vec3D normal_right(0.0,-1.0,0.0);
        Vec3D origin(1.0,1.0,0.0);
        B0.set(normal_left,normal_right,origin);
        boundaryHandler.copyAndAddObject(B0);

        SphericalParticle P0;
        P0.setSpecies(speciesHandler.getObject(0));
        P0.setPosition(Vec3D(9.0,0.5,0.0));      
        P0.setVelocity(Vec3D(0.0,1.0,0.0));
        P0.setOrientation(Vec3D(0.0,0.0,1.0));
        P0.setAngularVelocity(Vec3D(0.0,0.0,1.0));
        P0.setRadius(0.5);
        B0.distance(P0);
        //particleHandler.copyAndAddObject(P0);
        B0.shiftPosition(&P0); //test shifting
        particleHandler.copyAndAddObject(P0);   
        B0.shiftPosition(&P0); //test shifting
        particleHandler.copyAndAddObject(P0);    

        writeRestartFile();
    }

};


class AngledPeriodicBoundaryTestB : public DPMBase {
    /////////////////////////////////////////////////////////////////////////////////////////////////////
    ///This is were the walls are implemented
    /////////////////////////////////////////////////////////////////////////////////////////////////////
    void setupInitialConditions() override
    {
        setXMax(10.0);
        setYMax(10.0);

        auto S=speciesHandler.copyAndAddObject(LinearViscoelasticSlidingFrictionSpecies());
        S->setDissipation(1.0);
        S->setSlidingFrictionCoefficient(1.0);
        S->setSlidingStiffness(S->getStiffness());
        S->setSlidingDissipation(S->getDissipation());

        AngledPeriodicBoundary B0;
        Vec3D normal_left(1.0,-0.5,0.0);
        Vec3D normal_right(0.0,-1.0,0.0);
        Vec3D origin(1.0,1.0,0.0);
        B0.set(normal_left,normal_right,origin);
        boundaryHandler.copyAndAddObject(B0);

        SphericalParticle P0;
        P0.setSpecies(speciesHandler.getObject(0));
        P0.setPosition(Vec3D(9.0,0.5+1e-8,0.0));      
        P0.setVelocity(Vec3D(0.0,1.0,0.0));
        P0.setOrientation(Vec3D(0.0,0.0,1.0));
        P0.setAngularVelocity(Vec3D(0.0,0.0,1.0));
        P0.setRadius(0.5);
        B0.distance(P0);
        B0.shiftPosition(&P0); //test shifting
        B0.shiftPosition(&P0); //test shifting
        particleHandler.copyAndAddObject(P0);    

        P0.setPosition(Vec3D(9.0,1.5,0.0));      
        particleHandler.copyAndAddObject(P0);    

        //a standard contact in comparison
        //P0.setPosition(Vec3D(7.0,5.0,0.0));      
        //particleHandler.copyAndAddObject(P0);    
        //P0.setPosition(Vec3D(7.0,4.0+1e-8,0.0));      
        //particleHandler.copyAndAddObject(P0);    

    }

public:
    void output() {
        //this requires the user to non-privatize the function below
        //checkAndDuplicatePeriodicParticles();
        writeRestartFile();
    }

};

class AngledPeriodicBoundaryTestC : public DPMBase {
    /////////////////////////////////////////////////////////////////////////////////////////////////////
    ///This is were the walls are implemented
    /////////////////////////////////////////////////////////////////////////////////////////////////////
    void setupInitialConditions() override
    {
        setXMax(10.0);
        setYMax(10.0);

        setGravity(Vec3D(0.0,0.0,0.0));
        auto S=speciesHandler.copyAndAddObject(LinearViscoelasticSlidingFrictionSpecies());
        S->setSlidingFrictionCoefficient(0.5);
        S->setStiffness(50000);
        S->setDissipation(0);
        S->setSlidingStiffness(0.4* S->getStiffness());
        S->setSlidingDissipation(0.4* S->getDissipation());
        S->setDensity(6.0/constants::pi);

        AngledPeriodicBoundary B0;
        Vec3D normal_left(1.0,-0.5,0.0);
        Vec3D normal_right(0.0,-1.0,0.0);
        Vec3D origin(0.0,0.0,0.0);
        B0.set(normal_left,normal_right,origin);
        boundaryHandler.copyAndAddObject(B0);

        SphericalParticle P0;
        P0.setSpecies(speciesHandler.getObject(0));
        P0.setOrientation(Vec3D(0.0,0.0,1.0));
        P0.setAngularVelocity(Vec3D(0.0,0.0,2.0));

        //collision before one of the particles passes the periodic wall
        P0.setPosition(Vec3D(9.0,1.5,0.0));      
        P0.setVelocity(Vec3D(0.3,-2.0,0.0));
        P0.setRadius(0.5);
        particleHandler.copyAndAddObject(P0);    
        P0.setPosition(Vec3D(9.0,0.25,0.0));      
        P0.setVelocity(Vec3D(0.0,0.0,0.0));
        P0.setRadius(0.5);
        particleHandler.copyAndAddObject(P0);    

        //collision with particles on different sides of the periodic wall
        P0.setPosition(Vec3D(7.0,0.75,0.0));      
        P0.setVelocity(Vec3D(0.3,-2.0,0.0));
        P0.setRadius(0.5);
        particleHandler.copyAndAddObject(P0);    
        P0.setPosition(Vec3D(7.0,-0.5,0.0));      
        P0.setVelocity(Vec3D(0.0,0.0,0.0));
        P0.setRadius(0.5);
        particleHandler.copyAndAddObject(P0);    

        //collision away from periodic wall
        P0.setPosition(Vec3D(5.0,-0.0,0.0));      
        P0.setVelocity(Vec3D(0.3,-2.0,0.0));
        P0.setRadius(0.5);
        particleHandler.copyAndAddObject(P0);    
        P0.setPosition(Vec3D(5.0,-1.25,0.0));      
        P0.setVelocity(Vec3D(0.0,0.0,0.0));
        P0.setRadius(0.5);
        particleHandler.copyAndAddObject(P0);    

        //a standard contact in comparison
        //P0.setPosition(Vec3D(7.0,5.0,0.0));      
        //particleHandler.copyAndAddObject(P0);    
        //P0.setPosition(Vec3D(7.0,4.0+1e-8,0.0));      
        //particleHandler.copyAndAddObject(P0);    

    }
};


int main() {
    AngledPeriodicBoundaryTestA APBa;
    APBa.setSaveCount(1);
    APBa.setName("AngledPeriodicBoundaryTestA");
    APBa.setFileType(FileType::ONE_FILE);
    APBa.setTime(0);
    APBa.setTimeStep(1e-20);
    APBa.setTimeMax(0.5e-20);
    APBa.setupInitialConditions();

    AngledPeriodicBoundaryTestB APBb;
    APBb.setSaveCount(1);
    APBb.setName("AngledPeriodicBoundaryTestB");
    APBa.setFileType(FileType::ONE_FILE);
    APBb.setTime(0);
    APBb.setTimeStep(1e-20);
    APBb.setTimeMax(0.5e-20);
    APBb.solve();
    APBb.output();

    AngledPeriodicBoundaryTestC APBc;
    APBc.setSaveCount(30);
    APBc.setName("AngledPeriodicBoundaryTestC");
    APBa.setFileType(FileType::ONE_FILE);
    APBc.setTime(0);
    APBc.setTimeStep(1e-4);
    APBc.setTimeMax(1.0);
    APBc.solve();
}
