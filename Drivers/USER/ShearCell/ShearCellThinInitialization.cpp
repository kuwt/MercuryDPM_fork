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
#include <Species/LinearPlasticViscoelasticSlidingFrictionSpecies.h>
#include "Chute.h"
#include "Boundaries/AngledPeriodicBoundary.h"
#include "Walls/AxisymmetricIntersectionOfWalls.h"

class ShearCellMDCLR : public Chute {
public:
    // void computeInternalForces(BaseParticle *PI, BaseParticle *PJ) {
    //     computePlasticInternalForces(PI, PJ);
    // }

    ShearCellMDCLR ()
    {
        //set angular velocity of the outer walls
        AngularVelocity = 0.01;
        
        //set Species
        auto species =speciesHandler.copyAndAddObject(LinearPlasticViscoelasticSlidingFrictionSpecies());
        species->setDensity(2000.0);
        species->setLoadingStiffness(1e2);
        //setUnloadingStiffnessMax(1.1e2);
        species->setCohesionStiffness(0.0);
        species->setSlidingStiffness(1.2e1);
        species->setSlidingFrictionCoefficient(0.01);
        species->setDissipation(2e-3);
        species->setSlidingDissipation(0.5e-3);
        
        //set time step
        setTimeMax(100.0); //150 is the time of the loaded data file
        setTimeStep(5e-6);
        setSaveCount(2e4);
        eneFile.setSaveCount(2e3);
        
        //set gravity
        setGravity(Vec3D(0.0,0.0,-9.81));
        
        //load particles (setting mu!=0 loads particles as TSP)
		//particleHandler.setStorageCapacity(37012);
        //readDataFile(filename.c_str());
        
        //set default name
        setName("ShearCellMDCLR");
        
        //begin: set geometry
        
        //set domain size
        InnerRadius=0.01466667;
        Mdouble SplitRadius=0.08500000;
        OuterRadius=0.11000000;
        Mdouble DomainHeight=0.06891868232999998;
        setXMin(0.0);
        setYMin(0.0);
        setZMin(0.0);
        setXMax(OuterRadius);
        setYMax(OuterRadius);
        setZMax(DomainHeight);
        
        //set boundaries
        //quarter cell
        AngledPeriodicBoundary B0;
        Vec3D normal_left(1.0,0.0,0.0);
        Mdouble th = 90.0/7.0*constants::pi/180.0;
        Vec3D normal_right(cos(th),-sin(th),0.0);
        Vec3D origin(0.0,0.0,0.0);
        B0.set(normal_left,normal_right,origin);
        boundaryHandler.copyAndAddObject(B0);
        
        //set axisymmetric walls
        AxisymmetricIntersectionOfWalls w1;
        Vec3D AxisDirection(0.,0.,1.);
        Vec3D PointOnAxis(0.,0.,0.);
        Vec3D PrismAxis(0.,-1.,0.);
        w1.setPosition(PointOnAxis);
        w1.setOrientation(AxisDirection);
        //w1.setSpecies(speciesHandler.getObject(0));
        
        //add points in anti-clockwise direction around the prism axis
        //inner side wall
        std::vector<Vec3D> Points(2);
        Points[0] = Vec3D(InnerRadius,0.0,getZMin());
        Points[1] = Vec3D(InnerRadius,0.0,getZMax());
        w1.createOpenPrism(Points,PrismAxis);
        wallHandler.copyAndAddObject(w1);
        
        //inner base wall
        Points.resize(3);
        Points[0] = Vec3D(SplitRadius,0.0,getZMin()-1.0);
        Points[1] = Vec3D(SplitRadius,0.0,getZMin());
        Points[2] = Vec3D(InnerRadius,0.0,getZMin());
        w1.createOpenPrism(Points,PrismAxis);
        wallHandler.copyAndAddObject(w1);
        
        //give outer walls velocity
        w1.setVelocity(Vec3D(0.0,AngularVelocity,0.0));
        
        //outer side wall
        Points.resize(2);
        Points[0] = Vec3D(OuterRadius,0.0,getZMax());
        Points[1] = Vec3D(OuterRadius,0.0,getZMin());
        w1.createOpenPrism(Points,PrismAxis);
        wallHandler.copyAndAddObject(w1);
        
        //outer base wall
        Points.resize(3);
        Points[0] = Vec3D(OuterRadius,0.0,getZMin());
        Points[1] = Vec3D(SplitRadius,0.0,getZMin());
        Points[2] = Vec3D(SplitRadius,0.0,getZMin()-1.0);
        w1.createOpenPrism(Points,PrismAxis);
        wallHandler.copyAndAddObject(w1);
    }
    
private:
    /////////////////////////////////////////////////////////////////////////////////////////////////////
    ///This is were the walls are implemented
    /////////////////////////////////////////////////////////////////////////////////////////////////////
    void setupInitialConditions()
    {
    }
    
    void actionsBeforeTimeStep () {
        static int counter=0;
        std::cout << ++counter << std::endl;
    }
    
    void actionsAfterTimeStep () {
    }
    
    Mdouble OuterRadius;
    Mdouble InnerRadius;
    Mdouble AngularVelocity;
    std::vector<unsigned int> MPIndex;
    std::vector<Mdouble> MPRadius;
    std::vector<Mdouble> MPAngle;
    std::vector<Mdouble> MPHeight;
};


int main() {
    
    ShearCellMDCLR SC;
    SC.setXBallsAdditionalArguments("-v0 -solidf");
    SC.setHGridMaxLevels(2);
    SC.setTimeMax(50.5* SC.getTimeStep());
    SC.dataFile.setFileType(FileType::NO_FILE);
    SC.fStatFile.setFileType(FileType::NO_FILE);
    SC.restartFile.setFileType(FileType::MULTIPLE_FILES);
    SC.solve();
    return 0;
}
