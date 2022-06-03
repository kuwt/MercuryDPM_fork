//
// Created by irana on 9/13/17.
//


#include <sstream>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <Walls/InfiniteWall.h>
#include <Boundaries/CubeInsertionBoundary.h>
#include "Mercury3D.h"
#include "Species/LinearViscoelasticSlidingFrictionSpecies.h"
#include "Hopper.h"
#include "../BidispersedChute/BidispersedChute.h"


class HopperForAvalanche : public BidispersedChute
{
public:
    HopperForAvalanche(BidispersedChuteParameters parameters) : BidispersedChute(parameters)
    {
        setName("HopperIrana");
        setTimeStep(1e-3);
        setTimeMax(100);
        setMin(0,0,-10);
        setMax(100,10,100);
        setSaveCount(500);
        setRoughBottomType(RoughBottomType::MULTILAYER);
        setChuteLength(100);
        makeLong();
    }
    
    void setupInitialConditions() override
    {
        BidispersedChute::setupInitialConditions();
    
        setMin(-100,0,-10);
        Hopper hopper;
        hopper.setHopperLength(40);
        hopper.setHopperExitLength(10);
        hopper.setHopperHeight(40);
        hopper.setHopperAngle(45./180 * constants::pi);
        hopper.setHopperLift(-7);
        hopper.setHopperLowestPoint(5);
        hopper.setHopperShift(-13);
        hopper.setHopperDimension(1);
        hopper.makeHopper(wallHandler);
        setWallsWriteVTK(FileType::ONE_FILE);
        writeVTKFiles();
        setWallsWriteVTK(FileType::NO_FILE);
        InfiniteWall w0;
        w0.set({-1, 0, 0}, {-25, 0, 0});
    
        for (BaseWall* wall : wallHandler)
        {
            wall->setSpecies(speciesHandler.getObject(0));
        }
    
        SphericalParticle p0;
        p0.setSpecies(speciesHandler.getObject(0));
        p0.setRadius(0.5);
        CubeInsertionBoundary insertionBoundary;
        insertionBoundary.set(&p0, 50, {-20, 0, 30}, {5, 10, 35}, {0, 0, 0}, {0, 0, 0});
        boundaryHandler.copyAndAddObject(insertionBoundary);
    
        setMin(0,0,-10);
        setMax(100,10,100);
    }
};

// A quasi-2D inclined plane with hopper inflow conditions,
// and deletion of particles when they exit the domain.
//copied from HopperSelfTest
int main()
{
    
    //Problem parameters
    BidispersedChuteParameters parameters(0, 22, 0);
    HopperForAvalanche problem(parameters);
    //solve
    problem.solve();
    problem.write(std::cout, false);
}
