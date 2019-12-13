//
// Created by paolo on 21-9-19.
//

#ifndef MERCURY_SINGLECLUSTERCOMPRESSION_H
#define MERCURY_SINGLECLUSTERCOMPRESSION_H

#endif //MERCURY_SINGLECLUSTERCOMPRESSION_H


#include <iostream>
#include <Species/LinearPlasticViscoelasticFrictionSpecies.h>
#include <Mercury3D.h>
#include "BaseCluster.h"
#include <Walls/InfiniteWall.h>

class SingleClusterCompression : public Mercury3D{

public:

    SingleClusterCompression();

    SingleClusterCompression(int nP, Mdouble rP, Mdouble sDP, Mdouble dP, Mdouble rC, Mdouble pDM,
                             Mdouble sFC, Mdouble lS, Mdouble uSM, Mdouble cS, Mdouble rTS);

    ~SingleClusterCompression() override;


private:


    void setupInitialConditions() override;

    void actionsAfterTimeStep() override;

    void actionsAfterSolve() override;

    void printTime() const override;

    void makeWalls();

    void makeIsoCFile();

    void moveWalls();

    void computeData();

    void writeIsoCFile();

private:

    //! Geometry
    Mdouble radiusParticle;
    Mdouble boxSize;
    InfiniteWall* wall1;
    InfiniteWall* wall2;

    //! Species and cluister
    LinearPlasticViscoelasticFrictionSpecies* species;
    Mdouble restitutionCoefficient;
    Mdouble penetrationDepthMax;
    Mdouble densityParticle;
    Mdouble slidingFrictionCoefficient;
    Mdouble loadingStiffness;
    Mdouble unLoadingStiffnessMax;
    Mdouble cohesionStiffness;
    Mdouble relativeTangentialStiffness;
    Mdouble smallestMass = densityParticle * 4 * constants::pi * pow(radiusParticle*(1-sizeDispersityParticle), 3) / 3;
    Mdouble collisionTimeSmallestMass = sqrt(smallestMass * ( pow(constants::pi, 2) + pow(log(restitutionCoefficient), 2) ) / ( 2 * loadingStiffness ) );

    //BaseClusterInsertionBoundary cluster1;
    int nParticles;
    Mdouble sizeDispersityParticle;

    //! Name and output
    std::ostringstream name;
    std::ofstream isoCFile;
    int fileOutputCount;
    Mdouble fileOutputTime;

    //! Data of interest
    Mdouble meanForceOnInteraction_;
    Mdouble relativeOverlap_;
    Mdouble minRelativeOverlap_;
    Mdouble meanRelativeOverlap_;
    Mdouble maxRelativeOverlap_;
    Mdouble previousForceValue;
    Mdouble kEq;
    Mdouble forceL;
    Mdouble forceU;
    std::vector<Mdouble> record;
    Mdouble t0;

    //! Computation
    Mdouble wallRelativePosition;
    Mdouble initialPositionWall;
    Mdouble recordPosition;
    bool firstTime;

};
