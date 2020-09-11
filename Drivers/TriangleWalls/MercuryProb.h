//
// Created by cheng on 2/13/19.
//

#ifndef PROJECT_MERCURYPROB_H
#define PROJECT_MERCURYPROB_H

// Generic MercuryDPM header
#include <Mercury3D.h>
#include <Species/LinearViscoelasticSpecies.h>
#include "Walls/TriangleWall.h"

//=======start_namespace==========================================
/// Units:: mass:    1 ug = 1e-9 kg
///         length:  1 mm = 1e-3 m
///         time:    1 ms = 1e-3 s
///         Stress:  1 ug/(mm*ms^2) = 1 Pa
///         Density: 1 ug/mm^3 = 1 kg/m^3
//================================================================

namespace units {
    // mass: 1 ug =  kg
    Mdouble mUnit = 1e-9;
    // length: 1 mm = 1e-3 m
    Mdouble lUnit = 1e-3;
    // time: 1 ms = 1e-3 s
    Mdouble tUnit = 1e-3;

    // force: 1 uN = 1e-6 N
    inline Mdouble fUnit() { return mUnit * lUnit / pow(tUnit, 2); }

    // stiffness: 1 N/m = 1e-3 uN/mm
    inline Mdouble kUnit() { return fUnit() / lUnit; }

    // stress: ug/(mm*ms^2) = 1 Pa
    inline Mdouble sigUnit() { return mUnit / (lUnit * pow(tUnit, 2)); }

    // density: 1 ug/mm^3 = 1 kg/m^3
    inline Mdouble rhoUnit() { return mUnit / pow(lUnit, 3); }

    // velocity
    inline Mdouble velUnit() { return lUnit / tUnit; }

    // acceleration
    inline Mdouble accUnit() { return lUnit / pow(tUnit, 2); }

    // simulation name
    std::string name;
}

//=============begin_Mercury_problem====================================
/// Problem class for a single particle bouncing on a "beam" structure.
//======================================================================

class MercuryProblem : public Mercury3D {
public:

    MercuryProblem() = default;

    void setupMercuryProblem(const char *name, const unsigned &dim, const double &rveSize, const unsigned &rve,
                             const Mdouble &g, const Mdouble &tMax, const unsigned &nWrite) {
        units::name = name;
        setName(name);
        setSystemDimensions(dim);
        nParticlesPerRVE1D = rve;
        RVESize = rveSize;
        setGravity(Vec3D(0.0, 0.0, g));
        setXMin(-0.5 / units::lUnit);
        setXMax(0.5 / units::lUnit);
        setYMin(-5.0 / units::lUnit);
        setYMax(5.0 / units::lUnit);
        setZMin(-0.5 / units::lUnit);
        setZMax(5.0 / units::lUnit);
        setSpeciesProperties(0);
        setTimeMax(tMax);
//        setSaveCount(1);
        setSaveCount(unsigned(getTimeMax() / getTimeStep() / nWrite));
        setParticlesWriteVTK(true);
        setWallsWriteVTK(true);
    }

    void setSpeciesProperties(const unsigned &flag) {
        /// set default Young's modulus, Poisson's ratio, friction coefficient and dissipation
        radius = 0.25 * RVESize / nParticlesPerRVE1D;
        eps = 1.0;

        /// set material properties
        const Mdouble density = 2500.0 / units::rhoUnit();
        auto species = speciesHandler.copyAndAddObject(LinearViscoelasticSpecies());
        species->setDensity(density);
        Mdouble mass = species->getMassFromRadius(radius);
        /// k_n = E_n * A / L
        species->setStiffnessAndRestitutionCoefficient(1e8 * pow(2 * radius, 2) / (2 * radius), eps, mass);
        tc = species->getCollisionTime(mass);
        setTimeStep(tc * 0.02);
    }

    void setupInitialConditions() override;

    void createWalls();

    TriangleWall *createTriangleWall(std::array<Vec3D, 3> vertex) {
        TriangleWall wall;
        auto species = speciesHandler.getObject(0);
        wall.setSpecies(species);
        wall.setVertices(vertex[0], vertex[1], vertex[2]);
        auto w = wallHandler.copyAndAddObject(wall);
        return w;
    }

    void actionsAfterTimeStep() override {
/*
        for (auto inter : interactionHandler) {
            std::cout << inter->getI()->getGroupId() << inter->getI()->getName() << " " << inter->getP()->getGroupId()
                      << inter->getP()->getName() << std::endl;
        }
*/
        // check if equilibrium is broken
        Vec3D f = interactionHandler.getLastObject()->getForce();
        Vec3D a = particleHandler.getLastObject()->getForce();
        if (a.getZ() > 1.0e-5*f.getZ())
        {
            logger(INFO,"interaction force: % % % , body force: % % %",f.getX(),f.getY(),f.getZ(),a.getX(),a.getY(),a.getZ());
        }

        // check if each particle-wall interaction is unique
        unsigned n = 0;
        for (auto w : wallHandler) {
            for (auto inter : w->getInteractions())
            {
                Vec3D f = inter->getForce();
                logger(INFO, "wallID#% has force %, % , %", w->getId(), f.getX(), f.getY(), f.getZ());
            }
            n += w->getInteractions().size();
        }
        if (n>particleHandler.getNumberOfObjects()) logger(INFO, "% walls have interactions with the particle",n);
    }

private:

    Mdouble radius = 0.5;
    Mdouble eps = 1;
    Mdouble tc = 0.01;

    double RVESize = 0;
    unsigned nParticlesPerRVE1D = 1;
    Mdouble posZ0 = 0;
    std::vector<TriangleWall *> listOfTriangleWalls;

};

#endif //PROJECT_MERCURYPROB_H

