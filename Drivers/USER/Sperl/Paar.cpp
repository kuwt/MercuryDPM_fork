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

#include <Walls/AxisymmetricIntersectionOfWalls.h>
#include <Boundaries/CubeInsertionBoundary.h>
#include <Walls/TriangleWall.h>
#include "Mercury3D.h"
#include "Species/LinearViscoelasticSlidingFrictionSpecies.h"
#include "CMakeDefinitions.h"
using helpers::round;

class Paar : public Mercury3D {

    //simulation parameters

    ///height of the outer cylinder
    double cylinderHeight = 12e-2;
    ///filling height
    double fillingHeight = 5e-2;
    ///filling height
    double bladeHeight = 2e-2;
    ///radius of the outer cylinder
    double cylinderRadius = 5e-2;
    ///speed of shearing (mean shear rate: 2*pi*omega)
    double rpm = 6;
    ///number of revolutions
    double revolutions = 2.0;
    ///minimal particle radius; set in the main function
    double rMin;
    ///minimal particle radius; set in the main function
    double rMax;
    ///frame rate for output files (Hz)
    double frameRate = 50;

public:

    /** 
     * Constructor, sets all parameters except the particle and wall positions 
     * (this split is necessary for the code to work in parallel).
     * Note, since the parameters are set in the constructor, they can be overwritten in the main function.
     */
    Paar (double rMin_, double rMax_)
    : rMin(rMin_), rMax(rMax_) {
        setName("Paar");
        //define domain size
        setMin({-cylinderRadius,-cylinderRadius,0.0});
        setMax({cylinderRadius,cylinderRadius,fillingHeight});
        //define how to split up the domain if mpi is used
        splitDomain(DPMBase::DomainSplit::XYZ);
        //define species and time step
        setSpeciesAndTimeStep();
        //define gravity
        setGravity({0.0, 0.0, -9.8});
        //set final time
        setTimeMax(revolutions/rpm*60.0);
        // save data at given frame rate
        setSaveCount(1.0/getTimeStep()/frameRate);
        // write paraview output
        if (NUMBER_OF_PROCESSORS==1) setParticlesWriteVTK(true);
        wallHandler.setWriteVTK(FileType::MULTIPLE_FILES_PADDED);
    }

    /**
     * Set initial conditions for particles and walls
     */
    void setupInitialConditions() override {
        //set domain to the size that should be shown in the visualisation
        setMax({cylinderRadius,cylinderRadius,fillingHeight});
        //define mixer geometry
        setGeometry();
        //define particle properties
        setParticles();
    }

    /** 
     * Set species (material properties) and time step
     */
    void setSpeciesAndTimeStep() {
        //define a new material
        LinearViscoelasticSlidingFrictionSpecies s;
        //set density
        s.setDensity(2000);
        //set collision time much smaller than gravitational time scale of smallest particle, tg=sqrt(d/g)
        double mass = s.getMassFromRadius(rMin,speciesHandler);
        double collisionTime = 0.01*sqrt(2.0*rMin/9.8);
        //set stiffness and dissipation based on collision time and restitution coeff.
        s.setCollisionTimeAndRestitutionCoefficient(collisionTime,0,mass);
        //set friction
        s.setSlidingDissipation(2./7.*s.getDissipation());
        s.setSlidingStiffness(2./7.*s.getStiffness());
        s.setSlidingFrictionCoefficient(0.5);
        //add the material type to the handler
        speciesHandler.copyAndAddObject(s);
        //set the time step accordingly
        setTimeStep(collisionTime/25.);
        logger(INFO,"Time step: % s",getTimeStep());
    }

    /**
     * Set walls, i.e. the mixer geometry
     */
    void setGeometry() {
        // Set outer cylinder
        {
            AxisymmetricIntersectionOfWalls w;
            w.setSpecies(speciesHandler.getObject(0));
            w.setAxis({0,0,1});
            w.addObject({1, 0, 0}, {cylinderRadius, 0, 0});
            wallHandler.copyAndAddObject(w);
        }
        // Set flat base
        {
            InfiniteWall w;
            w.setSpecies(speciesHandler.getObject(0));
            w.set({0, 0, -1}, {0, 0, 0});
            wallHandler.copyAndAddObject(w);
        }
        // Add screw from STL-file; move and rotate it to the middle of the cylinder, at blade height; start rotation
        {
            Mdouble scaleFactor = 1e-3;
            Vec3D centerOfRotation {0,0.018,0.014};
            //Vec3D centerOfRotation {0,0,0};
            Vec3D velocity {0,0,0};
            Vec3D angularVelocity {0,0,2.0*constants::pi*rpm/60.0};
            //the STL file is in the source directory
            unsigned g = wallHandler.readTriangleWall(
                    getMercurySourceDir()+"/Drivers/USER/Sperl/159324-F.STL",
                    speciesHandler.getObject(0),
                    scaleFactor,
                    centerOfRotation,
                    velocity,
                    angularVelocity);
            for (auto& w: wallHandler) {
                if (w->getGroupId()==g) {
                    const Vec3D quarterTurn = {0,-2,0};
                    w->move(-centerOfRotation);
                    w->rotate(quarterTurn);
                    w->move({0,0,bladeHeight});
                }
            }
            for (auto w = wallHandler.begin(); w!=wallHandler.end(); w++) {
                if ((*w)->getGroupId()==g) {
                    auto t = static_cast<TriangleWall *>(*w);
                    if (t->getVertex(0).Z > cylinderHeight
                     && t->getVertex(1).Z > cylinderHeight
                     && t->getVertex(2).Z > cylinderHeight) {
                        w--;
                        wallHandler.removeObject(t->getIndex());
                    }
                }
            }
        }
        logger(INFO,"Inserted % walls",wallHandler.getSize());
    }

    /**
     * Set initial particle properties
     */
    void setParticles() {
        // Add an inner cylinder to avoid insertion into the rotor
        {
            AxisymmetricIntersectionOfWalls w;
            w.setSpecies(speciesHandler.getObject(0));
            w.setAxis({0,0,1});
            w.addObject({-1, 0, 0}, {3e-3, 0, 0});
            wallHandler.copyAndAddObject(w);
        }
        // Use spherical particles, of material-type 0
        SphericalParticle p;
        p.setSpecies(speciesHandler.getObject(0));
        // Use an insertion boundary to insert a specified volume of particles in a specified region.
        // Particles are instered at random locations, while avoiding overlaps.
        CubeInsertionBoundary b;
        b.setHandler(&boundaryHandler);
        // Why is there a loop here?
        // Initially, we try to insert all particles in the domain [xMin,xMax]x[yMin,yMax]x[zMin,zMax].
        // If not all particles could be inserted the zMax value is increased until all particles have been inserted.
        Vec3D min = getMin();
        Vec3D dif = getMax()-getMin();
        dif.Z = 0.1*(getMax().Z-getMin().Z);
        double insertedFraction;
        do {
            b.set(p, 10000, min, min + dif, {0, 0, 0}, {0, 0, 0});
            // Choose a volume of particles that fills the cylinder to the filling height when settled
            // (assuming a solid volume fraction of 0.64)
            b.setInitialVolume(0.9*0.64 * constants::pi * mathsFunc::square(cylinderRadius) * fillingHeight);
            // Insert particles
            b.checkBoundaryBeforeTimeStep(this);
            // Check whether all particles have been inserted. If not, increase zMax
            insertedFraction = b.getVolumeOfParticlesInserted()/b.getInitialVolume();
            logger(INFO,"Inserted % particles, max height %, insertedFraction %",
                    particleHandler.getNumberOfRealObjects(),min.Z+dif.Z,insertedFraction);
            min.Z += 0.5*dif.Z;
        } while (insertedFraction<1);
        // Do a final check whether all particles have been inserted.
        logger.assert_always(insertedFraction>1,"Could not insert all particles (fill fraction %)",insertedFraction);
        //remove that inner cylinder again
        wallHandler.removeLastObject();
    }

    void printTime() const override {
        logger(INFO,"t=%\ttmax=%\tene=%",round(getTime(),3),getTimeMax(),round(getKineticEnergy()/getElasticEnergy(),3));
    }

    void actionsAfterTimeStep() override {
        outputTorque();
    }

    /*
     * writes the torque applied by the blade into a file
     */
    void outputTorque() {
        // open a file for output
        static std::ofstream file(getName()+".torque");
        // write a header line
        if (getNumberOfTimeSteps()==0)
            file << "time\t torque\n";
        // collect the torque values from all triangle walls
        if (getNumberOfTimeSteps()%50!=0) return;
        Vec3D torque {0,0,0};
        unsigned g = wallHandler.getLastObject()->getGroupId();
        for (auto& w: wallHandler) {
            if (w->getGroupId()==g)
                torque += w->getTorque();
        }
        //write to file
        file << getTime() << '\t' << -torque.Z << '\n';
    }


};

/**
 * The main function, , 
 */
int main (int argc, char *argv[]) {
    // Read parameters from the command line, i.e. use './Paar -rMin 0.001' to set the particle radii to 1-1.6 mm 
    double rMin = helpers::readFromCommandLine<double>(argc,argv,"-rMin",1e-3);
    double rMax = 1.5*rMin;
    // Call the constructor
    Paar dpm(rMin,rMax);
    // removes old files with the same root name ("Paar.*"); be careful not to delete data you need!
    dpm.removeOldFiles();
    // Call the solve routine (starts with setupInitialConditions, then a time loop)
    dpm.solve();
    // Finish the programme
    return EXIT_SUCCESS;
}