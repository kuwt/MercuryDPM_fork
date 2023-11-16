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

#include "Math/Helpers.h"
#include "Species/LinearViscoelasticReversibleAdhesiveSpecies.h"
#include "Oomph/ScaleCoupling/ScaleCoupledSolid.h"
#include "Oomph/RefineableQDPVDElement.h"
#include "Oomph/SolidProblem.h"

/**
 * Define a coupled problem
 */
class ScaleCoupledBeam : public ScaleCoupledSolid<RefineableQDPVDElement<3, 2>>
{
    double length = 20.0; // length of beam
    double distance = 0.4; // distance between particles
    double bulkDensity = 1309; // bulk density (assuming a cubic packing)
    double elasticModulus = 1e8; // =stiffness/distance
    double overlapLength = 10.0*distance; // length of overlap zone
    double penalty = 8e9;
    double velocity = 1e-1;

public:
    
    ScaleCoupledBeam () {
        //set name
        setName("ScaleCoupledBeam");
        //remove existing output files
        removeOldFiles();
        
        // setup steps
        setupOomph();
        setupMercury(); //sets timeMax

        // set weight function
        auto couplingWeight = [this] (double x,double y,double z) {
            if (x<0.5*(length-overlapLength)) {
                // DEM region
                return 1.0;
            } else if (x>0.5*(length+overlapLength)) {
                // FEM region
                return 0.0;
            } else {
                // coupled region
                //return (0.5*(length+overlapLength)-x)/overlapLength;
                double unit = (0.5*(length+overlapLength)-x)/overlapLength;
                return 0.5+0.4*std::sin(constants::pi*(unit-0.5));
            }
        };
        setCouplingWeight(couplingWeight);
        setPenalty(penalty);

        // setup time
        double waveSpeed = sqrt(elasticModulus/bulkDensity);
        double propagationTime = length/waveSpeed;
        setTimeMax(1.25*propagationTime);
        setOomphTimeStep(getTimeStep());
        //setOomphTimeStep(propagationTime/100.0);
        setSaveCount(getTimeMax()/getTimeStep()/100.0);
        logger(INFO,"TimeMax: %, nTimeSteps %", getTimeMax(), getTimeMax()/getTimeStep());

        // solve
        writeToVTK();
        solveScaleCoupling();
        saveSolidMesh();
    }
    
    void setupOomph() {
        setElasticModulus(elasticModulus);
        setDensity(bulkDensity);
        setSolidCubicMesh(round(0.5*(length+overlapLength)/distance), 2, 2,
                          0.5*(length-overlapLength), length, -distance, distance, -distance, distance);
        //pinBoundaries({Beam::Boundary::X_MIN});
        setNewtonSolverTolerance(1e-2);
        prepareForSolve();
        linear_solver_pt()->disable_doc_time();
        //disable_info_in_newton_solve();
    }

    void setupMercury()
    {
        // positive overlap means adhesive forces
        double overlap = 0.0*distance; //overlap between particles
        double radius = 0.5*(distance+overlap); // particle radius
        double particleDensity = bulkDensity*mathsFunc::cubic(distance) / (constants::pi/6.0*mathsFunc::cubic(2.0*radius));
        double stiffness = elasticModulus*distance; // particle stiffness

        setDomain(Vec3D(-radius,-radius,0),Vec3D(radius,radius,0.5*(length+overlapLength)));
        setParticlesWriteVTK(true);

        auto species = speciesHandler.copyAndAddObject(LinearViscoelasticReversibleAdhesiveSpecies());
        species->setDensity(particleDensity);
        species->setStiffness(stiffness);
        species->setAdhesionStiffness(stiffness);
        species->setAdhesionForceMax(stiffness*overlap);

        SphericalParticle p(species);
        p.setRadius(radius);
        //p.setVelocity(Vec3D(1,0,0));
        double dx = distance;
        auto n = (unsigned)(0.5*(length+overlapLength)/dx);
        for (int i = 0; i<n; ++i) {
            for (int j = 0; j<2; ++j) {
                for (int k = 0; k<2; ++k) {
                    p.setPosition(Vec3D(1+2*i, j==0?-1:1, k==0?-1:1)*0.5*distance);
                    particleHandler.copyAndAddObject(p);
                    if (i==0) {
                        particleHandler.getLastObject()->setVelocity(Vec3D(velocity,0,0));
                        //particleHandler.getLastObject()->fixParticle();
                        //particleHandler.getLastObject()->setPrescribedVelocity(velocityAtBoundary);
                    }
                }
            }
        }

        setTimeStep(0.5*species->computeTimeStep(p.getMass()));
    }

    /// Write header of output file
    void actionsBeforeSolve() override {
        helpers::writeToFile(getName()+".gnu", "set key autotitle columnheader\n"
                                               "p 'SolidBeamUnsteady.out' u 1:3 w l, '' u 1:4 w l, '' u 1:5 w l, '' u 1:($3+$4+$5) t 'totalEnergy' w l");
        out.open(getName()+".out");
        out << "time deflection elasticEnergy kineticEnergy gravEnergy\n";
    }

    /// Write header of output file
    void actionsAfterSolve() override {
        out.close();
    }

    /// Each time step, compute deflection, elastic, kinetic and gravitational energy, and write to output file
    void actionsBeforeOomphTimeStep() override {
        double mass, elasticEnergy, kineticEnergy;
        Vector<double> com(3), linearMomentum(3), angularMomentum(3);
        getMassMomentumEnergy(mass, com, linearMomentum, angularMomentum, elasticEnergy, kineticEnergy);
        static double comZ0 = com[2];
        double gravEnergy = 9.8*mass*(com[2]-comZ0);
        out << getOomphTime() << ' ' << getBeamDeflection() << ' ' << elasticEnergy << ' ' << kineticEnergy << ' ' << gravEnergy << std::endl;
        std::cout << getOomphTime() << ' ' << getBeamDeflection() << ' ' << elasticEnergy << ' ' << kineticEnergy << ' ' << gravEnergy << '\n';
    }

    /// Computes beam deflection
    double getBeamDeflection() const {
        std::array<double, 3> min, max;
        getDomainSize(min, max);

        Vector<double> xi(3);
        xi[0] = max[0];
        xi[1] = 0.5 * (max[1] + min[1]);
        xi[2] = 0.5 * (max[2] + min[2]);
        return getDeflection(xi, 2);
    }

    /// output file stream
    std::ofstream out;
};

/**
 * Measure bag deformation loaded by body force
 */
int main()
{
    ScaleCoupledBeam problem;
    return 0;
}
