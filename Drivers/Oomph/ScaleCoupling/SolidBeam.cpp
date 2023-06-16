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

#include "Oomph/SolidProblem.h"
#include "Math/Helpers.h"

/// Defines a SolidProblem of element type RefineableQDPVDElement<3,2>. Add functionality to write output
class Beam : public SolidProblem<RefineableQDPVDElement<3, 2>>
{
    double length = 20.0; // length of beam
    double distance = 0.4; // distance between particles
    double bulkDensity = 1309; // bulk density (assuming a cubic packing)
    double elasticModulus = 1e8; //
    std::function<Vec3D(double)> velocityAtBoundary =
        [] (double) {return Vec3D(0.1,0,0);}; // velocity of particle at left boundary

public:
    Beam() {
        setName("SolidBeam");
        removeOldFiles();
        setElasticModulus(elasticModulus);
        setDensity(bulkDensity);
        setSolidCubicMesh(length/distance, 2, 2,
                          0, length, -distance, distance, distance, 3.0*distance);
        pinBoundaries({Beam::Boundary::X_MIN});
        setNewtonSolverTolerance(1e-2);
        prepareForSolve();
        linear_solver_pt()->disable_doc_time();
        // simulate two oscillations, check energy balance
        double waveSpeed = sqrt(getElasticModulus()/getDensity());
        double propagationTime = length/waveSpeed;
        double timeMax = 1.25*propagationTime;

        solveUnsteady(timeMax, 0.01*timeMax, 1);
    }

private:

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

        // set new position
        static Vec3D positionAtBoundary {0,0,0};
        Vec3D velocity = velocityAtBoundary(getOomphTime());
        positionAtBoundary += velocity * getOomphTimeStep();
        for (int n = 0; n < mesh_pt()->nboundary_node(Boundary::X_MIN); ++n) {
            SolidNode* node_pt = dynamic_cast<SolidNode*>(mesh_pt()->boundary_node_pt(Boundary::X_MIN,n));
            node_pt->x(0) = positionAtBoundary.X;
        }
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
 * Simulates a cubic rod (20mm x 2mm x 2mm) of density 2500 kg/m^3, elastic modulus 100 MPa, fixed on left side, bent by gravity
 */
int main()
{
    // Solve the problem
    Beam();
    //system("head -n 100 SolidBeam.out");
    return 0;
}