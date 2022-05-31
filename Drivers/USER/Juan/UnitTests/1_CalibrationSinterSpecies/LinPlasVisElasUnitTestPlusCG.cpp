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

#include <Mercury3D.h>
#include "Walls/InfiniteWall.h"
#include <Species/ThermalSinterLinFrictionReversibleAdhesiveSpecies.h>
#include "CG/Functions/Gauss.h"
#include "StatisticsVector.h"
#include "CG/CG.h"
#include <iostream>
#include <iomanip>
///This code tests the linear plastic viscoelastic behavior of two particles in sintering.
class calibrationUnitTest : public DPMBase {
public:

    explicit calibrationUnitTest(Mdouble radius) {
        //-----------------nu
        //Global parameters
        std::string r = helpers::to_string(radius);
        setName("LinearPlasticViscoElasticTestCG");

        setFileType(FileType::ONE_FILE);
        setSaveCount(10000);
        setTimeMax(7);
        //setParticlesWriteVTK(true);

        //setParticlesWriteVTK(true);
        //setWallsWriteVTK(FileType::MULTIPLE_FILES);

        setParticleDimensions(3);
        setSystemDimensions(3);
        setGravity(Vec3D(0.0, 0.0,0.0));
 
        //-------------------
        //Boundary parameters
        setXMax(2.0 * radius);
        setYMax(radius);
        setZMax(radius);
        setXMin(-getXMax());
        setYMin(-getYMax());
        setZMin(-getZMax());

        //-------------------
        //Setup parameters:
        const Mdouble density = 1000;
        const Mdouble k1 = 0.01 * radius; //[Stiffness depends on particle radius]
        const Mdouble restitutionCoefficient = 0.1;
        const Mdouble pDepth = 0.85;

        const Mdouble YoungM = 2.95e9; //[Pa] Young's Modulus for polystyrene

        //-------------------
        //Species:
        ThermalSinterLinFrictionReversibleAdhesiveSpecies sf;
        //SinterSpecies sf;
        sf.setHandler(&speciesHandler);
        sf.setDensity(density);

        const Mdouble mass = sf.getMassFromRadius(radius);
        sf.setStiffnessAndRestitutionCoefficient(k1, restitutionCoefficient, mass);
        sf.setPlasticParameters(k1, 10.0 * k1, k1, pDepth);

        const Mdouble collisionTime = sf.getCollisionTime(mass);
        logger(INFO, "Collision time %", collisionTime);

        //-------------------
        //Sinter parameters
        sf.setSinterType(SINTER_APPROACH::VISCOELASTIC_CONTACT);  //FRENKEL OR VISCOELASTIC_CONTACT
//        sf.setSinterRate(0.16);
        sf.setSinterAdhesion(0.0016 * k1);

//        -------------------
//        Parameters to be used in contact growth driven by adhesive intersurface forces and accomodated by
//        viscoelastic deformation.
        sf.setComplianceZero(1/(2*YoungM)); //Book: Adhesive Particles
        sf.setSurfTension(0.047); // [N/m] or [J/m^2] surface energy density or surface tension. Paper:Particle size limits for sintering polymer colloids. S Mazur.
        sf.setFluidity(4.0e-6);
        sf.setSeparationDis(2.8e-5);

        auto species = speciesHandler.copyAndAddObject(sf);

        setTimeStep(0.00007 * collisionTime);
        //-------------------
        //Particle properties:
        SphericalParticle P0, P1;
        //SphericalParticle P0, P1;
        P0.setSpecies(species);
        P1.setSpecies(species);

        P0.setRadius(radius);
        P1.setRadius(radius);

        P0.setPosition(Vec3D(-(1 - 1e-15) * radius, 0, 0));
        P1.setPosition(-P0.getPosition());

        particleHandler.copyAndAddObject(P0);
        particleHandler.copyAndAddObject(P1);

        //------------------
        //set walls
        InfiniteWall w0;
        w0.setSpecies(speciesHandler.getObject(0));
        w0.set(Vec3D(-1, 0, 0), Vec3D(getXMin(), 0, 0));
        wallHandler.copyAndAddObject(w0);
        w0.set(Vec3D(1, 0, 0), Vec3D(getXMax(), 0, 0));
        wallHandler.copyAndAddObject(w0);
        w0.set(Vec3D(0, -1, 0), Vec3D(0, getYMin(), 0));
        wallHandler.copyAndAddObject(w0);
        w0.set(Vec3D(0, 1, 0), Vec3D(0, getYMax(), 0));
        wallHandler.copyAndAddObject(w0);
        w0.set(Vec3D(0, 0, -1), Vec3D(0, 0, getZMin()));
        wallHandler.copyAndAddObject(w0);
        w0.set(Vec3D(0, 0, 1), Vec3D(0, 0, getZMax()));
        wallHandler.copyAndAddObject(w0);

    }

    Mdouble getMeanRelativeContactRadius() const
    {
        Mdouble meanOverlap = interactionHandler.getMeanOverlap();
        Mdouble meanRadius = particleHandler.getMeanRadius();
        logger(INFO," a/r %",sqrt(2.0*meanOverlap/meanRadius));
        return sqrt(2.0*meanOverlap/meanRadius);

    }

    void printTime() const override {
        std::cout << "t=" << std::setprecision(3) << std::left << std::setw(6) << getTime()
                  << ", tmax=" << std::setprecision(3) << std::left << std::setw(6) << getTimeMax()
                  << ", r=" << std::setprecision(3) << std::left << std::setw(6)
                  << " Ene " << std::setw(12) << getKineticEnergy()/getElasticEnergy()
                  << " X/a " << std::setw(12) << getMeanRelativeContactRadius()
                  //<< ", Ekin=" << std::setprecision(3) << std::left << std::setw(6) << getKineticEnergy()/getElasticEnergy()
                  << std::endl;
        //std::cout.flush();
    }
};

int main(int argc UNUSED, char *argv[] UNUSED)
{
    calibrationUnitTest sp0(3.216e-5);

    // GrainLearning results over time: For example, to lot the stress over time, or density.
//    CG<CGCoordinates::XZ> cg;
//    cg.setNX(30);
//    cg.setNZ(30);
//    cg.setWidth(0.15);
//    cg.statFile.setSaveCount(1);
//    cg.statFile.setName(sp0.getName() + ".GaussXZ_.stat");
//    sp0.cgHandler.copyAndAddObject(cg);

    sp0.solve();

    //This helper is to see the neck growth
    std::cout << "Execute 'gnuplot LinearPlasticViscoElasticTestCG.gnu' to view output" << std::endl;
    helpers::writeToFile("LinearPlasticViscoElasticTestCG.gnu",
                         "set xlabel 'time [s]'\n"
                         "set ylabel 'a(t)/R'\n"
                         "plot 'LinearPlasticViscoElasticTestCG.fstat' u ($1):(sqrt($7/3.216e-5)) title 'DEM simulation'  with lines linestyle 2\n"
                         "replot 0.3*x**0.5 title 'Frenkel' with lines linestyle 1"
    );

    logger(INFO,"Execute 'source LinearPlasticViscoElasticTestCG.sh' to get coarse-grained statistics of the last time step");
    helpers::writeToFile("LinearPlasticViscoElasticTestCG.sh","../../../../MercuryCG/fstatistics LinearPlasticViscoElasticTestCG -stattype XZ -w 1.0e-5 -h 0.1e-5 -tmin 5.0 -tmax 5.1 -o LinearPlasticViscoElasticTestCG.XZ.stat\n");

    logger(INFO,"Run 'LinearPlasticViscoElasticTestCG.m' in MATLAB/octave to visualise the statistical output");
    helpers::writeToFile("LinearPlasticViscoElasticTestCG.m","addpath('../../../../MercuryCG/')\n"
                                   "data = loadStatistics('LinearPlasticViscoElasticTestCG.XZ.stat');\n"
                                   "colormap(1-gray)\n"
                                   "contourf(data.x,data.z,data.Density,20,'EdgeColor','none')\n"
                                   "c = colorbar\n"
                                   "c.Label.String = '\\rho';\n"
                                   "title('Density')\n"
                                   "xlabel('x')\n"
                                   "ylabel('z');\n"
                                   "axis equal\n"
                                   "%%\n"
                                   "particles=importdata('LinearPlasticViscoElasticTestCG.data',' ',12);\n"
                                   "Cell = particles{1}(1,1);\n"
                                   "NumPart= Cell.N;\n"
                                   "Pradius = Cell.Radius;\n"
                                   "Position = Cell.Position;\n"
                                   "a=linspace(0,2*pi,40);\n"
                                   "xCircle = sin(a);\n"
                                   "zCircle = cos(a);\n"
                                   "hold on;\n"
                                   "for i=1:length(Pradius)\n"
                                   "  plot(Position(i,1) + Pradius(i)*xCircle,Position(i,2)+Pradius(i)*zCircle,'Color',.8*[1 1 1])\n"
                                   "end\n"
                                   "hold off");

//    logger(INFO,"Run 'LinearPlasticViscoElasticTestCG.m' in MATLAB/octave to visualise the statistical output");
//    helpers::writeToFile("LinearPlasticViscoElasticTestCG.m", "close all; clear; clc;\n"
//           "addpath('../../../../MercuryCG/')\n"
//           "data = loadStatistics('LinearPlasticViscoElasticTestCG.stat');\n"
//           "colormap(1-gray);\n"
//           "X = data.x;\n"
//           "Z = data.z;\n"
//           "Field = data.Density;\n"
//           "contourf(X,Z,Field,20,'EdgeColor','none')\n"
//           "c = colorbar\n"
//           "c.Label.String = '\\rho';\n"
//           "title('Density')\n"
//           "xlabel('x')\n"
//           "ylabel('z');\n"
//           "axis equal\n"
//           "%%\n"
//           "addpath('../../../../../../Matlab/');\n"
//           "particles=read_data('LinearPlasticViscoElasticTestCG.data');\n"
//           "Cell = particles{1}(1,1);\n"
//           "NumPart= Cell.N;\n"
//           "Pradius = Cell.Radius;\n"
//           "Position = Cell.Position;\n"
//           "a=0:0.1:2*pi;\n"
//           "xCircle = sin(a);\n"
//           "zCircle = cos(a);\n"
//           "hold on;\n"
//           "for i=1:length(Pradius)\n"
//           "  plot(Position(i,1) + Pradius(i)*xCircle,Position(i,2)+Pradius(i)*zCircle,'Color',.8*[1 1 1])\n"
//           "end\n"
//           "hold off\n"
//           "figure(2)\n"
//           "%%\n"
//           "s = surf(X,Z,Field,'FaceAlpha',0.3);\n"
//           "s =colormap(bone)\n"
//           "colorbar;\n"
//           "%s.LineStyle= ':';\n"
//           "%s.FaceColor = 'interp';\n"
//           "%s.FaceLighting = 'gouraud';\n"
//           "title('Density field')\n"
//           "xlabel('x')\n"
//           "ylabel('z')\n"
//           "zlabel('\rho')\n"
//           "zlim([-0 3500])\n"
//           "ylim([-4 4]*10^-5)\n"
//           "xlim([-6 6]*10^-5)\n"
//           "v = [-3 -8 4];\n"
//           "[caz,cel] = view(v)"
//           );
    return 0;

}
///Todo{Add MercuryCG during simulation}