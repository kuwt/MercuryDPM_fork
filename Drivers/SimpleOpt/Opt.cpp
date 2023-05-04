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

// SimpleOpt feature, opt.cpp: Example of the driver routine being run by the
// external scipy optimizer

#include <Mercury3D.h>
#include <Species/Species.h>
#include <Species/LinearViscoelasticSpecies.h>
#include<CMakeDefinitions.h>
#include <cstdio>
#include <iostream>
#include <fstream>
#include <memory>
#include <stdexcept>
#include <string>
#include <array>
#include "Walls/InfiniteWall.h"
#include <cmath>

class SimpleOpt : public Mercury3D
{

public:
    
    void setupInitialConditions() override {

        setName("SimpleOpt");
        setSystemDimensions(3);
        setGravity(Vec3D(0.0, 0.0, -10.0));
        setXMax(1.0);
        setXMin(-1.0);
        setYMax(1.0);
        setYMin(-1.0);
        setZMax(1.0);
        setZMin(-0.05);
        setTimeMax(0.8);



    }
};

// Helper function to implement extra OS commands after the driver code is done
std::string ExecCommand(const char* cmd) {
    std::array<char, 128> buffer;
    std::string result;
    std::unique_ptr<FILE, decltype(&pclose)> pipe(popen(cmd, "r"), pclose);
    if (!pipe) {
        std::string ret = "Could not run an extra OS command: "; ret.append(cmd);
        throw std::runtime_error(ret);
    }
    while (fgets(buffer.data(), buffer.size(), pipe.get()) != nullptr) {
        result += buffer.data();
    }
    return result;
}

int main(int argc, char* argv[])
{
    // In this simple example we optimize two angular parameters (theta and phi in a standard spherical coord system)
    // to extremize the distance the projectile travels in a given direction

    int NumParams = 2;

    // Get parameters passed through the command line
    std::vector <Mdouble> params(0);
    for (int i = 0; i<NumParams; i++) {
        std::string a;
        Mdouble param = stod(helpers::readFromCommandLine(argc, argv, "-p"+std::to_string(i), a));
        params.push_back(param);
        std::cout<<params[i]<<std::endl;
    }

    // Problem setup
    SimpleOpt problem;

    // Species
    LinearViscoelasticSpecies species;
    species.setDensity(2500.0);
    species.setStiffness(258.5);
    species.setDissipation(0.0);
    problem.speciesHandler.copyAndAddObject(species);

    // Tabletop
    problem.wallHandler.clear();
    InfiniteWall w0;
    w0.setSpecies(problem.speciesHandler.getObject(0));
    w0.set(Vec3D(0.0, 0.0, -1.0), Vec3D(0, 0, problem.getZMin()));
    problem.wallHandler.copyAndAddObject(w0);

    // Projectile
    SphericalParticle p0;
    p0.setSpecies(problem.speciesHandler.getObject(0));
    p0.setRadius(0.05); // sets particle radius
    p0.setPosition(Vec3D(0,0,0));
    p0.setVelocity(Vec3D(0.1 * problem.getXMax(), 0.1 * problem.getYMax(), 0.1 * problem.getZMax())); // sets particle position

    Mdouble VelMag = 4;
    Mdouble theta = params[0]; // Controls the inclination angle of a cannon
    Mdouble phi = params[1];  // Controls azimutal orientation of a cannon

    // Set the initial velocity of the projectile
    p0.setVelocity(Vec3D(VelMag * sin(theta) * cos(phi), VelMag * sin(theta) * sin(phi), VelMag * cos(theta)));
    problem.particleHandler.copyAndAddObject(p0);

    // Time integration and log parameters
    problem.setSaveCount(100);
    problem.dataFile.setFileType(FileType::ONE_FILE);
    problem.restartFile.setFileType(FileType::ONE_FILE);
    problem.fStatFile.setFileType(FileType::NO_FILE);
    problem.eneFile.setFileType(FileType::NO_FILE);
    logger(INFO, "run number: %", problem.dataFile.getCounter());
    problem.setXBallsAdditionalArguments("-solidf -v0");
    problem.setTimeStep(0.005 / 50.0); // (collision time)/50.0
    problem.solve();

    // Paraview data
    std::cout<<ExecCommand("rm -rf paraview_SimpleOpt/")<<std::endl;
    std::cout<<ExecCommand("mkdir paraview_SimpleOpt/")<<std::endl;
    std::cout<<ExecCommand("../../Tools/data2pvd SimpleOpt.data paraview_SimpleOpt/SimpleOpt")<<std::endl;

    // Return the functional via the text file
    std::ofstream funct;  funct.open ("functional.txt");
    Vec3D pos = problem.particleHandler.getLastObject()->getPosition();
    Mdouble fun = 1 / (pos.X * pos.X); // Minimization of this optimizes for the longest projectile travel along X direction
    funct << fun;
    funct.close();

    return 0;
}
