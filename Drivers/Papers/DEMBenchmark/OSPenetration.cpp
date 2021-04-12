//Copyright (c) 2013-2021, The MercuryDPM Developers Team. All rights reserved.
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
#include "MercuryOS.h"
#include <Walls/TriangleWall.h>

/**
 * Penetration Test. Three subtests can be performed, with a bed made of 25K, 50K, or 100K particles
 * An initial placement of particles is read from the file Penetration/Particles*.txt. This placement is generated
 * by filling a box with 25000, 50000 and 100000 particles of material M1 with a diameter of 2 mm. Afterwards, a steel
 * particle with diameter of 20 mm and initial velocity of 5 m/s should be placed around 15 cm above the bottom plate.
 * The initial placement of the walls is read from the file Penetration/walls.txt. One can change to a smooth wall
 * implementation (non-triangulated) by setting useMercuryWalls_ to true.
 * The time-dependent vertical position of the steel particle (above the bottom plate) is outputted to the ene file.
 * The simulation time step should be equal to 5e-7 s for a time period of 0.1 seconds.
 */
class Penetration : public MercuryOS
{
    /// whether 25k, 50k, or 100k particles are used for the bed
    std::string size_ = "25K";

public:
    
    // returns the variable size_
    std::string getSize() const { return size_; }
    
    // sets the variable size_
    void setSize(std::string size)
    {
        logger.assert_always(size == "25K" || size == "50K" || size == "100K", "Size argument % not valid", size);
        size_ = size;
    }
    
    // used to set the initial conditions of the particles, walls, species, etc
    void setupInitialConditions() override
    {
        // name of the output files
        setName("Penetration" + getSize());
    
        // turn on gravity
        setGravity({0, 0, -9.8});
        
        // set time step and maximum simulation time
        setTimeStep(5e-7);
        setTimeMax(0.1);
    
        // set output frequency
        setSaveCount(static_cast<unsigned>(0.001 / getTimeStep()));
    
        // remove files from previous run
        removeOldFiles();
    
        // determine which output files to write
        if (writeOutput()) {
            // write ene, data, fstat, restart and vtu files
            setParticlesWriteVTK(true);
            setWallsWriteVTK(FileType::MULTIPLE_FILES);
            fStatFile.writeFirstAndLastTimeStep();
            restartFile.writeFirstAndLastTimeStep();
        } else {
            // only .ene files are written
            setFileType(FileType::NO_FILE);
            eneFile.setFileType(FileType::ONE_FILE);
        }
    
        //use a shortened simulation in test mode
        if (test()) {
            setSaveCount(80);
            setTimeMax(800*getTimeStep());
            logger(INFO,"Test mode, reduced timeMax to %",getTimeMax());
        }
        
        // set domain for visualisation
        setMax(Vec3D(0.05, 0.03, 0.2 - 0.090243));
        setMin(Vec3D(-0.05, -0.03, -0.090243));
    
        // define the material properties of M1, M2, steel (see MercuryOS.h)
        setMaterialProperties();
    
        // read particle positions and radii from file
        {
            // open file
            std::string fileName = "Penetration/" + getSize() + "Particles.txt";
            std::ifstream file(fileName);
            logger.assert_always(file.is_open(), "File % could not be opened", fileName);
            // read header line
            std::string header;
            getline(file, header);
            // read particle positions and radii
            SphericalParticle particle;
            particle.setSpecies(m1);
            double rad;
            Vec3D pos;
            while (file >> rad >> pos) {
                particle.setRadius(rad);
                particle.setPosition(pos);
                particleHandler.copyAndAddObject(particle);
            }
            //there is some weird last line in the file that does not contain a proper particle
            particleHandler.removeLastObject();
            logger(INFO, "Read % particles from %", particleHandler.getSize(), fileName);
            
            // add steel particle
            particle.setSpecies(steel);
            particle.setRadius(0.01);
            particle.setPosition(Vec3D(0, 0, 0.06));
            particle.setVelocity(Vec3D(0, 0, -5));
            particleHandler.copyAndAddObject(particle);
        }
        
        // add walls
        if (useMercuryWalls()) {
            // use smooth walls
            InfiniteWall wall;
            wall.setSpecies(steel);
            wall.set(Vec3D(1, 0, 0), getMax());
            wallHandler.copyAndAddObject(wall);
            wall.set(Vec3D(-1, 0, 0), getMin());
            wallHandler.copyAndAddObject(wall);
            wall.set(Vec3D(0, 1, 0), getMax());
            wallHandler.copyAndAddObject(wall);
            wall.set(Vec3D(0, -1, 0), getMin());
            wallHandler.copyAndAddObject(wall);
            wall.set(Vec3D(0, 0, -1), getMin());
            wallHandler.copyAndAddObject(wall);
        } else {
            // read triangle walls from file
            // open file
            std::string fileName = "Penetration/walls.txt";
            std::ifstream file(fileName);
            logger.assert_always(file.is_open(), "File % could not be opened", fileName);
            // read header line
            std::string header;
            getline(file, header);
            // read triangle vertices
            TriangleWall wall;
            wall.setSpecies(steel);
            Vec3D a, b, c;
            while (file >> a >> b >> c) {
                wall.setVertices(a, b, c);
                wallHandler.copyAndAddObject(wall);
            }
            logger(INFO, "Read % walls from %", wallHandler.getSize(), fileName);
        }
    }
    
    // Write requested output to the ene file
    void writeEneTimeStep(std::ostream &os) const override
    {
        if (eneFile.getCounter() == 1) os << "Time Position\n";
        // time and z-position of steel particle
        os << getTime() << ' ' << particleHandler.getLastObject()->getPosition().Z << std::endl;
    }
    
    // Also write the ene information to the screen
    void printTime() const override
    {
        writeEneTimeStep(std::cout);
    }

};

int main(int argc, char **argv)
{
    // create an instance of the Penetration class
    Penetration dpm;
    // command line arguments:
    dpm.setNumberOfOMPThreads(helpers::readFromCommandLine(argc, argv, "-omp", 1));
    // turn on additional output files for viewing/analysing the data
    dpm.test(helpers::readFromCommandLine(argc, argv, "-test"));
    // turn on additional output files for viewing/analysing the data
    dpm.writeOutput(helpers::readFromCommandLine(argc, argv, "-writeOutput"));
    // read how many particles should be read (25K, 50K or 100K)
    dpm.setSize(helpers::readFromCommandLine(argc, argv, "-size", std::string("25K")));
    // call the solve routine
    dpm.solve();
    return 0;
}
