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

// Tutorial 7

/*
** This file is annotated with DoxyFile comments in order to show the code on
** the documentation - This is not needed for your real drivers.
** Please ignore these comments.
**
** For full documentation of this code, go to http://docs.mercurydpm.org/Alpha/d0/db0/BeginnerTutorials.html#T7
*/

//! [T7:headers]
#include <Species/LinearViscoelasticSpecies.h>
#include <Mercury3D.h>
#include <Walls/InfiniteWall.h>
//! [T7:headers]

//! [T7:class]
class Tutorial7 : public Mercury3D
{
public:
    
    void setupInitialConditions() override {
        SphericalParticle p0;
        p0.setSpecies(speciesHandler.getObject(0));
        p0.setRadius(0.005);
        p0.setPosition(Vec3D(0.5 * getXMax(), 0.5 * getYMax(), 0.0));
        p0.setVelocity(Vec3D(1.2, 1.3, 0.0));
        particleHandler.copyAndAddObject(p0);
        
        //! [T7:infiniteWalls]
        InfiniteWall w0;
        
        w0.setSpecies(speciesHandler.getObject(0));
        w0.set(Vec3D(1.0, 0.0, 0.0), Vec3D(getXMax(), 0.0, 0.0));
        wallHandler.copyAndAddObject(w0);
        
        w0.set(Vec3D(-1.0, 0.0, 0.0), Vec3D(getXMin(), 0.0, 0.0));
        wallHandler.copyAndAddObject(w0);
        
        w0.set(Vec3D(0.0, 1.0, 0.0), Vec3D(0.0, getYMax(), 0.0));
        wallHandler.copyAndAddObject(w0);
        
        w0.set(Vec3D(0.0, -1.0, 0.0), Vec3D(0.0, getYMin(), 0.0));
        wallHandler.copyAndAddObject(w0);
        //! [T7:infiniteWalls]
    }
    
};
//! [T7:class]

//! [T7:main]
int main(int argc, char* argv[])
{
    // Problem setup
    Tutorial7 problem; // instantiate an object of class Tutorial 7
    
    problem.setName("Tutorial7");
    problem.setSystemDimensions(2);
    problem.setGravity(Vec3D(0.0, 0.0, 0.0));
    problem.setXMax(1);
    problem.setYMax(0.5);
    problem.setZMax(1.0);
    problem.setTimeMax(2.5);

    //! [T7:speciesProp]
    // The normal spring stiffness and normal dissipation is computed and set as
    // For collision time tc=0.005 and restitution coefficient rc=1.0,
    LinearViscoelasticSpecies species;
    species.setDensity(2500.0); //sets the species type_0 density
    species.setStiffness(258.5);//sets the spring stiffness.
    species.setDissipation(0.0); //sets the dissipation.
    problem.speciesHandler.copyAndAddObject(species);
    //! [T7:speciesProp]
    
    problem.setSaveCount(10);
    problem.dataFile.setFileType(FileType::ONE_FILE);
    problem.restartFile.setFileType(FileType::ONE_FILE);
    problem.fStatFile.setFileType(FileType::NO_FILE);
    problem.eneFile.setFileType(FileType::NO_FILE);

    problem.wallHandler.setWriteVTK(FileType::ONE_FILE);
    problem.setParticlesWriteVTK(true);
    
    problem.setXBallsAdditionalArguments("-solidf -v0 -s .85");
    
    problem.setTimeStep(0.005 / 50.0);// (collision time)/50.0
    problem.solve(argc, argv);
    
    return 0;
}
//! [T7:main]
