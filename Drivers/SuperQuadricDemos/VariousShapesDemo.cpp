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

#include "VTKWriter/SuperQuadricParticleVtkWriter.h"
#include "Species/HertzianViscoelasticMindlinSpecies.h"
#include "Walls/InfiniteWall.h"
#include "Mercury3D.h"
#include "Species/LinearViscoelasticSpecies.h"

/*!
 * Small class that shows the influence of varying the geometrical parameters
 * To test if your visualisation is working, you can first try the minimal example in VisualisationTest.
 */
class ShapesDemo : public Mercury3D
{
    void setupInitialConditions() override
    {
        setMin(-10, -10, 10);
        setMax(20, 20, 20);
        HertzianViscoelasticMindlinSpecies species;
        species.setEffectiveElasticModulusAndRestitutionCoefficient(20000, 0.8);
        species.setDensity(constants::pi / 6);
        speciesHandler.copyAndAddObject(species);
        
        SuperQuadricParticle p;
        //standard: sphere
        p.setSpecies(speciesHandler.getLastObject());
        p.setPosition({0,0,0});
        particleHandler.copyAndAddObject(p);
        
        //ellipsoid
        p.setAxes({1,1,2});
        p.setPosition({4, 0, 0});
        particleHandler.copyAndAddObject(p);
        
        //rounded beam
        p.setExponents(0.5, 0.5);
        p.setPosition({8, 0, 0});
        particleHandler.copyAndAddObject(p);
        
        //cylinder
        p.setAxes({1,1,1});
        p.setExponents(0.01, 1);
        p.setPosition({0, 0, -5});
        particleHandler.copyAndAddObject(p);
        
        //pillow
        p.setExponents(1, 0.01);
        p.setPosition({4,0,-5});
        particleHandler.copyAndAddObject(p);
        
        //cube
        p.setExponents(0.01, 0.01);
        p.setPosition({8, 0, -5});
        particleHandler.copyAndAddObject(p);
        
        
        setTimeStep(1e-5);
        setTimeMax(0);
        write(std::cout, true);
        SuperQuadricParticleVtkWriter vtkWriter(particleHandler);
        vtkWriter.writeVTK();
    }

private:
    
};

int main(int argc, char* argv[])
{
    ShapesDemo problem;
    problem.setName("ShapesDemo");
    problem.setSuperquadricParticlesWriteVTK(true);
    //problem.setWallsWriteVTK(FileType::ONE_FILE);
    problem.solve();
    return 0;
}

