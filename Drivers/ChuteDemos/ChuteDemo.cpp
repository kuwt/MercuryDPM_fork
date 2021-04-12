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

//! [ChuteDemo:include]
#include <iostream>
#include <Species/LinearViscoelasticSpecies.h>
#include "Chute.h"
//! [ChuteDemo:include]

// Creates a quasi-2D inclined plane with inflow conditions on the left boundary, 
// and deletion of particles when they exit the domain on the right.
int main() 
{
    
    //! [ChuteDemo:initial]
    // Problem parameters
    Chute problem;
    
    problem.setName("ChuteDemo");   // data output file name
    problem.setGravity({0,-1,0});
    problem.setSaveCount(102);      // number of time steps skipped between saves
    Mdouble tc = 2.5e-3;            // collision time
    problem.setTimeStep(0.02 * tc); // actual time step
    problem.setTimeMax(0.5);        // maximum time
                                    // NB: number of time steps saved to output files
                                    // is timeMax/(timeStep*saveCount)
    //! [ChuteDemo:initial]
    
    //! [ChuteDemo:particles]
    // Particle radii
    problem.setFixedParticleRadius(0.001);  // radius of fixed chute bottom particles
    problem.setInflowParticleRadius(0.001); // radius of (monodisperse) inflow particles

    // Particle species
    LinearViscoelasticSpecies species;              // initialise species
    species.setHandler(&problem.speciesHandler);    // assign problem species handler to species
    species.setDensity(2000);                       // particle density
    species.setCollisionTimeAndRestitutionCoefficient( tc, 
            0.8, species.getMassFromRadius(problem.getInflowParticleRadius())); // material properties
    problem.speciesHandler.copyAndAddObject(species);   // assign species to problem species handler
    
    //! [ChuteDemo:particles]

    //! [ChuteDemo:chute]
    // Chute properties
    problem.setChuteAngle(30.0);                    // set angle of chute relative to horizontal
    problem.setXMax(0.1);                           // chute length = 0.1
    problem.setYMax(2.0 * problem.getInflowParticleRadius()); // chute width = 1 particle diameter
    //! [ChuteDemo:chute]

    //! [ChuteDemo:inflow]
    // Inflow properties
    problem.setInflowHeight(0.1);                   // particle inflow between 0 <= Z <= 0.1
    problem.setInflowVelocity(0.1);                 // particle inflow mean velocity
    problem.setInflowVelocityVariance(0.02);        // particle inflow velocity variance (in ratio of the mean velocity)
    //Write paraview data
    //problem.setParticlesWriteVTK(true);
    //problem.setWallsWriteVTK(true);
    //! [ChuteDemo:inflow]


    /*problem.setParticlesWriteVTK(true);
    problem.setWallsWriteVTK(true);*/


    //solve
    problem.solve();
} // the end
