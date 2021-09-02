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

#include "Mercury3D.h"
#include "Species/HertzianViscoelasticSlidingFrictionParhamiMcMeekingSinterSpecies.h"
#include "Walls/InfiniteWall.h"
#include "Boundaries/PeriodicBoundary.h"
#include "Walls/AxisymmetricIntersectionOfWalls.h"


int main(int argc UNUSED, char *argv[] UNUSED) {
    //create new simulation; as no specialised functions are needed, 
    //Mercury3D is used directly instead of creating a derived class as usual. 
    Mercury3D ps;

    //add Species
    HertzianViscoelasticSlidingFrictionParhamiMcMeekingSinterSpecies s;
    s.setDensity(3950e18);
    s.setEffectiveElasticModulus(0); //?? This seems to correspond to dt~tc/50
    s.setDissipation(0); //?
    // mass 2e-21 => fg = 2e-20
    // k=1e-04 => del=fg/k = tiny
    Mdouble alpha = 4.5;
    Mdouble beta = 4;
    Mdouble atomicVolume = 8.47e-30; /*Omega*/
    Mdouble surfaceEnergy = 1.1; /*gamma_s*/
    Mdouble thicknessDiffusion = 1.3e-8; /*deltaB*D0B*/
    Mdouble activationEnergy = 475e3 /*QB*/;
    Mdouble temperature = 1473; /*T*/
    Mdouble pseudoSlidingFrictionCoefficient = 0.01; /*\etaPart*/
    s.set(alpha, beta, atomicVolume, surfaceEnergy, thicknessDiffusion, activationEnergy, temperature, pseudoSlidingFrictionCoefficient);
    ps.speciesHandler.copyAndAddObject(s);

    //add particle
    SphericalParticle p;
    p.setSpecies(ps.speciesHandler.getObject(0));
    p.setRadius(50e-9);
    p.setPosition(p.getRadius() * Vec3D(0, 0, 1.0-1e-6*1e-3));
    ps.particleHandler.copyAndAddObject(p);

    //add second particle
    p.setPosition(p.getRadius() * Vec3D(0, 0, -1.0));
    ps.particleHandler.copyAndAddObject(p);

    //set time-stepping and output parameters
    ps.setFileType(FileType::ONE_FILE);
    ps.setXBallsAdditionalArguments(" -v0 -solidf ");
    ps.setGravity(Vec3D(0,0,0));
    ps.setTimeStep(1e-4);
    ps.setSaveCount(1000);
    ps.setTimeMax(10);
    ps.setMin(50e-9*Vec3D(-1,-1,-2));
    ps.setMax(50e-9*Vec3D(1,1,2));
    ps.setName("ParhamiMcKeepingParticlePairSelfTest");

    //run
    ps.solve();

    std::cout << "Execute 'gnuplot ParhamiMcKeepingParticlePairSelfTest.gnu' to view output" << std::endl;
    helpers::writeToFile("ParhamiMcKeepingParticlePairSelfTest.gnu",
                         "set xlabel 't [s]'\n"
                         "set ylabel 'delta/d'\n"
                         "p 'ParhamiMcKeepingParticlePairSelfTest.fstat' u 1:($7/1e-7) every 5::3 w lp\n"
                         );


}
