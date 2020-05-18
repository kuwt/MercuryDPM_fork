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
#include "Boundaries/CubeInsertionBoundary.h"
#include "Boundaries/SubcriticalMaserBoundary.h"
#include "Species/LinearViscoelasticSpecies.h"


/*! \brief Test for the MaserBoundaryOldStyle: make a chute-like domain with a maser inflow boundary in the beginning.
 * \details A maser boundary is created at 0<x<5. Thus, all particles in 0<x<5 are inside the maser region,
 * and all particles in x>5 are in the maser outflow.
 * The difference between MaserBoundary and MaserBoundaryOldStyle is that in MaserBoundary, there is a periodic box and
 * a separate outflow domain, where the outflow domain does not influence the periodic box. On the other hand, in
 * MaserBoundaryOldStyle, the left side of the periodic box does not have any influence on the rest of the flow, but the
 * right side of the periodic box and the left side of the outflow domain interact.
 *
 * The MaserBoundary is most useful for fast (supercritical) flows, and for flows for which the flux across the boundary
 * needs to be controlled. The MaserBoundaryOldStyle is more useful for slow (subcritical) flows, as the MaserBoundary
 * might generate "pulse-waves" in those cases.
 *
 * Check out the visualisation to see the exact behaviour of this boundary.
 */
class SubcriticalMaserBoundarySelfTest : public Mercury3D
{
public:
    
    void setupInitialConditions() override {
        setName("SubcriticalMaserBoundarySelfTest");
        
        //set species properties: some standard values
        LinearViscoelasticSpecies species;
        species.setDensity(6.0 / constants::pi);
        species.setStiffness(2e5);
        species.setDissipation(25);
        speciesHandler.copyAndAddObject(species);
        
        //set time and file properties
        setTimeStep(species.getCollisionTime(1) / 50.0);
        setTimeMax(5.0);
        setSaveCount(2000);
        setParticlesWriteVTK(true);
        
        //set domain size
        setMin({0,0,0});
        setMax({5, 2, 50});
        
        //for testing purposes, just set the gravity to 0.
        setGravity({0,0,0});
        
        //Check if particle is copied correctly when moving
        SphericalParticle p0;
        p0.setSpecies(speciesHandler.getLastObject());
        p0.setPosition({1,0,0});
        p0.setVelocity({1,0,0});
        p0.setRadius(1);
        particleHandler.copyAndAddObject(p0);
        
        //maser particle moving to the left should just continue
        p0.setPosition({1,0,3});
        p0.setVelocity({-1,0,0});
        particleHandler.copyAndAddObject(p0);
    
        //outflow particle moving to the left should continue, without changing species
        p0.setPosition({6,0,6});
        p0.setVelocity({-1,0,0});
        particleHandler.copyAndAddObject(p0);
        
        //maser particle on the left side of the periodic box should not influence periodic particle on the right
        p0.setPosition({1,0,9});
        p0.setVelocity({-1,0,0});
        particleHandler.copyAndAddObject(p0);
        p0.setPosition({4,0,9});
        p0.setVelocity({0,0,0});
        particleHandler.copyAndAddObject(p0);
    
        //maser particle on the right side of the periodic box should influence periodic particle on the left
        p0.setPosition({1,0,12});
        p0.setVelocity({0,0,0});
        particleHandler.copyAndAddObject(p0);
        p0.setPosition({4,0,12});
        p0.setVelocity({1,0,0});
        particleHandler.copyAndAddObject(p0);
    
        //outflow particle should influence periodic particle on the right side of the periodic box
        p0.setPosition({7,0,15});
        p0.setVelocity({-1,0,0});
        particleHandler.copyAndAddObject(p0);
        p0.setPosition({4.5,0,15});
        p0.setVelocity({0,0,0});
        particleHandler.copyAndAddObject(p0);
    
        //periodic particle on the right should influence outflow particle
        p0.setPosition({5.5,0,18});
        p0.setVelocity({0,0,0});
        particleHandler.copyAndAddObject(p0);
        p0.setPosition({3,0,18});
        p0.setVelocity({1,0,0});
        particleHandler.copyAndAddObject(p0);
    
        //periodic particle on the left should not influence outflow particle
        p0.setPosition({5.5,0,21});
        p0.setVelocity({0,0,0});
        particleHandler.copyAndAddObject(p0);
        p0.setPosition({2,0,21});
        p0.setVelocity({-1,0,0});
        particleHandler.copyAndAddObject(p0);
    
        //set the maser boundary
        SubcriticalMaserBoundary* b0 = boundaryHandler.copyAndAddObject(SubcriticalMaserBoundary());
        b0->set(Vec3D(1.0, 0.0, 0.0), 0.0, 5.0);
        write(std::cout, false);
    }
    
};

int main(int argc UNUSED, char* argv[] UNUSED)
{
    SubcriticalMaserBoundarySelfTest maserSelfTest;
    maserSelfTest.solve();
}
