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
#include "Boundaries/ConstantMassFlowMaserBoundary.h"
#include "Species/LinearViscoelasticSpecies.h"
#include "Walls/InfiniteWall.h"


/*! \brief Test for the MaserBoundary: make a chute-like domain with a maser inflow boundary in the beginning.
 * \details A maser boundary is created at 0<x<5. Thus, all particles in 0<x+gapSize<5 are inside the maser region,
 * and all particles in x>5 are in the maser outflow. The gapSize is 3 times the largest particle diameter.
 * All particles crossing the boundary of the maser region at x+gap=5 are duplicated in the maser outflow.
 *
 * To make this test different from the MaserSelfTest, we have now two different species, and change the mixed species
 * by hand. The mixed species after opening the maser should be the same as the ones set by the user.
 *
 * See MaserBoundary for more details.
 */
class ConstantMassFlowMaserBoundaryMixedSpeciesSelfTest : public Mercury3D
{
public:
    
    void printTime() const override {
        logger(INFO, "t = %, tmax = %", getTime(), getTimeMax());
    }
    
    void setupInitialConditions() override {
        setName("ConstantMassFlowMaserBoundaryMixedSpeciesSelfTest");
        
        //set species properties: a particle with diameter 1 should have mass 1
        LinearViscoelasticSpecies species;
        species.setDensity(6.0 / constants::pi);
        species.setStiffness(2e5);
        species.setDissipation(25);
        speciesHandler.copyAndAddObject(species);
        auto species2 = speciesHandler.copyAndAddObject(species);
        species2->setStiffness(1e5);
        
        auto mixedSpecies = dynamic_cast<LinearViscoelasticMixedSpecies*>(speciesHandler.getMixedObject(0,1));
        mixedSpecies->setStiffness(1.5e5);

        //set time and file properties
        setTimeStep(species.getCollisionTime(1) / 50.0);
        setTimeMax(4.0);
        setSaveCount(2000);
        
        //set domain properties, gravitational constant equals 1
        const double angle = 24.0 / 180.0 * constants::pi;
        setGravity(Vec3D(sin(angle), 0.0, -cos(angle)));
        
        setXMin(0.0);
        setYMin(0.0);
        setZMin(0.0);
        setXMax(12.0);
        setYMax(2.0);
        setZMax(2.0);
    
        SphericalParticle p;
        p.setSpecies(speciesHandler.getObject(0));
        p.setPosition(Vec3D(2.0,0.0,1.0));
        p.setVelocity(Vec3D(1.0,0.0,0.0));
        p.setRadius(0.5);
        particleHandler.copyAndAddObject(p);
        p.setPosition(Vec3D(4.0,0.0,1.0));
        p.setSpecies(speciesHandler.getObject(1));
        particleHandler.copyAndAddObject(p);
        
        //set the bottom of the chute: a solid (smooth) wall
        InfiniteWall w0;
        w0.setSpecies(speciesHandler.getObject(0));
        w0.set(Vec3D(0.0, 0.0, -1.0), Vec3D(0.0, 0.0, 0.0));
        wallHandler.copyAndAddObject(w0);
        
        //set the maser boundary
        ConstantMassFlowMaserBoundary* b0 = boundaryHandler.copyAndAddObject(ConstantMassFlowMaserBoundary());
        b0->set(Vec3D(1.0, 0.0, 0.0), 0.0, 5.0);
        write(std::cout, false);
    }
    
};

int main(int argc UNUSED, char* argv[] UNUSED)
{
    ConstantMassFlowMaserBoundaryMixedSpeciesSelfTest maserSelfTest;
    maserSelfTest.solve();
}


