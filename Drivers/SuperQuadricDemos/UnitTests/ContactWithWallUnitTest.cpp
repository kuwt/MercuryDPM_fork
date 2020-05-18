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

#include "Particles/SuperQuadricParticle.h"
#include "Species/LinearViscoelasticSpecies.h"
#include "Mercury3D.h"
#include "Walls/InfiniteWall.h"

///tests whether the radius of the bounding sphere for superquadrics is computed correctly
class ContactDetectionWithWallTester : public Mercury3D
{
public:
    void test()
    {
        testSpheresContact();
        testEllipsoidsContact();
        
        logger(INFO, "All tests pass.");
    }



private:
    void setupParticleAndWall()
    {
        LinearViscoelasticSpecies species;
        species.setStiffness(2e5);
        species.setDissipation(25);
        speciesHandler.copyAndAddObject(species);
        p0 = new SuperQuadricParticle;
        p0->setSpecies(speciesHandler.getObject(0));
        p0->setAxesAndExponents(1.0, 1.0, 1.0, 1.0, 1.0);
        p0->setInertia();
        
        p0->setPosition(Vec3D(0.0, 0.0, 0.0));
        p0->setVelocity(Vec3D(0.0, 0.0, 0.0));
        p0->setOrientationViaNormal({1.0, 0.0, 0.0});
        particleHandler.copyAndAddObject(p0);
        
        w0 = new InfiniteWall;
        w0->setSpecies(speciesHandler.getObject(0));
        w0->set(Vec3D(0.0,0.0,-1.0),Vec3D(0.0,0.0,0.0));
        wallHandler.copyAndAddObject(w0);
    }
    
    ///\todo For each contact, check position, overlap and normal
    void testSpheresContact()
    {
        setupParticleAndWall();
        p0->setPosition({0, 0, 0.5});
        BaseInteraction* C = w0->getInteractionWith(p0, 0, &interactionHandler);
        logger.assert_always(C!=nullptr, "particle should be in interaction with wall");
        logger.assert_always(mathsFunc::isEqual(C->getOverlap(), 0.5, 1e-5), "overlap with sphere-wall contact");
        logger.assert_always(mathsFunc::isEqual(C->getContactPoint(), Vec3D(0,0,-.25), 1e-5), "contact point with sphere-wall contact");
        p0->setOrientationViaNormal({0,1,0});
        C = w0->getInteractionWith(p0, 0, &interactionHandler);
        logger.assert_always(C!=nullptr, "rotated particle should be in interaction with wall");
        logger.assert_always(mathsFunc::isEqual(C->getOverlap(), 0.5, 1e-5), "overlap with rotated sphere-wall contact");
        logger.assert_always(mathsFunc::isEqual(C->getContactPoint(), Vec3D(0,0,-.25), 1e-5), "contact point with rotated sphere-wall contact");
        p0->setPosition({0, 0, 1.5});
        C = w0->getInteractionWith(p0, 0, &interactionHandler);
        logger.assert_always(C==nullptr, "Sphere should not be in contact with wall");
        cleanup();
    }
    
    ///\todo For each contact, check position, overlap and normal
    void testEllipsoidsContact()
    {
        setupParticleAndWall();
        p0->setAxes(2, 1, 1);
        p0->setPosition({0, 0, 0.5});
        BaseInteraction* C = w0->getInteractionWith(p0, 0, &interactionHandler);
        logger.assert_always(C!=nullptr, "particle should be in interaction with wall");
        logger.assert_always(mathsFunc::isEqual(C->getOverlap(), 0.5, 1e-5), "overlap with ellipsoid-wall contact");
        logger.assert_always(mathsFunc::isEqual(C->getContactPoint(), Vec3D(0,0,-.25), 1e-5), "contact point with ellipsoid-wall contact");
        p0->setAxes(1,1,2);
        p0->setPosition({0, 0, 1.5});
        C = w0->getInteractionWith(p0, 0, &interactionHandler);
        logger.assert_always(C!=nullptr, "rotated ellipsoid should be in interaction with wall");
        logger.assert_always(mathsFunc::isEqual(C->getOverlap(), 0.5, 1e-5), "overlap with rotated ellipsoid-wall contact");
        logger.assert_always(mathsFunc::isEqual(C->getContactPoint(), Vec3D(0,0,-.25), 1e-5), "contact point with rotated ellipsoid-wall contact");
        p0->setOrientationViaNormal({0,0,1});
        C = w0->getInteractionWith(p0, 0, &interactionHandler);
        logger.assert_always(C==nullptr, "rotated ellipsoid should not be in contact with wall");
        cleanup();
    }
    
    void cleanup()
    {
        delete p0;
    }
    
    SuperQuadricParticle* p0;
    InfiniteWall* w0;
    
};

int main()
{
    ContactDetectionWithWallTester test;
    test.test();
    return 0;
}
