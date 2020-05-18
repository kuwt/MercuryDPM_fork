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

///tests whether the radius of the bounding sphere for superquadrics is computed correctly
class ContactDetectionTester : public Mercury3D
{
public:
    void test()
    {
        testSpheresContact();
        testEllipsoidsContact();
    
        logger(INFO, "All tests pass.");
    }
    


private:
    void setupParticles()
    {
        LinearViscoelasticSpecies species;
        species.setStiffness(2e5);
        species.setDissipation(25);
        speciesHandler.copyAndAddObject(species);
        p0 = new SuperQuadricParticle;
        p0->setSpecies(speciesHandler.getObject(0));
        p0->setAxesAndExponents(1.0, 1.0, 1.0, 1.0, 1.0);
        p0->setInertia();
        p1 = p0->copy();
        
        p0->setPosition(Vec3D(0.0, 0.0, 0.0));
        p1->setPosition(Vec3D(5.0, 0.0, 0.0));
        p0->setVelocity(Vec3D(0.0, 0.0, 0.0));
        p1->setVelocity(Vec3D(0.0, 0.0, 0.0));
        p0->setOrientationViaNormal({1.0, 0.0, 0.0});
        p1->setOrientationViaNormal({1.0, 0.0, 0.0});
        particleHandler.addObject(p0);
        particleHandler.addObject(p1);
    }
    
    ///\todo For each contact, check position, overlap and normal
    void testSpheresContact()
    {
        setupParticles();
        auto C = p0->getInteractionWith(p1, 0, &interactionHandler);
        logger.assert_always(C==nullptr, "Spheres far away should not touch");
        p1->setPosition({1.99, 0, 0});
        C = p0->getInteractionWith(p1, 0, &interactionHandler);
        logger.assert_always(C!=nullptr, "Spheres close together should touch");
        p1->setOrientationViaNormal({0, 1.0, 0});
        C = p0->getInteractionWith(p1, 0, &interactionHandler);
        logger.assert_always(C!=nullptr, "Rotated spheres close together should touch");
        p1->setPosition({3.0, 0, 0});
        C = p0->getInteractionWith(p1, 0, &interactionHandler);
        logger.assert_always(C==nullptr, "Rotated spheres far away should not touch");
        p1->setVelocity({-3.0, 0, 0});
        C = p0->getInteractionWith(p1, 0, &interactionHandler);
        logger.assert_always(C==nullptr, "Moving spheres far away should not touch");
        p1->setPosition({1.99, 0, 0});
        C = p0->getInteractionWith(p1, 0, &interactionHandler);
        logger.assert_always(C!=nullptr, "Moving spheres close together should touch");
        cleanup();
    }
    
    ///\todo For each contact, check position, overlap and normal
    void testEllipsoidsContact()
    {
        setupParticles();
        p0->setAxes(2, 1, 1);
        p1->setAxes(2, 1, 1);
    
        auto C = p0->getInteractionWith(p1, 0, &interactionHandler);
        logger.assert_always(C==nullptr, "Ellipsoids far away should not touch");
        logger.assert_always(!p0->isInContactWith(p1), "isInContactWith: Ellipsoids far away");
        p1->setPosition({3.99, 0, 0});
        C = p0->getInteractionWith(p1, 0, &interactionHandler);
        logger.assert_always(C!=nullptr, "Ellipsoids close together should touch");
        logger.assert_always(p0->isInContactWith(p1), " isInContactWith: Ellipsoids close together should touch");
        p1->setOrientationViaNormal({0,1,0});
        C = p0->getInteractionWith(p1, 0, &interactionHandler);
        logger.assert_always(C==nullptr, "Rotated ellipsoid with normal ellipsoid far away should not touch");
        p1->setPosition({2.99, 0, 0});
        C = p0->getInteractionWith(p1, 0, &interactionHandler);
        logger.assert_always(C!=nullptr, "One rotated ellipsoid with normal ellipsoid close together should touch");
        p0->setOrientationViaNormal({0, 0, 1});
        C = p0->getInteractionWith(p1, 0, &interactionHandler);
        logger.assert_always(C==nullptr, "Rotated ellipsoids  far away should not touch");
        p1->setPosition({1.99, 0, 0});
        C = p0->getInteractionWith(p1, 0, &interactionHandler);
        logger.assert_always(C!=nullptr, "Rotated ellipsoids close together should touch");
        p0->setOrientationViaNormal({1,0,0});
        p1->setPosition({1.99, -1.99, 0});
        C = p0->getInteractionWith(p1, 0, &interactionHandler);
        logger.assert_always(C!=nullptr, "One rotated ellipsoid with normal ellipsoid close together should touch");
        cleanup();
    }
    
    void cleanup()
    {
        particleHandler.removeObject(1);
        particleHandler.removeObject(0);
    }
    
    SuperQuadricParticle* p0;
    SuperQuadricParticle* p1;
    
};

int main()
{
    ContactDetectionTester test;
    test.test();
    return 0;
}
