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
/// At the moment it only checks for ellipsoids and spheres using simplified expressions in the SuperQuadricParticle::getMaxInteractionRadius() function.
/// Definitely needs generalisation for complex shapes.

class BoundingRadiusTester : public Mercury3D
{
public:
    void test()
    {
        SuperQuadricParticle p;
        auto species = new LinearViscoelasticSpecies();
        speciesHandler.addObject(species);
        p.setSpecies(species);
        
        //spheres and ellipsoids
        
        p.setAxesAndExponents(1, 1, 1, 1, 1);
        logger.assert_always(mathsFunc::isEqual(p.getMaxInteractionRadius(), 1, 1e-5),
                             "interaction radius of sphere");
        
        p.setAxesAndExponents(2, 2, 2, 1, 1);
        logger.assert_always(mathsFunc::isEqual(p.getMaxInteractionRadius(), 2, 1e-5),
                             "interaction radius of sphere");
        
        //ellipsoids
        p.setAxesAndExponents(2, 1, 1, 1, 1);
        logger.assert_always(mathsFunc::isEqual(p.getMaxInteractionRadius(), 2, 1e-5),
                             "interaction radius of ellipsoid");
        
        p.setAxesAndExponents(0.5, .2, .1, 1, 1);
        logger.assert_always(mathsFunc::isEqual(p.getMaxInteractionRadius(), 0.5, 1e-5),
                             "interaction radius of ellipsoid");
        
        //same axes but other epsilon1, epsilon2
        p.setAxesAndExponents(1, 1, 1, 1, 0.5);
        logger.assert_always(mathsFunc::isEqual(p.getMaxInteractionRadius(), 1.1892, 1e-2),
                             "interaction radius of epsilon2=0.5, equal axes. "
                                     "Expected % got %", 1.1892, p.getMaxInteractionRadius());
        //
        p.setAxesAndExponents(3, 2, 1, 1, 0.125);
        logger.assert_always(mathsFunc::isEqual(p.getMaxInteractionRadius(), 3.4711, 1e-2),
                             "interaction radius of epsilon2=0.125, unequal axes. "
                                     "Expected % got %", 3.4711, p.getMaxInteractionRadius());

        
        
        p.setAxesAndExponents(1, 1, 1, 0.5, 1);
        logger.assert_always(mathsFunc::isEqual(p.getMaxInteractionRadius(), 1.1892, 1e-2),
                             "interaction radius of epsilon1=0.5, equal axes");
        
        
        p.setAxesAndExponents(1, 1, 1, 0.5, 0.5);
        logger.assert_always(mathsFunc::isEqual(p.getMaxInteractionRadius(), 1.3161, 1e-2),
                             "interaction radius of epsilon1 = 0.5, epsilon2=0.5, equal axes."
                                     " Expected % got %", 1.3161, p.getMaxInteractionRadius());
        
        
        p.setAxesAndExponents(3, 2, 1, 0.5, 0.5);
        logger.assert_always(mathsFunc::isEqual(p.getMaxInteractionRadius(), 3.1461, 1e-2),
                             "interaction radius of epsilon1 = 0.5, epsilon2=0.5, unequal axes."
                                     " Expected % got %", 3.1461, p.getMaxInteractionRadius());
        
        
        p.setAxesAndExponents(1, 1, 1, 0.25, 0.5);
        logger.assert_always(mathsFunc::isEqual(p.getMaxInteractionRadius(), 1.4283, 1e-2),
                             "interaction radius of epsilon1 = 0.25, epsilon2=0.5, equal axes."
                                     " Expected % got %", 1.4283, p.getMaxInteractionRadius());
        
        
        p.setAxesAndExponents(3, 2, 1, 0.25, 0.5);
        logger.assert_always(mathsFunc::isEqual(p.getMaxInteractionRadius(), 3.1928, 1e-2),
                             "interaction radius of epsilon1 = 0.25, epsilon2=0.5, unequal axes."
                                     " Expected % got %", 3.1928, p.getMaxInteractionRadius());
        
        logger(INFO, "All tests pass.");
    }
    
};

int main()
{
    BoundingRadiusTester test;
    test.test();
    return 0;
}
