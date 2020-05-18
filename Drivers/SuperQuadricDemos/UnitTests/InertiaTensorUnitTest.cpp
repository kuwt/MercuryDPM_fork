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
class InertiaTensorTester : public Mercury3D
{
public:
    void test()
    {
        auto species = new LinearViscoelasticSpecies();
        species->setDensity(1);
        speciesHandler.addObject(species);
        SuperQuadricParticle p;
        p.setSpecies(species);
        //spheres and ellipsoids
        p.setAxesAndExponents(1, 1, 1, 1, 1);
        particleHandler.copyAndAddObject(p);
        MatrixSymmetric3D inertia;
        inertia.XX = 8 * constants::pi / 15;
        inertia.YY = 8 * constants::pi / 15;
        inertia.ZZ = 8 * constants::pi / 15;
        helpers::check(p.getInertia(), inertia, 1e-5,"inertia tensor of sphere radius 1");

        p.setAxesAndExponents(2, 2, 2, 1, 1);
        inertia.XX = 256 * constants::pi / 15;
        inertia.YY = 256 * constants::pi / 15;
        inertia.ZZ = 256 * constants::pi / 15;
        helpers::check(p.getInertia(), inertia, 1e-5,"inertia tensor of sphere radius 2");

        p.setAxesAndExponents(2, 1, 1, 1, 1);
        inertia.XX = 16 * constants::pi / 15;
        inertia.YY = 8 * constants::pi / 3;
        inertia.ZZ = 8 * constants::pi / 3;
        helpers::check(p.getInertia(), inertia, 1e-5,"inertia tensor of (2,1,1)-ellipsoid");

        p.setAxesAndExponents(1, 2, 3, 1, 1);
        inertia.XX = 104 * constants::pi / 5;
        inertia.YY = 16 * constants::pi;
        inertia.ZZ = 8 * constants::pi;
        helpers::check(p.getInertia(), inertia, 1e-5,"inertia tensor of (1,2,3)-ellipsoid");
        logger(INFO, "All tests pass.");
    }
    
};

int main()
{
    InertiaTensorTester test;
    test.test();
    return 0;
}
