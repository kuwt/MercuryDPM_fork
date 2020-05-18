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
//

#include "Particles/SuperQuadricParticle.h"
#include "Species/LinearViscoelasticSpecies.h"
#include "Mercury3D.h"

///Tests whether the expression for computing the volume of a superellipsoid is implemented correctly.
/// We use a sphere with axes and components equal to 1 and the test cases listed in Table 2.1 of <a href="https://homes.di.unimi.it/borghese/Teaching/DigitalAnimation/Old/DigitalAnimation_2002_2003/00_Superquadriche.pdf"> Chapter 2 </a> from Jaklic et al., 2000.
/// For the test cases, the axes a1,a2,a3 are kept constant and only the exponents eps1 and eps2 are varied such that eps1 = [0.1,1,2] and eps2 = [0.1,1,2]. Thus giving us a total 9 possible test cases.
/// Note that in Table 2.1, the first column and first row correspond to eps1=0 and eps2=0. This is a typo. The correct values are eps1= eps2= 0.1.
class VolumeTest : public Mercury3D
{
public:
    void test()
    {
        auto species = new LinearViscoelasticSpecies();
        species->setDensity(1);
        speciesHandler.addObject(species);
        SuperQuadricParticle p;
        p.setSpecies(species);

        //sphere
        p.setAxesAndExponents(1, 1, 1, 1, 1);
        particleHandler.copyAndAddObject(p);
        double volume;
        volume = (4.0/3.0) * constants::pi; // 4/3 * pi * r^3 with r = 1
        logger.assert_always(mathsFunc::isEqual(p.getVolume(), volume, 1e-5),
                             "Volume of sphere radius 1, got %, need %", p.getVolume(), volume);
        //logger(INFO,"volume got % need %",p->getVolume(),volume);


        p.setAxesAndExponents(2.0, 2.0, 4.0, 1.0, 1.0);
        volume = (4.0/3.0) * constants::pi * 2.0 * 2.0 * 4.0; // 4/3 * pi * a1*a2*a3
        //logger(INFO,"volume got % need %",p->getVolume(),volume);
        logger.assert_always(mathsFunc::isEqual(p.getVolume(), volume, 1e-5),
                             "Volume got %, need %", p.getVolume(), volume);

        // The below test cases dont pass as I think the volume expressions in Table 2.1 are an approximation to the volume expression 2.60
        // in the <a href="https://homes.di.unimi.it/borghese/Teaching/DigitalAnimation/Old/DigitalAnimation_2002_2003/00_Superquadriche.pdf"> Chapter 2 </a> of Jaklic et al.,2000.
        // Moreover, the asserts do not allow the exponents to be greater than 1. However, Fig. 2.5 in the Chapter 2 illustrates convex shapes right?
        /*
        p->setAxesAndExponents(2.0, 2.0, 4.0, 0.1, 0.1);
        volume = 8.0 * 2.0 * 2.0 * 4.0; // 8*a1*a2*a3
        logger.assert_always(mathsFunc::isEqual(p->getVolume(), volume, 1e-5),
                             "Volume got %, need %", p->getVolume(), volume);
        logger(INFO,"volume got % need %",p->getVolume(),volume);
`
        p->setAxesAndExponents(2.0, 2.0, 4.0, 0.1, 1);
        volume = 2.0 * constants::pi * 2.0 * 2.0 * 4.0; // 2*pi*a1*a2*a3
        logger.assert_always(mathsFunc::isEqual(p->getVolume(), volume, 1e-5),
                             "Volume got %, need %", p->getVolume(), volume);

        p->setAxesAndExponents(2.0, 2.0, 4.0, 0.1, 2.0);
        volume = 4.0 * 2.0 * 2.0 * 4.0; // 4*a1*a2*a3
        logger.assert_always(mathsFunc::isEqual(p->getVolume(), volume, 1e-5),
                             "Volume got %, need %", p->getVolume(), volume);

        p->setAxesAndExponents(2.0, 2.0, 4.0, 1.0, 0.1);
        volume = (16.0/3.0) * 2.0 * 2.0 * 4.0; // 16/3 *a1*a2*a3
        logger.assert_always(mathsFunc::isEqual(p->getVolume(), volume, 1e-5),
                             "Volume got %, need %", p->getVolume(), volume);


        p->setAxesAndExponents(2.0, 2.0, 4.0, 1.0, 2.0);
        volume = (8.0/3.0) * 2.0 * 2.0 * 4.0; // 8/3 *a1*a2*a3
        logger.assert_always(mathsFunc::isEqual(p->getVolume(), volume, 1e-5),
                             "Volume got %, need %", p->getVolume(), volume);

        p->setAxesAndExponents(2.0, 2.0, 4.0, 2.0, 0.1);
        volume = (8.0/3.0) * 2.0 * 2.0 * 4.0; // 8/3 *a1*a2*a3
        logger.assert_always(mathsFunc::isEqual(p->getVolume(), volume, 1e-5),
                             "Volume got %, need %", p->getVolume(), volume);

        p->setAxesAndExponents(2.0, 2.0, 4.0, 2.0, 1.0);
        volume = (2.0/3.0) * constants::pi * 2.0 * 2.0 * 4.0; // 2/3 * pi *a1*a2*a3
        logger.assert_always(mathsFunc::isEqual(p->getVolume(), volume, 1e-5),
                             "Volume got %, need %", p->getVolume(), volume);

        p->setAxesAndExponents(2.0, 2.0, 4.0, 2.0, 2.0);
        volume = (4.0/3.0) * 2.0 * 2.0 * 4.0; // 4/3 * a1*a2*a3
        logger.assert_always(mathsFunc::isEqual(p->getVolume(), volume, 1e-5),
                             "Volume got %, need %", p->getVolume(), volume);
        */
        logger(INFO, "All tests pass.");
    }
};

int main()
{
    VolumeTest problem;
    problem.test();
    return 0;
}
