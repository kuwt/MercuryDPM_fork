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
#include "Math/ExtendedMath.h"
#include "Species/LinearViscoelasticSpecies.h"
#include "Boundaries/PeriodicBoundary.h"
#include "CG/CG.h"
//number of particles in all directions
const unsigned n = 5;

class Packing : public DPMBase
{
public:
    void setupInitialConditions() override
    {
        //stiffness and distance are set such that all forces are 1
        Mdouble distance = 0.99;

        LinearViscoelasticSpecies species;
        species.setStiffnessAndRestitutionCoefficient(100.0,0.5,1.0);
        species.setDensity(6.0 / constants::pi);
        speciesHandler.copyAndAddObject(species);

        //put particles into the domain in a 3d grid of mesh size distance
        SphericalParticle p;
        p.setSpecies(speciesHandler.getObject(0));
        p.setRadius(0.5);
        for (unsigned i=0; i<n; i++) {
            for (unsigned j = 0; j < n; j++) {
                for (unsigned k = 0; k < n; k++) {
                    p.setPosition(distance * Vec3D(i + .5, j + .5, k + .5));
                    particleHandler.copyAndAddObject(p);
                }
            }
        }

        // define the domain
        setMin(0, 0, 0);
        setMax(n*distance, n*distance, n*distance);

        //create periodic boundaries around the domain
        PeriodicBoundary b;
        b.set(Vec3D(1,0,0),getMin(),getMax());
        boundaryHandler.copyAndAddObject(b);
        b.set(Vec3D(0,1,0),getMin(),getMax());
        boundaryHandler.copyAndAddObject(b);
        b.set(Vec3D(0,0,1),getMin(),getMax());
        boundaryHandler.copyAndAddObject(b);

        dataFile.setFileType(FileType::NO_FILE);
        fStatFile.setFileType(FileType::NO_FILE);
        setSaveCount(20);
        setTimeStep(0.02*species.getCollisionTime(1.0));
        setTimeMax(10*getTimeStep());
    }

    void test()
    {
        setName("CG3DPackingUnitTest");
        auto cg = cgHandler.copyAndAddObject(CG<CGCoordinates::O>());
        cg->statFile.setSaveCount(100);
        solve();
        const Mdouble density = cg->getPoint(0).getDensity();
        const Matrix3D contactStress = cg->getPoint(0).getContactStress();
        const Vec3D domainSize = getMax()-getMin();
        const Mdouble volume = domainSize.X*domainSize.Y*domainSize.Z;
        const Mdouble realDensity = n*n*n/volume;
        const Mdouble realPressure = n*n*n*0.99/volume;

        logger.assert_always(mathsFunc::isEqual(density,realDensity,1e-10),
                             "Wrong density: % (should be %)",density,realDensity);
        logger.assert_always(mathsFunc::isEqual(contactStress.ZZ,realPressure,1e-10),
                             "Wrong contactStress.XX: % (should be %)",contactStress.XX,realPressure);
        logger.assert_always(mathsFunc::isEqual(contactStress.ZZ,realPressure,1e-10),
                             "Wrong contactStress.YY: % (should be %)",contactStress.YY,realPressure);
        logger.assert_always(mathsFunc::isEqual(contactStress.ZZ,realPressure,1e-10),
                             "Wrong contactStress.ZZ: % (should be %)",contactStress.ZZ,realPressure);
        logger(INFO,"Test successful");
    }

};



int main()
{
    logger(INFO,"Create cubic packing and checking the resulting forces and density.");
    Packing problem;
    problem.test();
    return 0;
}
