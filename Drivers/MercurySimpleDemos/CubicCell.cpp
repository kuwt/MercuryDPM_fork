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

#include <Species/Species.h>
#include <Species/LinearViscoelasticSpecies.h>
#include <Mercury3D.h>
#include <Particles/BaseParticle.h>
#include <Boundaries/CubeInsertionBoundary.h>
#include <Boundaries/PeriodicBoundary.h>

class CubicCell : public Mercury3D
{
public:

    void actionsAfterTimeStep() override
    {
        for (BaseParticle* p : particleHandler)
        {
            p->setRadius((1.0 + getTimeStep() * growthRate) * p->getRadius());
        }
        ///\todo IFCD: I uncommented the next lines because they prevent a lot of warnings. Please check if uncommenting is correct.
        if (hGridNeedsRebuilding())
        {
            hGridRebuild();
        }
    }

    void setupInitialConditions() override
    {
        LinearViscoelasticSpecies species; // defines linear-viscoelastic contact law
        species.setDensity(6.0 / constants::pi); // sets the particle density
        species.setStiffness(2e3); // sets the linear contact stiffness
        species.setDissipation(25.0); // sets the dissipation
        speciesHandler.copyAndAddObject(species);
        logger(INFO, "Restitution coefficient for a particle with species 0, mass 1: %",
               species.getRestitutionCoefficient(1.0));
            
        //add periodic walls
        PeriodicBoundary b0;
		b0.set(Vec3D(1,0,0), getXMin(), getXMax());
		boundaryHandler.copyAndAddObject(b0);
		b0.set(Vec3D(0,1,0), getYMin(), getYMax());
		boundaryHandler.copyAndAddObject(b0);
		b0.set(Vec3D(0,0,1), getZMin(), getZMax());
		boundaryHandler.copyAndAddObject(b0);

        //insert the particles in a cubic domain
        SphericalParticle p; // defines a new particle
        p.setSpecies(speciesHandler.getObject(0)); // sets particle species

        CubeInsertionBoundary b;
        b.set(&p, 100,
              Vec3D(getXMin(), getYMin(), getZMin()),
              Vec3D(getXMax(), getYMax(), getZMax()),
              Vec3D(0, 0, 0), Vec3D(0, 0, 0), .5, 1.0);
        ///\todo Periodic walls are still not checked in checkParticleForInteraction!
        
        BaseParticle* q;
        for (unsigned int i = 0; i < n * n * n; ++i)
        {
            q = b.generateParticle(random);
            if (checkParticleForInteraction(*q))
            {
                particleHandler.copyAndAddObject(q);
            }
        }

    }

    Mdouble getHGridTargetMaxInteractionRadius() const override 
    {
        return 1.1 * MercuryBase::getHGridTargetMaxInteractionRadius();
    };
    
    unsigned int n = 10;
    Mdouble growthRate = 4e-2;
};

int main(int argc, char* argv[])
{
    CubicCell problem;
    problem.setName("CubicCell");
    problem.setSystemDimensions(3);
    problem.setGravity(Vec3D(0.0, 0.0, 0.0));
    problem.setXMax(problem.n);
    problem.setYMax(problem.n);
    problem.setZMax(problem.n);
    problem.setTimeStep(1e-3);
    problem.setTimeMax(10.0);
    problem.fStatFile.setFileType(FileType::NO_FILE);
    problem.solve();
    return 0;
}
