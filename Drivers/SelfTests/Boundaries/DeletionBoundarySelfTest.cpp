//Copyright (c) 2013-2023, The MercuryDPM Developers Team. All rights reserved.
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

#include <iostream>
#include <Species/LinearViscoelasticSpecies.h>
#include <Mercury3D.h>
#include <Boundaries/DeletionBoundary.h>

class DeletionBoundarySelfTest : public Mercury3D
{
public:

    void setupInitialConditions() override {
        setName("DeletionBoundarySelfTest");
        setGravity(Vec3D(0,0,-1));
        setDomain(Vec3D(-1,-1,0),Vec3D(1,1,5));
        setTimeStep(1e-3);
        setTimeMax(3.0);
        setSaveCount(100);
        fStatFile.setFileType(FileType::NO_FILE);
        restartFile.setFileType(FileType::NO_FILE);
        eneFile.setFileType(FileType::NO_FILE);

        LinearViscoelasticSpecies s;
        s.setDensity(6.0/constants::pi);
        s.setStiffness(2e5);
        auto species = speciesHandler.copyAndAddObject(s);

        SphericalParticle p(species);
        p.setRadius(.5);
        p.setPosition(Vec3D(0,0,4));
        particleHandler.copyAndAddObject(p);

        //All particles whose position (x,y,z) satisfies (0,0,-1).(x,y,z)<-1, or z<1, get deleted
        DeletionBoundary b;
        b.set(Vec3D(0,0,-1),-1);
        boundaryHandler.copyAndAddObject(b);
    }

    void printTime() const override
    {
        if (particleHandler.getSize()!=0) logger(INFO,"t=%, Z=%", getTime(), particleHandler.getLastObject()->getPosition().Z);
    }
};

int main(int argc , char *argv[] )
{
    logger(INFO,"DeletionBoundarySelfTest: One particle falls under gravity, gets deleted at z=1.");
    DeletionBoundarySelfTest problem;
    problem.solve(argc, argv);
}
