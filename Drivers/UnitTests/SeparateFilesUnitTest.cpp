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

#include "DPMBase.h"
#include "Particles/SphericalParticle.h"
#include <iostream>
#include <Species/LinearViscoelasticSpecies.h>

class SeparateFilesSelfTest : public DPMBase{
public:

	void setupInitialConditions() override {
        setXMax(1.0);
        setYMax(1.0);
        setZMax(2.0);
        setSystemDimensions(3);
        setParticleDimensions(3);
        setGravity(Vec3D(0.0,0.0,0.0));

        LinearViscoelasticSpecies* species = speciesHandler.copyAndAddObject(LinearViscoelasticSpecies());
        species->setDensity(constants::pi/6);
        species->setStiffness(200000);

        particleHandler.clear();
        SphericalParticle p;
        p.setSpecies(speciesHandler.getObject(0));
        p.setPosition(Vec3D(0.5,0.5,0.5));
        p.setVelocity(Vec3D(0.0,0.0,0.0));
        p.setRadius(0.5);
        particleHandler.copyAndAddObject(p);
        SphericalParticle q;
        q.setSpecies(speciesHandler.getObject(0));
        q.setPosition(Vec3D(0.5,0.5,1.499));
        q.setVelocity(Vec3D(0.0,0.0,0.0));
        q.setRadius(0.5);
        particleHandler.copyAndAddObject(q);
        
        setTimeStep(1e-5);
		setTimeMax(6e-3);
        setSaveCount(helpers::getSaveCountFromNumberOfSavesAndTimeMaxAndTimeStep(2, getTimeMax(), getTimeStep()));
	}

};

int main(int argc UNUSED, char *argv[] UNUSED)
{
	SeparateFilesSelfTest problem;

	problem.setName("NoFiles");
    problem.setFileType(FileType::NO_FILE);
	problem.solve();
	
	problem.setName("CombinedFiles");
	problem.setFileType(FileType::ONE_FILE);
	problem.solve();
	
	problem.setName("SeparateFiles");
	problem.setFileType(FileType::MULTIPLE_FILES_PADDED);
	problem.solve();
}
