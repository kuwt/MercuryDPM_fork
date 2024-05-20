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


//#include <Species/LinearViscoelasticFrictionSpecies.h>
#include <Mercury3D.h>
#include <Walls/InfiniteWall.h>
//#include "Particles/SuperQuadricParticle.h"
#include "Species/LinearViscoelasticFrictionJKRAdhesiveSpecies.h"

class SlidingAndRollingParticlesTEst : public Mercury3D
{
public:

    void setupInitialConditions() override {
/*        SphericalParticle p0;

        //sets the particle to species type-1
        p0.setSpecies(speciesHandler.getObject(0));
        p0.setRadius(20e-6);//);//0.005
        p0.setPosition(Vec3D(250e-6, 250e-6, getZMin() + p0.getRadius()));
        p0.setVelocity(Vec3D(0.0, 0.0, 0.0));
        particleHandler.copyAndAddObject(p0);

        //sets the particle to species type-2
        p0.setSpecies(speciesHandler.getObject(0));
        p0.setRadius(37e-6);//);//0.005*2
        p0.setPosition(Vec3D(250e-6, 500e-6, getZMin() + p0.getRadius()));
        p0.setVelocity(Vec3D(0.0, 0.0, 0.0));
        particleHandler.copyAndAddObject(p0);

        //sets the particle to species type-3
        p0.setSpecies(speciesHandler.getObject(0));
        p0.setRadius(54e-6);//);//0.005*4
        p0.setPosition(Vec3D(250e-6 , 750e-6, getZMin() + p0.getRadius()));
        p0.setVelocity(Vec3D(0.0, 0.0, 0.0));
        particleHandler.copyAndAddObject(p0);*/
        //
        //sets the particle to species type-3
       SuperQuadricParticle p1;
        p1.setAxesAndExponents(54e-6,54e-6,75e-6,1,1);
        p1.setPosition(Vec3D(250e-6 , 750e-6, 200e-6));//80e-6));
        p1.setVelocity(Vec3D(0.0, 0.0, 0.0));
        p1.setSpecies(speciesHandler.getObject(0));
        particleHandler.copyAndAddObject(p1);
        //p1.setRadius(0.005);
        //p1.setExponents(1, 1);
        //p1.setAxes(0.5, 0.5, 0.5);
        //p1.setAxes(1.8, std::sqrt(0.2), std::sqrt(0.2));

        //
        InfiniteWall w0;

        w0.setSpecies(speciesHandler.getObject(0));
        w0.set(Vec3D(0.0, 0.0, -1.0), Vec3D(0.0, 0.0, getZMin()));
        wallHandler.copyAndAddObject(w0);

    }
    //
    void writeFstatHeader(std::ostream& os) const override
    {
        for (std::vector<BaseInteraction*>::const_iterator it = interactionHandler.begin(); it != interactionHandler.end(); ++it)
        {
            //(*it)->writeToFStat(os);
            //
/*            << " " << (*it)->getOverlap()
                   << " " << (*it)->getContactPoint()
                   << " " << (*it)->getForce()
                   << " " << (*it)->getRelativeVelocity()*/
            //
            os << (*it)->getTimeStamp()
               << " " << (*it)->getP()->getIndex()
               << " " << (*it)->getI()->getIndex()
               << " " << (*it)->getTorque()
               << std::endl;
        }
        //
/*        for (const auto &q: particleHandler) {
            auto p0 = dynamic_cast<SphericalParticle *>(q);
            os << p0->getRadius()
               << " " << p0->getVelocity()
               << " " << p0->getAngularVelocity()
               << " " << p0->getTorque()
               << std::endl;
        }*/
    }
};


int main(int argc, char* argv[])
{

    // Problem setup
    SlidingAndRollingParticlesTEst problem;

    double angle = constants::pi / 180.0 * 45.0;

    problem.setName("SlidingAndRollingParticlesTest");
    problem.setSystemDimensions(3);
    problem.setGravity(Vec3D(0.0,0.0,0.0));//sin(angle), 0.0, -cos(angle)) * 9.81);
    problem.setXMax(10000e-6);//0.3
    problem.setYMax(1000e-6);//0.3
    problem.setZMax(250e-6);//0.05
    problem.setTimeMax(0.05);

    Mdouble MatDensity = 4430;
    Mdouble tc = 90e-7; //calculated as tc = 0.005-0.01 * tg : tg=sqrt(d50/g) => 0.005*sqrt(35e-6/9.81)=
    Mdouble timeStep = 0.02*tc;
    Mdouble restitutionCoeff = 0.4;//0.1;
    Mdouble minRaduis = 12.0e-6;///2.0;//MinRadius;//ALL PREVIOUS CASES THIS WAS USED -> 12.0e-6;   //15.0e-6/2.0;//13.183e-6/2.0;
    //
    auto species = problem.speciesHandler.copyAndAddObject(LinearViscoelasticFrictionJKRAdhesiveSpecies());

    species->setDensity(MatDensity);
    Mdouble mass = species->getMassFromRadius(minRaduis);
    species->setCollisionTimeAndRestitutionCoefficient(tc,restitutionCoeff,mass);
    //
    Mdouble adhesionStiffness = 0.5 * species->getStiffness();
    Mdouble surfaceEnergy = 0.0;//0.1e-3;
    //species->setAdhesionStiffness(adhesionStiffness);
    species->setSurfaceEnergy(surfaceEnergy);
    //
    species->setSlidingStiffness(2./7.*species->getStiffness());
    species->setSlidingDissipation(2./7.*species->getDissipation());
    species->setSlidingFrictionCoefficient(0.4);
    species->setRollingStiffness(2./5.*species->getStiffness());
    species->setRollingDissipation(2./5.*species->getDissipation());
    species->setRollingFrictionCoefficient(0.1);

    //
    //Properties index 0
    //LinearViscoelasticFrictionSpecies species0;
    //LinearViscoelasticFrictionJKRAdhesiveSpecies species0;


/*    species0.setDensity(2500.0);//sets the species type_0 density
    //species0.setStiffness(259.018);//sets the spring stiffness.
    //species0.setDissipation(0.0334);//sets the dissipation.
    species0.setSlidingStiffness(2.0 / 7.0 * species0.getStiffness());
    species0.setRollingStiffness(2.0 / 5.0 * species0.getStiffness());
    species0.setSlidingFrictionCoefficient(0.1);
    species0.setRollingFrictionCoefficient(0.1);
    //
    Mdouble adhesionStiffness = 0.5 * species0.getStiffness();
    Mdouble surfaceEnergy = 0.4e-3;
    species0.setAdhesionStiffness(adhesionStiffness);
    species0.setSurfaceEnergy(surfaceEnergy);
    //
    auto ptrToSp0=problem.speciesHandler.copyAndAddObject(species0);*/


    //Properties index 1
   /* LinearViscoelasticFrictionSpecies species1;

    species1.setDensity(2500.0);//sets the species type-1 density
    species1.setStiffness(259.018);//sets the spring stiffness
    species1.setDissipation(0.0334);//sets the dissipation
    species1.setSlidingStiffness(2.0 / 7.0 * species1.getStiffness());
    species1.setRollingStiffness(2.0 / 5.0 * species1.getStiffness());
    species1.setSlidingFrictionCoefficient(0.0);
    species1.setRollingFrictionCoefficient(0.1);
    auto ptrToSp1=problem.speciesHandler.copyAndAddObject(species1);

    //Combination of properties index 0 and index 1
    auto species01 = problem.speciesHandler.getMixedObject(ptrToSp0,ptrToSp1);

    species01->setStiffness(259.018);//sets the spring stiffness
    species01->setDissipation(0.0334);//sets the dissipation
    species01->setSlidingStiffness(2.0 / 7.0 * species01->getStiffness());
    species01->setRollingStiffness(2.0 / 5.0 * species01->getStiffness());
    species01->setSlidingFrictionCoefficient(0.0);
    species01->setRollingFrictionCoefficient(0.1);

    //Properties index 2
    LinearViscoelasticFrictionSpecies species2;

    species2.setDensity(2500.0);//sets the species type-2 density
    species2.setStiffness(258.5);//sets the spring stiffness
    species2.setDissipation(0.0);//sets the dissipation
    species2.setSlidingStiffness(2.0 / 7.0 * species2.getStiffness());
    species2.setRollingStiffness(2.0 / 5.0 * species2.getStiffness());
    species2.setSlidingFrictionCoefficient(0.0);
    species2.setRollingFrictionCoefficient(0.1);
    auto ptrToSp2 = problem.speciesHandler.copyAndAddObject(species2);

    //Combination of properties index 0 and index 2
    auto species02 = problem.speciesHandler.getMixedObject(ptrToSp0, ptrToSp2);

    species02->setStiffness(259.018);//sets the stiffness
    species02->setDissipation(0.0334);//sets the dissipation
    species02->setSlidingStiffness(2.0 / 7.0 * species02->getStiffness());
    species02->setRollingStiffness(2.0 / 5.0 * species02->getStiffness());
    species02->setSlidingFrictionCoefficient(0.0);
    species02->setRollingFrictionCoefficient(0.1);*/

    //Output
    problem.setSaveCount(1000);
    problem.dataFile.setFileType(FileType::ONE_FILE);
    problem.restartFile.setFileType(FileType::ONE_FILE);
    problem.fStatFile.setFileType(FileType::ONE_FILE);
    problem.eneFile.setFileType(FileType::NO_FILE);

    //problem.setXBallsAdditionalArguments("-solidf -v0 -s 8");

    problem.setParticlesWriteVTK(true);
    //problem.setSuperquadricParticlesWriteVTK(true);
    problem.setWallsWriteVTK(FileType::MULTIPLE_FILES);

    problem.setTimeStep(timeStep);//0.005 / 50.0);
    problem.solve(argc, argv);

    return 0;
}
