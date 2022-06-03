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

// T_protectiveWall

/* This tutorial presents the simulation of multiple particles rolling through a slope and then a protective
 * wall stop the flowing.
 * For full documentation of this code, go to http://docs.mercurydpm.org/Alpha/d0/db0/BeginnerTutorials.html#T5
*/
//

//! [AT_PW:headers]
#include <Mercury3D.h>
#include <Species/LinearViscoelasticSpecies.h>
#include "Walls/InfiniteWall.h"
#include "Walls/IntersectionOfWalls.h"
#include "Boundaries/CubeInsertionBoundary.h"
#include "Boundaries/CubeDeletionBoundary.h"
//! [AT_PW:headers]

class T_protectiveWall : public Mercury3D
{
public:

    T_protectiveWall(int Nump, Mdouble pRadius, Mdouble height, Mdouble width,Mdouble Length, Mdouble slopeAngle)
    {
        //! [AT_PW:globalConditions]
        setNumParticles = Nump;
        setGlobalRadius = pRadius;
        setParticleHeight = height;
        setSlopeAngle = constants::pi / 180.0 * slopeAngle;

        setName("T_protectiveWall");
        setWallsWriteVTK(FileType::MULTIPLE_FILES);
        setParticlesWriteVTK(true);

        setSystemDimensions(3);
        setGravity(Vec3D(0.0, 0.0, -9.81));

        setSaveCount(10);
        setTimeStep(0.005 / 50.0); // (collision time)/50.0

        double Hipo = Length/cos(setSlopeAngle);

        setXMax(Length);
        setYMax(width);
        setZMax(Hipo*sin(setSlopeAngle));
        //! [AT_PW:globalConditions]

        //! [AT_PW:slopeSpecies]
        slopeSpecies = speciesHandler.copyAndAddObject(LinearViscoelasticSpecies());
        slopeSpecies->setDensity(2500.0); //sets the species type_0 density
        slopeSpecies->setStiffness(20058.5);//sets the spring stiffness.
        slopeSpecies->setDissipation(0.01); //sets the dissipation.
        //! [AT_PW:slopeSpecies]

        //! [AT_PW:particleSpecies]
        particleSpecies = speciesHandler.copyAndAddObject(LinearViscoelasticSpecies());
        particleSpecies->setDensity(2500.0); //sets the species type_0 density
        particleSpecies->setStiffness(258.5);//sets the spring stiffness.
        particleSpecies->setDissipation(0.5); //sets the dissipation.
        //! [AT_PW:particleSpecies]

        speciesHandler.getMixedObject(slopeSpecies,particleSpecies)->mixAll(slopeSpecies,particleSpecies); //particle-wall interactions

        //! [AT_PW:particles]
        insb = new CubeInsertionBoundary;
        SphericalParticle particles;
        particles.setSpecies(particleSpecies);

        insb->set(
                &particles,
                1,
                Vec3D(0.9*getXMax(),  0.45*getYMax(), getZMax() + 0.97*setParticleHeight),
                Vec3D(getXMax(), 0.55*getYMax(), getZMax() + setParticleHeight),
                Vec3D(0, 0, 0),
                Vec3D(0, 0, 0)
        );
        PSD psd;
        psd.setDistributionUniform(pRadius*0.98, pRadius*1.1, 50);
        insb->setPSD(psd);
        insb = boundaryHandler.copyAndAddObject(insb);
        //! [AT_PW:particles]

        delb = new DeletionBoundary;
        delb->set(Vec3D(-1,0,0), getXMin());
        delb =  boundaryHandler.copyAndAddObject(delb);
    }

    LinearViscoelasticSpecies* slopeSpecies;
    LinearViscoelasticSpecies* particleSpecies;

    void setupInitialConditions() override
    {
        //! [AT_PW:walls]
        InfiniteWall slope,lateralwall1, lateralwall2, lateralwall3;

        //! [AT_PW:slope]
        slope.setSpecies(slopeSpecies);
        slope.setPosition(Vec3D(getXMin(),getYMin(),getZMin()));
        slope.setOrientation(Vec3D(sin(setSlopeAngle),0,-cos(setSlopeAngle)));
        wallHandler.copyAndAddObject(slope);
        //! [AT_PW:slope]

        //! [AT_PW:lateralWall]
        lateralwall1.setSpecies(slopeSpecies);
        lateralwall1.set(Vec3D(0.0, 1.0, 0.0), Vec3D(0.0,getYMax(),0.0));
        wallHandler.copyAndAddObject(lateralwall1);

        lateralwall2.setSpecies(slopeSpecies);
        lateralwall2.set(Vec3D(0.0, -1.0, 0.0), Vec3D(0.0,getYMin(),0.0));
        wallHandler.copyAndAddObject(lateralwall2);

        IntersectionOfWalls roarWall, proctWall;
        roarWall.addObject(Vec3D(1.0, 0.0, 0.0), Vec3D(getXMax(), 0.0, 0.0));
        roarWall.addObject(Vec3D(0.0, 0.0, -1.0), Vec3D(getXMax(), 0.0, getZMax()+setParticleHeight));
        wallHandler.copyAndAddObject(roarWall);

        //! [AT_PW:lateralWall]

        //! [AT_PW:protectiveWall]
        heightProtWall = setGlobalRadius*10.0;
        proctWall.setSpecies(slopeSpecies);
        proctWall.addObject(Vec3D(-1.0, 0.0, 0.0), Vec3D(getXMin(), 0.0, 0.0));
        proctWall.addObject(Vec3D(0.0, 0.0, -1.0), Vec3D(getXMin(), 0.0, heightProtWall));
        wallHandler.copyAndAddObject(proctWall);
        //! [AT_PW:protectiveWall]

        //! [AT_PW:walls]
    }

    void actionsAfterTimeStep() override {

        num_inserted = insb->getNumberOfParticlesInserted();
        vol_inserted = insb->getVolumeOfParticlesInserted();
        partDel = delb->getNumberOfParticlesDeleted();

        if(num_inserted == setNumParticles && !removed_insb){
            boundaryHandler.removeObject(insb->getIndex());
            removed_insb = true;
        }
        wallForce = fabs(wallHandler.getObject(4)->getForce().X);
        wallPressure = wallForce / (getYMax()*heightProtWall);
    }

    bool continueSolve() const override
    {
        static unsigned int counter = 0;
        if (++counter>100)
        {
            counter=0;
            if (getKineticEnergy()<0.01*getElasticEnergy())
                return false;
        }
        return true;
    }

    //! [AT_PW:printFunction]
    void printTime()const override
    {
        std::cout << "t=" << std::setprecision(3) << std::left<< getTime()
        << ", tmax=" << std::setprecision(3) << std::left<< getTimeMax()
        << ", # Particles inserted=" << std::setprecision(3) << std::left << setNumParticles - partDel
        << ", # Particles deleted=" << std::setprecision(3) << std::left << partDel
        << ", Volume inserted=" << std::setprecision(3) << std::left << std::setw(6)<< vol_inserted
        << ", WallForce=" << std::setprecision(3) << std::left << std::setw(6) << wallForce
        << ", WallPressure=" << std::setprecision(3) << std::left << std::setw(6) << wallPressure
//        << ", Kin=" << std::setprecision(3) << std::left << getKineticEnergy()
//        << ", Elast=" << std::setprecision(3) << std::left << getElasticEnergy()
        << std::endl;
    }
    //! [AT_PW:printFunction]

private:
    Mdouble setGlobalRadius = 1.0; // This radius defines the wall depth.
    Mdouble setParticleHeight = 0.1; // This height defines the particle insertion height, by default 0.1
    Mdouble setSlopeAngle = 20.0; // By defaults is 20 degrees
    Mdouble heightProtWall = 1.0;

    CubeInsertionBoundary* insb;
    DeletionBoundary* delb;
    int partDel = 0;
    int num_inserted =1;
    int setNumParticles;

    Mdouble wallForce = 0.0;
    Mdouble wallPressure = 0.0;
    Mdouble vol_inserted = 0.0;

    bool removed_insb = false;
};

int main(int argc, char* argv[])
{
    std::cout << "Write './T_protectiveWall -Np 100 -r 0.05 -h 0.5 -w 1.0 -l 5.0 -s 20.0 -t 2.0' to run the program" << std::endl;

    //! [AT_PW:setUp]
    int Nump = helpers::readFromCommandLine(argc,argv,"-Np",500); //100 partilces inserted by default
    Mdouble pRadius = helpers::readFromCommandLine(argc,argv,"-r",0.01); //by default is 2.0
    Mdouble height = helpers::readFromCommandLine(argc,argv,"-h",0.1); //by default is 2.0
    Mdouble width = helpers::readFromCommandLine(argc,argv,"-w",0.25);  //by default is 1.0
    Mdouble length = helpers::readFromCommandLine(argc,argv,"-l",1.0); //by default is 4.0
    Mdouble slopeAngle = helpers::readFromCommandLine(argc,argv,"-s",15.0); //by default is 20 degrees

    T_protectiveWall problem(Nump,pRadius,height,width,length,slopeAngle);

    Mdouble simTime = helpers::readFromCommandLine(argc,argv,"-t",20.0); //by default is 5.0
    problem.setTimeMax(simTime);
    //! [AT_PW:setUp]

    problem.removeOldFiles();
    problem.solve();

    return 0;
}

