//Copyright (c) 2013-2020, The MercuryDPM Developers Team. All rights reserved.
//For the list of developers, see <http://www.MercuryDPM.org/Team>.
//
//Redistribution and use in source and MarbleRun forms, with or without
//modification, are permitted provided that the following conditions are met:
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//  * Redistributions in MarbleRun form must reproduce the above copyright
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

#include<iostream>
#include <Species/LinearViscoelasticFrictionSpecies.h>
#include <Walls/TriangleWall.h>
#include "Mercury3D.h"
#include "Walls/InfiniteWall.h"

//#define DEBUG_OUTPUT

class MarbleRun : public Mercury3D 
{
    // pointer to the particle
    SphericalParticle* particle;
    // pointer to the species storing the particle/contact properties
    LinearViscoelasticFrictionSpecies* species;
    // restitution is stored separately, because it is not a direct contact property
    double restitutionCoefficient;

public:

    // create default particle and species (properties to be filled by the user)
    MarbleRun() {
        species = speciesHandler.copyAndAddObject(LinearViscoelasticFrictionSpecies());
        particle = particleHandler.copyAndAddObject(SphericalParticle(species));
    }

    void setParticlePosition(Vec3D pos) {
        particle->setPosition(pos);
    }

    void setParticleRadius(double radius) {
        particle->setRadius(radius);
    }

    void setParticleDensity(double density) {
        species->setDensity(density);
    }

    void setSlidingFrictionCoefficient(double coeff) {
        species->setSlidingFrictionCoefficient(coeff);
    }

    void setRollingFrictionCoefficient(double coeff) {
        species->setRollingFrictionCoefficient(coeff);
    }

    void setTorsionFrictionCoefficient(double coeff) {
        species->setTorsionFrictionCoefficient(coeff);
    }

    double getParticleMass() {
        return species->getMassFromRadius(particle->getRadius());
    }
    
    void setRestitutionCoefficient(double coeff) {
        restitutionCoefficient = coeff;
    }
    
    void loadSTLFile(std::string stlFile) {
        // read assuming the file is written in mm units
        wallHandler.readTriangleWall(stlFile,speciesHandler.getLastObject(),1e-3);
    }
    
    /// groupId of the big seesaw (so we can find all the triangles belonging to it)
    unsigned bigSeesawId;
    /// first element in the wallHandler belonging to the big seesaw
    /// (so we can extract information that is equal for all seesaw triangles,
    /// like getOrientation, getAngularVelocity)
    BaseWall* bigSeesawFirstElement = nullptr;
    /// mass
    double bigSeesawMass = 139.49e-3;
    /// initial center of mass
    /// (use getOrientation.rotate(bigSeesawCOM0) to compute the current center of mass)
    Vec3D bigSeesawCOM0 = {305.03e-3, 186.95e-3, 30e-3};
    /// center of rotation
    Vec3D bigSeesawCOR = {302e-3,179e-3,30e-3};
    /// inertia around com
    ///\todo apply the parallel axis theorem to compute the inertia around the center of roation instead of the center of mass: https://en.wikipedia.org/wiki/Parallel_axis_theorem
    Matrix3D bigSeesawInvInertia = Matrix3D::inverse(
        Matrix3D(70476.55, -37641.24, 0, -37641.24, 329397.59, 0, 0, 0, 296830.02)*1e-9
        );
    /// limits on the angle in the xy-plane the seesaw can make
    /// (necessary to limit the seesaw movement, since inter-wall forces are not computed)
    double bigSeesawAngleZMax = 0;
    double bigSeesawAngleZMin = -23.5*constants::pi/180;
    
    /// Loads the stl file for the big seesaw,
    /// sets the center of rotation to bigSeesawCOR,
    /// sets bigSeesawId, bigSeesawFirstElement.
    void loadBigSeesaw() {
        std::string stlFile = "Big Seasaw.STL";
        Vec3D velocity {0,0,0};
        Vec3D angularVelocity {0,0,0};
        bigSeesawId = wallHandler.readTriangleWall(stlFile,speciesHandler.getLastObject(),1e-3,bigSeesawCOR, velocity, angularVelocity);
        for (auto w : wallHandler) {
            if (w->getGroupId() == bigSeesawId) {
                bigSeesawFirstElement = w;
                break;
            }
        }
    }
    
    /// groupId of the small seesaw (so we can find all the triangles belonging to it)
    unsigned smallSeesawId;
    /// first element in the wallHandler belonging to the small seesaw
    /// (so we can extract information that is equal for all seesaw triangles,
    /// like getOrientation, getAngularVelocity)
    BaseWall* smallSeesawFirstElement = nullptr;
    /// mass
    double smallSeesawMass = 145e-3;
    /// initial center of mass
    /// (use getOrientation.rotate(smallSeesawCOM0) to compute the current center of mass)
    Vec3D smallSeesawCOM0 = {572.02e-3,208.62e-3,30e-3};
    /// center of rotation
    Vec3D smallSeesawCOR = {572.02e-3,208.62e-3,30e-3};
    /// inertia around com
    ///\todo apply the parallel axis theorem to compute the inertia around the center of roation instead of the center of mass: https://en.wikipedia.org/wiki/Parallel_axis_theorem
    Matrix3D smallSeesawInvInertia = Matrix3D::inverse(
        Matrix3D(77503.61, -50757.9, 0, -50757.9, 353140.21, 0, 0, 0, 327073.81)*1e-9
    );
    /// limits on the angle in the xy-plane the seesaw can make
    /// (necessary to limit the seesaw movement, since inter-wall forces are not computed)
    double smallSeesawAngleZMax = 6.0*constants::pi/180;
    double smallSeesawAngleZMin = 0;
    
    /// Loads the stl file for the small seesaw,
    /// sets the center of rotation to smallSeesawCOR,
    /// sets smallSeesawId, smallSeesawFirstElement.
    void loadSmallSeesaw() {
        std::string stlFile = "Small Seasaw.STL";
        Vec3D velocity {0,0,0};
        Vec3D angularVelocity {0,0,0};
        smallSeesawId = wallHandler.readTriangleWall(stlFile,speciesHandler.getLastObject(),1e-3,smallSeesawCOR, velocity, angularVelocity);
        for (auto w : wallHandler) {
            if (w->getGroupId() == smallSeesawId) {
                smallSeesawFirstElement = w;
                break;
            }
        }
    }

    void includeInDomain(Vec3D pos) {
        if (pos.X<getXMin()) {
            setXMin(pos.X);
        } else if (pos.X>getXMax()) {
            setXMax(pos.X);
        }
        if (pos.Y<getYMin()) {
            setYMin(pos.Y);
        } else if (pos.Y>getYMax()) {
            setYMax(pos.Y);
        }
        if (pos.Z<getZMin()) {
            setZMin(pos.Z);
        } else if (pos.Z>getZMax()) {
            setZMax(pos.Z);
        }
    }

    void setupInitialConditions() override
    {
        // set domain such that both particle and walls are included
        setDomain(particle->getPosition(),particle->getPosition());
        for (const auto wall : wallHandler) {
            auto triangle = static_cast<TriangleWall*>(wall);
            includeInDomain(triangle->getVertices()[0]);
            includeInDomain(triangle->getVertices()[1]);
            includeInDomain(triangle->getVertices()[2]);
        }
        logger(INFO,"Simulation domain set to [%,%]x[%,%],[%,%]",
               getXMin(),getXMax(),getYMin(),getYMax(),getZMin(),getZMax());
        // set restitution and contact time
        double gravityTimeScale = sqrt(2.0*particle->getRadius()/getGravity().getLength());
        double collisionTime = gravityTimeScale/20; //ensures small overlaps
        species->setCollisionTimeAndRestitutionCoefficient(collisionTime,restitutionCoefficient, getParticleMass());
        species->setSlidingDissipation(2./7.*species->getDissipation());
        species->setRollingDissipation(2./5.*species->getDissipation());
        species->setTorsionDissipation(2./5.*species->getDissipation());
        species->setSlidingStiffness(2./7.*species->getStiffness());
        species->setRollingStiffness(2./5.*species->getStiffness());
        species->setTorsionStiffness(2./5.*species->getStiffness());
        // set time step
        setTimeStep(collisionTime/25.0);
        logger(INFO,"Simulating for %s with a time step of %s",getTimeMax(),getTimeStep());
        // output
        removeOldFiles();
        setSaveCount(0.02/getTimeStep());
        fStatFile.setFileType(FileType::NO_FILE);
        restartFile.writeFirstAndLastTimeStep();
        setParticlesWriteVTK(true);
        setWallsWriteVTK(FileType::MULTIPLE_FILES);
    }
    
    /// Computes the movement of the two seesaws
    void actionsAfterTimeStep() override
    {
        // if bigSeesaw is set
        if (bigSeesawFirstElement) {
            // branch vector of the COM
            Vec3D branch = bigSeesawCOM0-bigSeesawCOR;
            // \todo gravity torque should act on rotated branch, but something's wrong
            //bigSeesawFirstElement->getOrientation().rotate(branch);
            // torque due to gravity
            Vec3D torque = Vec3D::cross(branch, bigSeesawMass * 0.1 * getGravity());
            Vec3D torque0 = torque;
            // add torques due to particle contacts
            for (auto w : wallHandler) {
                if (w->getGroupId() == bigSeesawId) {
                    torque += w->getTorque();
                }
            }
            Vec3D angularAcceleration = bigSeesawInvInertia * torque;
            Vec3D angularVelocity = bigSeesawFirstElement->getAngularVelocity()
                + angularAcceleration * getTimeStep();
    
            //limit angle around the z-axis [-pi,pi]
            double angleZ = bigSeesawFirstElement->getOrientation().getAngleZ();
            if (angleZ<bigSeesawAngleZMin && angularVelocity.Z>0) {
                angularVelocity = -0.1*angularVelocity;
                //logger(INFO,"Bump off left");
            } else if (angleZ>bigSeesawAngleZMax && angularVelocity.Z<0) {
                angularVelocity = -0.1*angularVelocity;
                //logger(INFO,"Bump off right");
            }
            //logger(INFO, "T % % V % A % B %",torque.Z-torque0.Z, torque0.Z, angularVelocity.Z, angleZ, branch);
            // apply angular velocity
            for (auto w : wallHandler) {
                if (w->getGroupId() == bigSeesawId) {
                    w->setAngularVelocity(angularVelocity);
                }
            }
        }
    
        // if smallSeesaw is set
        if (smallSeesawFirstElement) {
            // branch vector of the COM
            Vec3D branch = smallSeesawCOM0-smallSeesawCOR;
            // \todo gravity torque should act on rotated branch, but something's wrong
            //smallSeesawFirstElement->getOrientation().rotate(branch);
            // torque due to gravity
            Vec3D torque = Vec3D::cross(branch, smallSeesawMass * 0.1 * getGravity());
            Vec3D torque0 = torque;
            // add torques due to particle contacts
            for (auto w : wallHandler) {
                if (w->getGroupId() == smallSeesawId) {
                    torque += w->getTorque();
                }
            }
            Vec3D angularAcceleration = smallSeesawInvInertia * torque;
            Vec3D angularVelocity = smallSeesawFirstElement->getAngularVelocity()
                                    + angularAcceleration * getTimeStep();
        
            //limit angle around the z-axis [-pi,pi]
            double angleZ = smallSeesawFirstElement->getOrientation().getAngleZ();
            if (angleZ<smallSeesawAngleZMin && angularVelocity.Z>0) {
                angularVelocity = -0.1*angularVelocity;
                //logger(INFO,"Bump off left");
            } else if (angleZ>smallSeesawAngleZMax && angularVelocity.Z<0) {
                angularVelocity = -0.1*angularVelocity;
                //logger(INFO,"Bump off right");
            }
            //logger(INFO, "T % % V % A % B %",torque.Z-torque0.Z, torque0.Z, angularVelocity.Z, angleZ, branch);
            // apply angular velocity
            for (auto w : wallHandler) {
                if (w->getGroupId() == smallSeesawId) {
                    w->setAngularVelocity(angularVelocity);
                }
            }
        }
    }
};

int main() {

    // Set up a problem of type MarbleRun
    MarbleRun dpm;
    // Set name of output files
    dpm.setName("MarbleRun");
    // Set name of output files
    dpm.loadSTLFile("KnikkerbaanV12_compleet_without_seasaw.STL");
    dpm.loadBigSeesaw();
    dpm.loadSmallSeesaw();
    // Set physical particle properties
    dpm.setParticlePosition(Vec3D(0.6,0.8,0.03));
    //dpm.setParticlePosition(Vec3D(0.24,0.23,0.03)); //position above big seesaw
    //dpm.setParticlePosition(Vec3D(0.52,0.25,0.03)); //position above small seesaw
    dpm.setParticleRadius(0.01);
    // Set material particle properties
    dpm.setParticleDensity(1000);
    // Set contact properties
    dpm.setSlidingFrictionCoefficient(0.5);
    dpm.setRollingFrictionCoefficient(1e-4);
    dpm.setTorsionFrictionCoefficient(0.0);
    dpm.setRestitutionCoefficient(0.5);
    // set gravity direction
    dpm.setGravity(Vec3D(0,-9.8,0));
    // Set simulation time
    dpm.setTimeMax(2.5);
    // start the solver
    dpm.solve();
}
