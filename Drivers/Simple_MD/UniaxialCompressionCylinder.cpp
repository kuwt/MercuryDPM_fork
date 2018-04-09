//Copyright (c) 2013-2018, The MercuryDPM Developers Team. All rights reserved.
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
#include "CylindricalWall.h"


///This class performs and uniaxial compression stage in a cyclinder.
///It actually comprisses of 5 simulation stages
///Relaxation stage 1, particles are created at random position in an initial domain and allowed to relax during TimeSteps1 at PackingFraction1
///Transion stage 1->2, in TimeSteps12 the system is compacted from PackingFraction1 to PackingFraction2 
///Relaxation stage 2, the system is relaxed during TimeSteps2 at PackingFraction2
///Transion stage 2->3, in TimeSteps23 the system is compacted from PackingFraction2 to PackingFraction3
///Relaxation stage 3, the system is relaxed during TimeSteps3 at PackingFraction3

class UniaxialCompressionCylinder : public Mercury3D
{
public:

    void setupInitialConditions()
    {
        //Set the initial positions and velocities of the particles
        particleHandler.clear();
        BaseParticle p0;
        for (unsigned int i=0;i<numberOfParticles_;i++)
        {
            p0.setRadius(random.getRandomNumber(minParticleRadius_,maxParticleRadius_));
            p0.setPosition(Vec3D(random.getRandomNumber(-domainRadius_+p0.getRadius(),domainRadius_-p0.getRadius()),random.getRandomNumber(-domainRadius_+p0.getRadius(),domainRadius_-p0.getRadius()),random.getRandomNumber(getZMin()+p0.getRadius(),getZMax()-p0.getRadius())));
            p0.setVelocity(Vec3D(random.getRandomNumber(-maxInitialVelocity_,maxInitialVelocity_),random.getRandomNumber(-maxInitialVelocity_,maxInitialVelocity_),random.getRandomNumber(-maxInitialVelocity_,maxInitialVelocity_)));
            //Make sure the particle fits in the initial domain by creating initail positions untill it fits into the domain.
            while( pow(p0.getPosition().X,2)+pow(p0.getPosition().Y,2)<pow(domainRadius_-p0.getRadius(),2))
            {
                p0.setPosition(Vec3D(random.getRandomNumber(-domainRadius_+p0.getRadius(),domainRadius_-p0.getRadius()),random.getRandomNumber(-domainRadius_+p0.getRadius(),domainRadius_-p0.getRadius()),random.getRandomNumber(getZMin()+p0.getRadius(),getZMax()-p0.getRadius())));
            }
            p0.computeMass(Species);
            particleHandler.copyAndAddObject(p0);
        }
        
        //Set the initial position of the walls
        wallHandler.clear();
        CylindricalWall cylindricalWall;
        cylindricalWall.set(domainRadius_);
        wallHandler.copyAndAddObject(cylindricalWall);
        InfiniteWall topWall,bottomWall;
        bottomWall.set(Vec3D( 0, 0,-1), -getZMin());
        wallHandler.copyAndAddObject(bottomWall);
        topWall.set(Vec3D( 0, 0, 1),  getZMax());
        wallHandler.copyAndAddObject(topWall);
    }
    
    
    void actionsBeforeTimeStep()
    {
        //Makes the top wall move according to in which state the process is in
        if(getTime()/getTimeStep()<=timeSteps1_)
        {
            //Initial relaxtation
            //double q = (getTime()/getTimeStep())/timeSteps1_;
            setZMax(z1_);
            wallHandler.getObject(2)->move(getZMax());
            return;
        }
        else if(getTime()/getTimeStep()<=timeSteps1_+timeSteps12_)
        {
            //First compression stage
            double q = (getTime()/getTimeStep()-timeSteps1_)/timeSteps12_;
            setZMax(0.5*((z1_+z2_)+(z1_-z2_)*cos(q*constants::pi)));
            wallHandler.getObject(2)->move(getZMax());
            return;            
        }
        else if(getTime()/getTimeStep()<=timeSteps1_+timeSteps12_+timeSteps2_)
        {
            //Relaxation at middle packing fraction
            //double q = (getTime()/getTimeStep()-(timeSteps1_+timeSteps12_))/timeSteps2_;
            setZMax(z2_);
            wallHandler.getObject(2)->move(getZMax());
            return;
        }
        else if(getTime()/getTimeStep()<=timeSteps1_+timeSteps12_+timeSteps2_+timeSteps23_)
        {
            //Second compression stage
            double q = (getTime()/getTimeStep()-(timeSteps1_+timeSteps12_+timeSteps2_))/timeSteps23_;
            setZMax(0.5*((z2_+z3_)+(z2_-z3_)*cos(q*constants::pi)));
            wallHandler.getObject(2)->move(getZMax());
            return;
        }
        else
        {
            //Final relaxation
            //double q = (getTime()/getTimeStep()-(timeSteps1_+timeSteps12_+timeSteps2_+timeSteps23_))/timeSteps3_;
            setZMax(z3_);
            wallHandler.getObject(2)->move(getZMax());
            return;    
        }
    }
    
    //Set functions for all variables
    void setNumberOfParticles           (int numberOfParticles)             {numberOfParticles_=numberOfParticles;}
    void setMinParticleRadius           (double minParticleRadius)          {minParticleRadius_=minParticleRadius;}
    void setMaxParticleRadius           (double maxParticleRadius)          {maxParticleRadius_=maxParticleRadius;}
    void setMaxInitialVelocity          (double maxInitialVelocity)         {maxInitialVelocity_=maxInitialVelocity;}
    void setDomainRadius                (double domainRadius)               {domainRadius_=domainRadius;}
    void setPackingFraction1            (double packingFraction1)           {packingFraction1_=packingFraction1;}
    void setPackingFraction2            (double packingFraction2)           {packingFraction2_=packingFraction2;}
    void setPackingFraction3            (double packingFraction3)           {packingFraction3_=packingFraction3;}
    void setTimeSteps1                  (int timeSteps1)                    {timeSteps1_=timeSteps1;}
    void setTimeSteps12                 (int timeSteps12)                   {timeSteps12_=timeSteps12;}
    void setTimeSteps2                  (int timeSteps2)                    {timeSteps2_=timeSteps2;}
    void setTimeSteps23                 (int timeSteps23)                   {timeSteps23_=timeSteps23;}
    void setTimeSteps3                  (int timeSteps3)                    {timeSteps3_=timeSteps3;}
    void setCoefficientOfRestitution    (double coefficientOfRestitution)   {coefficientOfRestitution_=coefficientOfRestitution;}
    void setCollisionTime               (double collisionTime)              {collisionTime_=collisionTime;}
    
    void applySettings()
    {
        setCollisionTimeAndRestitutionCoefficient(collisionTime_,coefficientOfRestitution_,get_Mass_from_Radius(minParticleRadius_));
        setTimeStep_by_mass(get_Mass_from_Radius(minParticleRadius_));
        setGravity(Vec3D(0,0,0));
        
        double averageParticleVolume=constants::pi/3.0*(minParticleRadius_+maxParticleRadius_)*(pow(minParticleRadius_,2)+pow(maxParticleRadius_,2));
        z1_=(averageParticleVolume*numberOfParticles_/(constants::pi*pow(domainRadius_,2)*packingFraction1_));
        z2_=(averageParticleVolume*numberOfParticles_/(constants::pi*pow(domainRadius_,2)*packingFraction2_));
        z3_=(averageParticleVolume*numberOfParticles_/(constants::pi*pow(domainRadius_,2)*packingFraction3_));
        setXMax( domainRadius_);
        setXMin(-domainRadius_);
        setYMax( domainRadius_);
        setYMin(-domainRadius_);
        setZMax( z1_);
        setZMin( 0.0);
        
        setTimeMax((timeSteps1_+timeSteps12_+timeSteps2_+timeSteps23_+timeSteps3_)*getTimeStep());            
    }


private:
    unsigned int numberOfParticles_;
    double minParticleRadius_, maxParticleRadius_, maxInitialVelocity_, domainRadius_;
    double packingFraction1_,packingFraction2_,packingFraction3_;
    int timeSteps1_,timeSteps12_,timeSteps2_,timeSteps23_,timeSteps3_;
    double z1_,z2_,z3_;
    double coefficientOfRestitution_;
    double collisionTime_;
    
};

int main(int /*argc*/, char **/*argv[]*/)
{
    UniaxialCompressionCylinder problem;
    
    //Set the cylinder size
    problem.setDomainRadius(5.0);
    
    //Set the 3 different packingfractions
    problem.setPackingFraction1(0.3);
    problem.setPackingFraction2(0.6);
    problem.setPackingFraction3(1.0);
    
    //Set the number of timesteps each stage of the simulation takes
    problem.setTimeSteps1 (1e3);
    problem.setTimeSteps12(1e4);
    problem.setTimeSteps2 (1e4);
    problem.setTimeSteps23(1e4);
    problem.setTimeSteps3 (1e3);
    
    //Set particle properties
    //Coefficent of restitution and collision time are set for the smalles particle, note however the internally the stiffness and damping coefficients are stored
    problem.setNumberOfParticles(100);
    problem.setMinParticleRadius(1.0);
    problem.setMaxParticleRadius(2.0);
    problem.setMaxInitialVelocity(1.0);
    problem.setDensity(1000);
    problem.setCoefficientOfRestitution(0.8);
    problem.setCollisionTime(1e-5);
    
    //Set the saving interval for all data
    problem.setSaveCount(100);

    problem.set_name("UniaxialCompressionCylinder");        
    problem.applySettings();
    problem.solve();
}
