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

//based on /storage2/usr/people/sluding/MDCC/C3DshearXL30/MU0_LONG2
#include "Mercury3D.h"
#include "Species/LinearViscoelasticSlidingFrictionSpecies.h"
#include "Boundaries/PeriodicBoundary.h"
#include "Walls/AxisymmetricIntersectionOfWalls.h"

class CSCWalls : public Mercury3D {
public:
    CSCWalls (Mdouble width, Mdouble length, Mdouble height, Mdouble wallThickness)
    {
        eneTolerance = 1.0;
        sizeDistribution = cbrt(2.0);

        setXMin(-width/2-wallThickness);
        setXMax( width/2+wallThickness);
        setYMin(0);
        setYMax(length);
        setZMin(-wallThickness);
        setZMax(height);

        std::cout << "Creating flat-wall box of"
        " width ["  << getXMin() << ", " << getXMax() << "],"
        " length [" << getYMin() << ", " << getYMax() << "],"
        " height [" << getZMin() << ", " << getZMax() << "]." << std::endl;

        //set default name
        setName("CSCWalls");
        setXBallsAdditionalArguments("-v0 -solidf -3dturn 0");

        //set gravity (extra strong since we have few particle layers)
        setGravity(Vec3D(0.0,0.0,-1));

        //create new species
        species =speciesHandler.copyAndAddObject(LinearViscoelasticSlidingFrictionSpecies());
        species->setDensity(6.0/constants::pi);

        //set inter-particle contact properties
        species->setStiffness(2e5/100.0);
        species->setDissipation(25.0/10.0);
        species->setSlidingStiffness(2.0/7.0* species->getStiffness());
        species->setSlidingDissipation(species->getDissipation());
        species->setSlidingFrictionCoefficient(0.5);

        setTimeStep(0.02*species->getCollisionTime(species->getMassFromRadius(0.5)));
        std::cout << "recommended time step: " << getTimeStep() << std::endl;
        setTimeStep(5.0 * getTimeStep());
        std::cout << "time step use for creating initial conditions: " 
            << getTimeStep() << std::endl;
 
        //set time step
        setTimeMax(1e20);
        setSaveCount(2000);

        //set boundaries
        PeriodicBoundary b;
        b.set(Vec3D(0,1,0),getYMin(),getYMax());
        boundaryHandler.copyAndAddObject(b);

        //set walls
        InfiniteWall w;
        w.set(Vec3D(1,0,0), getMax());
        wallHandler.copyAndAddObject(w);
        w.set(Vec3D(-1,0,0), getMin());
        wallHandler.copyAndAddObject(w);
        w.set(Vec3D(0,0,-1), getMin());
        wallHandler.copyAndAddObject(w);

        //remove wall thickness
        setXMin(-width/2);
        setXMax( width/2);
        setZMin(0);
    }

    bool continueSolve() const
    {
        static unsigned int counter = 0;
        if (++counter>100)
        {
            counter=0;
            if (getKineticEnergy()<eneTolerance*getElasticEnergy())
                return false;
        }
        return true;
    }

    void printTime() const
    {
        std::cout << "t=" << getTime() 
            << " Ene " << getKineticEnergy()/getElasticEnergy() << std::endl;
    }

    //add flow particles
    void setupInitialConditions() 
    {
        //number of particles to be inserted
        unsigned int n=(getXMax()-getXMin())*(getYMax()-getYMin())
            *1.8*(getZMax()-getZMin());
        std::cout << "Inserting " << n << " particles" << std::endl;
        //try to find new insertable particles
        unsigned int i = 0;
        SphericalParticle p;
        p.setSpecies(species);
        Mdouble s = sizeDistribution;
        Mdouble rMin=cbrt(0.5/(s*s + 1.0)/(s + 1.0));
        p.setRadius(s*rMin);
        std::cout << "Particle sizes from " << rMin << " to " << s*rMin 
            << " (sizeDistribution " << sizeDistribution << ")" << std::endl;
        Vec3D position;
        Mdouble a=0.0;
        while (i<n){
            position.X = random.getRandomNumber(getXMin()+p.getRadius(), 
                                                getXMax()-p.getRadius());
            position.Y = random.getRandomNumber(getYMin()+p.getRadius(), 
                                                getYMax()-p.getRadius());
            position.Z = random.getRandomNumber(getZMin()+p.getRadius(), 
                                                a*getZMax()-p.getRadius());
            p.setPosition(position);
            if (checkParticleForInteraction(p)) {
                //std::cout << i << std::endl;
                particleHandler.copyAndAddObject(p);
                p.setRadius(random.getRandomNumber(rMin, s*rMin));
                i++;
                if (particleHandler.getNumberOfObjects()%100==0) { 
                    std::cout << " " << particleHandler.getNumberOfObjects()/100; 
                    std::flush(std::cout);
                }
                
            } else { 
                a += 0.00001;
            }
        }
        std::cout << std::endl <<"Inserted " << n << " particles within z<" << a*getZMax() << std::endl;
    }

    void saveWalls()
    {
        std::cout << "Creating flat-wall box of"
        " width ["  << getXMin() << ", " << getXMax() << "],"
        " length [" << getYMin() << ", " << getYMax() << "],"
        " height [" << getZMin() << ", " << getZMax() << "]." << std::endl;

        //only keep wall particles
        for (BaseParticle* p: particleHandler)
        {
            Mdouble x = p->getPosition().X;
            Mdouble z = p->getPosition().Z;
            if ( x > getXMax() || x < getXMin() || z < getZMin() ) 
            {
                if (z<getZMax())
                {
                    p->fixParticle();
                }
            }
        }
        for (int i = particleHandler.getNumberOfObjects()-1; i>=0; i--)
        {
            if (!particleHandler.getObject(i)->isFixed())
                particleHandler.removeObject(i);
        }
        interactionHandler.clear();

        std::cout << "kept " << particleHandler.getNumberOfObjects() 
            << " fixed particles" << std::endl;

        //save data and restart file of last time step
        dataFile.open();
        outputXBallsData(dataFile.getFstream());
        dataFile.close();
        
        restartFile.open();
        writeRestartFile();
        restartFile.close();
    }

    LinearViscoelasticSlidingFrictionSpecies* species;
    Mdouble eneTolerance;
    Mdouble sizeDistribution;
};

int main(int argc, char *argv[]) {
    Mdouble width = 30; //in particle diameter
    Mdouble length = 20;//in particle diameter
    Mdouble height = 30;//in particle diameter
    Mdouble wallThickness = 1;//in particle diameter
    CSCWalls SC(width, length, height, wallThickness);
    SC.solve(argc, argv);
    SC.saveWalls();
    return 0;
}
