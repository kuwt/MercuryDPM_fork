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

/// if you restart this code the third argument will be used as the number of large particles to add and the forth the number of small.

#include <sstream>
#include <iostream>
#include <iomanip>
#include <cmath>
#include "Species/LinearViscoelasticFrictionSpecies.h"
#include "Chute.h"
#include "Walls/InfiniteWall.h"
#include "Boundaries/PeriodicBoundary.h"
#include "Boundaries/CubeInsertionBoundary.h"
#include <CG/CG.h>

class Chutebelt : public Chute
{
public:
    Chutebelt() : radius_s(0.5)
    {
    }
    
    void setupInitialConditions()
    {      
        //Number of small particles
        int Ns = num_small;
        int Nl = num_large;

        setZMax(200.0*radius_l);


        //Set up live statistics
        //CG<CGCoordinates::Z> cg; //declare a cg object
        //cg.setN(50); //set the CG to be 100 by 100.
        //cg.setWidth(5.0*radius_s);
        //cg.statFile.setSaveCount(200); //set parameters such as the output frequency
        //cgHandler.copyAndAddObject(cg); // add the CG object to the cgHandler*/
        
        //setInflowParticleRadius(radius_s, radius_l);
        
        double rho = 6.0 / constants::pi;
        
        double tc = 1.0 / 200.0;
        double r = 0.1;
        
        double mass_small = 4.0 / 3.0 * constants::pi * pow(radius_s, 3.0) * rho;
        double mass_large = 4.0 / 3.0 * constants::pi * pow(radius_l, 3.0) * rho;
        
        auto species0 = speciesHandler.copyAndAddObject(LinearViscoelasticFrictionSpecies());
        auto species1 = speciesHandler.copyAndAddObject(species0);
        auto species2 = speciesHandler.copyAndAddObject(species0);
        auto species01 = speciesHandler.getMixedObject(species0, species1);
        auto species02 = speciesHandler.getMixedObject(species0, species2);
        auto species12 = speciesHandler.getMixedObject(species1, species2);
        
        
        //  Set the contact time (tc), resitution coefficeient (r) and density (rho) for small for all particles
        
        /////////
        //
        species0->setDensity(rho);
        species0->setCollisionTimeAndRestitutionCoefficient(tc, r, mass_small); //MD::setCollisionTimeAndRestitutionCoefficient(tc, r,mass_small);
        species0->setSlidingDissipation(species0->getDissipation()); //  Set the tangential dissipation equal to the normal disipation for small-small collsions
        species0->setSlidingStiffness(species0->getStiffness()*2.0/7.0);
        species0->setSlidingFrictionCoefficient(0.5);
        //////
        //
        species1->setDensity(rho);
        species1->setCollisionTimeAndRestitutionCoefficient(tc, r, mass_small);
        species1->setSlidingDissipation(species1->getDissipation()); //  Set the tangential dissipationequal to the normal disipation for large-large collision
        species1->setSlidingStiffness(species1->getStiffness()*2.0/7.0);
        species1->setSlidingFrictionCoefficient(0.5);
        //
        //////////
        species2->setDensity(rho);
        species2->setCollisionTimeAndRestitutionCoefficient(tc, r, mass_large);
        species2->setSlidingDissipation(species2->getDissipation());
        species2->setSlidingStiffness(species2->getStiffness()*2.0/7.0);
        species2->setSlidingFrictionCoefficient(0.5);
        

        
   
        
      
        //
        //////////
        //
        species01->setCollisionTimeAndRestitutionCoefficient(tc, r, mass_small, mass_small);
        species01->setSlidingDissipation(species01->getDissipation()); //  Set the tangential dissipation equal to the normal disipation for mixed collision
        species01->setSlidingFrictionCoefficient(0.5);
        species01->setSlidingStiffness(species01->getStiffness()*2.0/7.0);
                
        species02->setCollisionTimeAndRestitutionCoefficient(tc, r, mass_small, mass_large);
        species02->setSlidingDissipation(species01->getDissipation()); //  Set the tangential dissipation equal to the normal disipation
        species02->setSlidingFrictionCoefficient(0.5);
        species02->setSlidingStiffness(species02->getStiffness()*2.0/7.0);
                
        species12->setCollisionTimeAndRestitutionCoefficient(tc, r, mass_small, mass_large);
        species12->setSlidingDissipation(species01->getDissipation()); //  Set the tangential dissipation equal to the normal disipation
        species12->setSlidingFrictionCoefficient(0.5);
        species12->setSlidingStiffness(species12->getStiffness()*2.0/7.0);
                
        std::cout << "Number of large particles:" << Nl << std::endl;
        std::cout << "Number of small particles:" << Ns << std::endl;
        
        Chute::setupInitialConditions();

        
        // Remove the two existing boundaries insertion and periodic and put the peroidic back if perodic else put
        boundaryHandler.clear();
        

        
        //Now add two extra solid walls at the end
        //InfiniteWall w0;
        gateLeft_ = new InfiniteWall;
        //w0.set(Vec3D(-1.0, 0.0, 0.0), Vec3D(getXMin() + 2 * getFixedParticleRadius(),0,0));
        //wallHandler.copyAndAddObject(w0);
        gateLeft_->set(Vec3D(1.0, 0.0, 0.0), Vec3D(getXMax(),0,0));
        gateLeft_->setSpecies(species0);
        wallHandler.addObject(gateLeft_);

        gateRight_ = new InfiniteWall;
        gateRight_->set(Vec3D(-1.0, 0.0, 0.0), Vec3D(getXMin() ,0,0));
        gateRight_->setSpecies(species0);
        wallHandler.addObject(gateRight_);
        //wallHandler.copyAndAddObject(w0);

        if (getIsPeriodic())
        {
            //Periodic in y
            PeriodicBoundary b0;
            b0.set(Vec3D(0.0, 1.0, 0.0), getYMin(), getYMax());
            boundaryHandler.copyAndAddObject(b0);

            //Periodic in x
           // b0.set(Vec3D(1.0, 0.0, 0.0), getXMin(), getXMax());
           // boundaryHandler.copyAndAddObject(b0);
        }

        
        
        // CREATE THE PARTICLES

        SphericalParticle p0;
        smallInsertionBoundary = new CubeInsertionBoundary();
        largeInsertionBoundary = new CubeInsertionBoundary();

        double volumeToInsert=num_small*4.0/3.0*constants::pi*pow(1.0,3);

        logger(INFO,"Need to insert a volume of % of each particle type",volumeToInsert);


        smallInsertionBoundary->set(p0,1000000,Vec3D(getXMin(),getYMin(),getZMax()/2.0),Vec3D(getXMax(),getYMax(),getZMax()),Vec3D(0.0,0.0,0.0),Vec3D(0.0,0.0,0.0),1.0*0.95,1.0*1.05);
        smallInsertionBoundary->setInitialVolume(volumeToInsert);
        smallInsertionBoundary->deactivate();
        boundaryHandler.addObject(smallInsertionBoundary);

        largeInsertionBoundary->set(p0,1000000,Vec3D(getXMin(),getYMin(),getZMin()),Vec3D(getXMax(),getYMax(),getZMax()/2.0),Vec3D(0.0,0.0,0.0),Vec3D(0.0,0.0,0.0),2.0*0.95,2.0*1.05);
        largeInsertionBoundary->setInitialVolume(volumeToInsert);
        largeInsertionBoundary->deactivate();
        boundaryHandler.addObject(largeInsertionBoundary);
       /* while ((Ns > 0) || (Nl > 0))
        {
            
            //random to see if want to generate a large or small particles, helps makes the initial conditions homogenious
            if (random.getRandomNumber(1.0, Nl + Ns) > Nl)
            {
                P0.setRadius(radius_s);
                P0.setSpecies(species1);
                Ns--;
            }
            else
            {
                P0.setRadius(radius_l);
                P0.setSpecies(species2);
                Nl--;
            }
            //randomise particle position, zero intial velocity
	    do 
                {
                P0.setPosition(
                    Vec3D(random.getRandomNumber(getXMin() + radius_l + 2 * getFixedParticleRadius(),
                            getXMax() - radius_l - 2 * getFixedParticleRadius()),
                            random.getRandomNumber(getYMin() + radius_l, getYMax() - radius_l),
                            random.getRandomNumber(getZMin() + radius_l, getZMax() - radius_l)));
                P0.setVelocity(Vec3D(0.0, 0.0, 0.0));
            
                 }
             while (!checkParticleForInteraction(P0));
	     particleHandler.copyAndAddObject(P0);

        }*/
        
        std::cout << "Finished creating particles" << std::endl;
        
        for (int i = 0; i < particleHandler.getNumberOfObjects(); i++)
        {
            BaseParticle* P0 = particleHandler.getObject(i);
            if (P0->getIndSpecies() == 0)
            {
                BedParticles.push_back(P0); 
            }
        }
        std::cout << "Finished storing bed particles" << std::endl;
    }
    
    void actionsOnRestart()
    {
        
        int Ns=num_restart_small;
        int Nl=num_restart_large;
        
        auto species0=speciesHandler.getObject(0);
        auto species1=speciesHandler.getObject(1);
        auto species2=speciesHandler.getObject(2);
        
        std::cout << "Restarting and adding " << Ns << " small particles and " << Nl <<" large particles" << std::endl;
        
        
        SphericalParticle P0;
        while ((Ns > 0) || (Nl > 0))
        {
            
            //random to see if want to generate a large or small particles, helps makes the initial conditions homogenious
            if (random.getRandomNumber(1.0, Nl + Ns) > Nl)
            {
                P0.setRadius(radius_s);
                P0.setSpecies(species1);
                Ns--;
            }
            else
            {
                P0.setRadius(radius_l);
                P0.setSpecies(species2);
                Nl--;
            }
            //randomise particle position, zero intial velocity
            do
            {
                P0.setPosition(
                               Vec3D(random.getRandomNumber(getXMin() + radius_l + 2 * getFixedParticleRadius(),
                                                            getXMax() - radius_l - 2 * getFixedParticleRadius()),
                                     random.getRandomNumber(getYMin() + radius_l, getYMax() - radius_l),
                                     random.getRandomNumber(getZMin() + radius_l, getZMax() - radius_l)));
                P0.setVelocity(Vec3D(0.0, 0.0, 0.0));
                
            }
            while (!checkParticleForInteraction(P0));
            particleHandler.copyAndAddObject(P0);
            
        }
        
        std::cout << "Finished adding new particles. Now simulation will restart" << std::endl;
        
        
    }

    void actionsAfterTimeStep()
    {

        setParticlesWriteVTK(true);

        switch (step)
        {
            case 0 :
                logger(INFO,"Finished create base now moving to insert large particles");
                setChuteAngle(0.0);
                step++;
                break;

            case 1 :
                largeInsertionBoundary->activate();
                /// \todo Work out why the volume inserted is not what I expext
                logger(DEBUG,"Inserted volume of % of larger particles, tagget is %",largeInsertionBoundary->getMassOfParticlesInserted(),largeInsertionBoundary->getInitialVolume());

                if (getTime()>50.0)
                {
                    logger(INFO,"Moving on to insert small particles");
                    largeInsertionBoundary->deactivate();
                    smallInsertionBoundary->activate();
                    step++;
                }

                break;

            case 2:

                if (getTime()>100.0)
                {
                    logger(INFO,"Now starting the simulations");
                    step++;
                    setTime(0.0);
                    setZMax(100.0*radius_l);
                    smallInsertionBoundary->deactivate();
                    setChuteAngle(25.0);
                    //To not clear hte boundary is this includes the periodic boundaries.
                    //boundaryHandler.clear();
                    //wallHandler.clear();

                    wallHandler.removeObject(gateLeft_->getIndex());
                    if (getIsPeriodic())
                    {
                        wallHandler.removeObject(gateRight_->getIndex());
                        //Periodic in y
                        PeriodicBoundary b0;
                        //b0.set(Vec3D(0.0, 1.0, 0.0), getYMin(), getYMax());
                        //boundaryHandler.copyAndAddObject(b0);

                        //Periodic in x
                        b0.set(Vec3D(1.0, 0.0, 0.0), getXMin(), getXMax());
                        boundaryHandler.copyAndAddObject(b0);
                    }
                }
                break;
        }

    }

    // MOVING BED
    void actionsBeforeTimeStep()
    {
        if (BedParticles.empty())
        {
            for (int i = 0; i < particleHandler.getNumberOfObjects(); i++)
            {
                BaseParticle* P0 = particleHandler.getObject(i);
                if (P0->getIndSpecies() == 0)
                {
                    P0->setVelocity(Vec3D(0.0,0.0,0.0));
                    BedParticles.push_back(P0);
                }
            }
            std::cout << "Finished storing bed particles" << std::endl;
        }
        
        if (getTime() > 1.0)
        {
            for (int i = 0; i < BedParticles.size(); i++)
            {
                
                BaseParticle* P0 = BedParticles[i];
                if (P0->getIndSpecies() == 0)
                {
                    
                    Vec3D position;
                    position = P0->getPosition();
                    
                    double dt = getTimeStep();
                    position.X = position.X - beltSpeed * dt;
                    
                    if (position.X < getXMin())
                    {
                        position.X = position.X - getXMin() + getXMax();
                    }
                    P0->setPosition(position);                  
                }
            }
        }
    }
    
    void set_beltSpeed(double new_speed)
    {
        beltSpeed = new_speed;
    }
    void set_particle_numbers(int new_num_small, int new_num_large)
    {
        
        if (new_num_small > 0)
        {
            num_small = new_num_small;
        }
        else
        {           
            std::cerr << "Please give a positive numnber if small particles" << std::endl;
        }
        
        if (new_num_large > 0)
        {
            num_large = new_num_large;
        }
        else
        {           
            std::cerr << "Please give a positive numnber if small particles" << std::endl;
        }
        
    }
    void set_particle_numbers(int new_num_small)
    {
        if (new_num_small > 0)
        {          
            num_small = new_num_small;
            num_large = std::pow(radius_s / radius_l, 3) * num_small * particleVolRatio;          
        }
        else
        {
            std::cerr << "Please give a positve number of particels" << std::endl;
        }
        
    }
    
    void set_radiusLarge(double new_large_radius)
    {
        
        if (new_large_radius > 0)
        {           
            radius_l = new_large_radius;    
        }
        else
        {           
            std::cerr << "Radius must be greater than zero" << std::endl;
        }
    }
    
    void set_particle_number_volRatio(double new_volume_ratio)
    {
        particleVolRatio = new_volume_ratio;
    }
    
unsigned int num_restart_small;
unsigned int num_restart_large;
    
private:
    std::vector<BaseParticle*> BedParticles;
    double beltSpeed;
    double radius_l;
    const double radius_s;
    int num_small;
    int num_large;
    double particleVolRatio;

    InfiniteWall* gateLeft_;
    InfiniteWall* gateRight_;

    CubeInsertionBoundary* smallInsertionBoundary;
    CubeInsertionBoundary* largeInsertionBoundary;
    /**
     * \brief Tracks where we are up to in setting the problem up
     * \details
     * Step 0 : Insert small particles at the bottom
     * Step 1 : Insert large particles at the top.
     * Both step 0 and 1 setup are done with a gate at the right but the correct chute angle
     */
    unsigned step=0;




};

int main(int argc, char *argv[])
{
    //Print description
    std::cout << std::endl << "Description: A quasi-2D moving-bed channel with walls on the left and right boundary." << std::endl;
    
    // Problem parameters
    Chutebelt problem;

    problem.autoNumber();
    problem.setName("Segregation");
    problem.setTimeMax(2000.0);
    problem.setTimeStep(1. / (200.0 * 50.0));
    
    problem.set_radiusLarge(2.0);
    problem.set_particle_number_volRatio(1.0); //volume ratio of large to small
    problem.set_particle_numbers(2500, 130); //Should be 5000 and 130

    problem.setChuteAngleAndMagnitudeOfGravity(25.0, 1.0);
    problem.set_beltSpeed(0.0);
    
    // Chute properties : Simply remove the first line to add side walls.
    problem.makeChutePeriodic();
    problem.setXMax(150.0);
    problem.setYMax(5.0);
    
    problem.setZMin(0.0);
    problem.setZMax(250.0);

    /// \todo change this to the mean particle size (not small particle size)
    problem.setFixedParticleRadius(0.75);
    
    //Swap the next two lines to swap between the different type of rought bottoms.
    problem.setRoughBottomType(MULTILAYER);
    //problem.setRoughBottomType(MONOLAYER_DISORDERED);

    problem.setSaveCount(helpers::getSaveCountFromNumberOfSavesAndTimeMaxAndTimeStep(2000, problem.getTimeMax(), problem.getTimeStep()));
    
    problem.setXBallsAdditionalArguments("-cmode 7");
    
    problem.readArguments(argc, argv);
    
    if (argc > 4)
    {
        problem.num_restart_large=atoi(argv[3]);
        problem.num_restart_small=atoi(argv[4]);
    }
    else
    {
        problem.num_restart_large=0;
        problem.num_restart_small=0;
    }
    
    std::cout << problem.num_restart_small <<" " <<problem.num_restart_large << std::endl;


    problem.solve();

    
}
