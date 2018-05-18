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

#include <Mercury3D.h>
#include "Species/LinearViscoelasticFrictionSpecies.h"
#include "Walls/InfiniteWall.h"
#include "Boundaries/PeriodicBoundary.h"
using namespace std;
/*
 *This is a shear box code i.e. the two walls on the left (xmin) and right (xmax) rotate there angle.
 * It have a five step filling process
 * Step 1 : Create the walls (very simple at momnet as flat simple walls)
 * Step 2 : Fill particles: insert particles (with step 3)
 * Step 3 : Let the particles settle so more can be inserted not return to step 2
 * Note, if you are doing normally graded insertion step 2 and 3 will be done seperate for large and small particles.
 * Step 4 : This is a full settle until the particles are almost stationary
 * Step 5 : This starts to move the wall (rotate them it is a shear box). Note, this step resets the clock to zero
 */


class ShearBox : public Mercury3D
{
public:
    
    /*!
     * \brief The only consrtuctors
     * \details Sets a few defauls so the code can be called with no parameters
     */
    ShearBox()
    {
        radius_s=0.5;
        radius_l=0.75;
        radius_l2 = radius_l*0.9;
        radius_s2 = radius_s*1.2;
        
        numSmallToBeInserted=100;
        numLargeToBeInserted=100;
        
        segregatedFilling=false;
        
        setShearParameters(0,1.0);
        
        setXMin(0.0);
        setYMin(0.0);
        setZMin(0.0);
        
        setXMax(10.0);
        setYMax(10.0);
        
        showFillingFlag=false;
        
        
    }
    
    
    /*!
     * \breif Required restarter does nothing other that what is required
     * \details  This is the restarted teh code has lambda functions for the moving walls and needs so pointer magic; therefore these need to recreated on restart. Note, in the new version of MercuryDPM this stuff should all be stored in the restart so this will not be used in the future.
     */
    void actionsOnRestart() override
    {
        

        speciesWall=dynamic_cast<LinearViscoelasticFrictionSpecies*>(speciesHandler.getObject(0));
        speciesS1=dynamic_cast<LinearViscoelasticFrictionSpecies*>(speciesHandler.getObject(1));
        speciesS2=dynamic_cast<LinearViscoelasticFrictionSpecies*>(speciesHandler.getObject(2));
        speciesL1=dynamic_cast<LinearViscoelasticFrictionSpecies*>(speciesHandler.getObject(3));
        speciesL2=dynamic_cast<LinearViscoelasticFrictionSpecies*>(speciesHandler.getObject(4));
        
        
        left=dynamic_cast<InfiniteWall*>(wallHandler.getObject(0));
        right=dynamic_cast<InfiniteWall*>(wallHandler.getObject(1));
        
        
        setupMovingWalls();
        
        
    }
    
    /*!
     * \brief Writes the required information for a restart
     * \details Lots of extra parameters that are required for a restart. This makes sure all these are written
     * \todo It does not write what the varibles are to the file. So the files are not human readable at the moment.
     */
    void write(std::ostream& os, bool print_all = false) const override
    {
        
        
        os << step << std::endl;
        os << checkTime << std::endl;
        os << numSmallAtBottom << std::endl;
        os << numSmallToBeInserted << std::endl;
        os << numLargeToBeInserted << std::endl;
        os << segregatedFilling << std::endl;
        os << radius_s << std::endl;
        os << radius_l << std::endl;
        os << radius_s2 << std::endl;
        os << radius_l2 << std::endl;
        os << zMinInsert << std::endl;
        os << shearBoxAngle << std::endl;
        os << shearBoxTimePeriod << std::endl;
        os << saveCountPerCycle << std::endl;
        os << showFillingFlag << std::endl;
        os << volfrac << std::endl;
        DPMBase::write(os, print_all);
        
        
    }
    
    /*!
     * \brief Reader for the custom restart files
     * \details This is the counterpart of ShearBox::read and is used when restarted. Changes in ShearBox::read should be mirroed here
     */
    void read(std::istream& is) override
    {
        
        
        is >> step;
        is >> checkTime;
        is >> numSmallAtBottom;
        is >> numSmallToBeInserted;
        is >> numLargeToBeInserted;
        is >> segregatedFilling;
        is >> radius_s;
        is >> radius_l;
        is >> radius_s2;
        is >> radius_l2;
        is >> zMinInsert;
        is >> shearBoxAngle;
        is >> shearBoxTimePeriod;
        is >> saveCountPerCycle;
        is >> showFillingFlag;
        is >> volfrac;
        std::cout << "The current step is " << step <<std::endl;
        std::cout << "The current time is " << checkTime <<std::endl;
        std::cout << "Small particle volume fraction is " << volfrac <<std::endl;
        
        DPMBase::read(is);
        
        
        
    }
    
    /*!
     * \brief Setups of the walls and particle properties i.e. Step 1
     * \details For this code the particles are fixed here (as determined in this funcitons) where as the wall properties and set by access functions before
     */
    void setupInitialConditions() override
    {
        
        
        //Define the particle properties tc, COR and density
        double tc = 1.0 / 200.0;
        double COR=0.9;
        double rho = 6.0 / constants::pi;
        
        //Set the time step to a senible value bases on the information about
        setTimeStep(tc/50.0);
        
        //Compute the mass of both the large and small particles used to compute spring and damping constants
        radius_l2 = radius_l*0.95;
        radius_s2 = radius_s*1.15;
        
        double mass_small = 4 / 3 * constants::pi * pow(radius_s, 3.0) * rho;
        double mass_small2 = 4 / 3 * constants::pi * pow(radius_s2, 3.0) * rho;
        double mass_large = 4 / 3 * constants::pi * pow(radius_l, 3.0) * rho;
        double mass_large2 = 4 / 3 * constants::pi * pow(radius_l2, 3.0) * rho;
        
        //Set up the Materials properties. Here everythink is the same and standard. Note, the COR is set for the small particles so will be slightly different for the large partilces. Maybe these should be changed in the future.
        
        
        
//        speciesAll = new LinearViscoelasticFrictionSpecies;
//        speciesAll->setDensity(6.0/constants::pi);
//        speciesAll->setCollisionTimeAndRestitutionCoefficient(tc, COR, mass_small);
//        speciesAll->setSlidingDissipation(speciesAll->getDissipation());
//        speciesAll->setSlidingStiffness(speciesAll->getStiffness()*2.0/7.0);
//        speciesAll->setSlidingFrictionCoefficient(0.5);
//        speciesHandler.addObject(speciesAll);
//        
//        speciesWall = new LinearViscoelasticFrictionSpecies;
//        speciesWall->setDensity(6.0/constants::pi);
//        speciesWall->setCollisionTimeAndRestitutionCoefficient(tc, COR, mass_small);
//        speciesWall->setSlidingDissipation(speciesWall->getDissipation());
//        speciesWall->setSlidingStiffness(speciesWall->getStiffness()*2.0/7.0);
//        speciesWall->setSlidingFrictionCoefficient(2.0);
//        speciesHandler.addObject(speciesWall);
        
        //
        speciesWall = new LinearViscoelasticFrictionSpecies;
        speciesWall->setDensity(rho);
        speciesWall->setCollisionTimeAndRestitutionCoefficient(tc, COR, mass_small);
        speciesWall->setSlidingDissipation(speciesWall->getDissipation());
        speciesWall->setSlidingStiffness(speciesWall->getStiffness()*2.0/7.0);
        speciesWall->setSlidingFrictionCoefficient(1.0);
        speciesHandler.addObject(speciesWall);
        
        
        
        speciesS1 = new LinearViscoelasticFrictionSpecies;
        speciesS1->setDensity(rho);
        speciesS1->setCollisionTimeAndRestitutionCoefficient(tc, COR, mass_small);
        speciesS1->setSlidingDissipation(speciesS1->getDissipation());
        speciesS1->setSlidingStiffness(speciesS1->getStiffness()*2.0/7.0);
        speciesS1->setSlidingFrictionCoefficient(0.5);
        speciesHandler.addObject(speciesS1);
        
        speciesS2 = new LinearViscoelasticFrictionSpecies;
        speciesS2->setDensity(rho);
        speciesS2->setCollisionTimeAndRestitutionCoefficient(tc, COR, mass_small2);
        speciesS2->setSlidingDissipation(speciesS2->getDissipation());
        speciesS2->setSlidingStiffness(speciesS2->getStiffness()*2.0/7.0);
        speciesS2->setSlidingFrictionCoefficient(0.5);
        speciesHandler.addObject(speciesS2);
        
        speciesL1 = new LinearViscoelasticFrictionSpecies;
        speciesL1->setDensity(rho);
        speciesL1->setCollisionTimeAndRestitutionCoefficient(tc, COR ,mass_large);
        speciesL1->setSlidingDissipation(speciesL1->getDissipation());
        speciesL1->setSlidingStiffness(speciesL1->getStiffness()*2.0/7.0);
        speciesL1->setSlidingFrictionCoefficient(0.5);
        speciesHandler.addObject(speciesL1);

        speciesL2 = new LinearViscoelasticFrictionSpecies;
        speciesL2->setDensity(rho);
        speciesL2->setCollisionTimeAndRestitutionCoefficient(tc, COR, mass_large2);
        speciesL2->setSlidingDissipation(speciesL2->getDissipation());
        speciesL2->setSlidingStiffness(speciesL2->getStiffness()*2.0/7.0);
        speciesL2->setSlidingFrictionCoefficient(0.5);
        speciesHandler.addObject(speciesL2);

        

        
        auto speciesWallS1 = speciesHandler.getMixedObject(speciesWall, speciesS1);
        auto speciesWallS2 = speciesHandler.getMixedObject(speciesWall, speciesS2);
        auto speciesWallL1 = speciesHandler.getMixedObject(speciesWall, speciesL1);
        auto speciesWallL2 = speciesHandler.getMixedObject(speciesWall, speciesL2);
        auto speciesS1S2 = speciesHandler.getMixedObject(speciesS1, speciesS2);
        auto speciesS1L1 = speciesHandler.getMixedObject(speciesS1, speciesL1);
        auto speciesS1L2 = speciesHandler.getMixedObject(speciesS1, speciesL2);
        auto speciesS2L1 = speciesHandler.getMixedObject(speciesS2, speciesL1);
        auto speciesS2L2 = speciesHandler.getMixedObject(speciesS2, speciesL2);
        auto speciesL1L2 = speciesHandler.getMixedObject(speciesL1, speciesL2);
        
        
        speciesWallS1->setCollisionTimeAndRestitutionCoefficient(tc, COR ,mass_small, mass_small);
        speciesWallS1->setSlidingDissipation(speciesWallS1->getDissipation());
        speciesWallS1->setSlidingFrictionCoefficient(1.0);
        speciesWallS1->setSlidingStiffness(speciesWallS1->getStiffness()*2.0/7.0);
        
        speciesWallS2->setCollisionTimeAndRestitutionCoefficient(tc, COR ,mass_small, mass_small2);
        speciesWallS2->setSlidingDissipation(speciesWallS2->getDissipation());
        speciesWallS2->setSlidingFrictionCoefficient(1.0);
        speciesWallS2->setSlidingStiffness(speciesWallS2->getStiffness()*2.0/7.0);
        
        speciesWallL1->setCollisionTimeAndRestitutionCoefficient(tc, COR ,mass_small, mass_large);
        speciesWallL1->setSlidingDissipation(speciesWallL1->getDissipation());
        speciesWallL1->setSlidingFrictionCoefficient(1.0);
        speciesWallL1->setSlidingStiffness(speciesWallL1->getStiffness()*2.0/7.0);
        
        speciesWallL2->setCollisionTimeAndRestitutionCoefficient(tc, COR ,mass_small, mass_large2);
        speciesWallL2->setSlidingDissipation(speciesWallL2->getDissipation());
        speciesWallL2->setSlidingFrictionCoefficient(1.0);
        speciesWallL2->setSlidingStiffness(speciesWallL2->getStiffness()*2.0/7.0);
        
        
        speciesS1S2->setCollisionTimeAndRestitutionCoefficient(tc, COR, mass_small, mass_small2);
        speciesS1S2->setSlidingDissipation(speciesS1S2->getDissipation());
        speciesS1S2->setSlidingFrictionCoefficient(0.5);
        speciesS1S2->setSlidingStiffness(speciesS1S2->getStiffness()*2.0/7.0);
        
        speciesS1L1->setCollisionTimeAndRestitutionCoefficient(tc, COR, mass_small, mass_large);
        speciesS1L1->setSlidingDissipation(speciesS1L1->getDissipation());
        speciesS1L1->setSlidingFrictionCoefficient(0.5);
        speciesS1L1->setSlidingStiffness(speciesS1L1->getStiffness()*2.0/7.0);
        
        speciesS1L2->setCollisionTimeAndRestitutionCoefficient(tc, COR, mass_small, mass_large2);
        speciesS1L2->setSlidingDissipation(speciesS1L2->getDissipation());
        speciesS1L2->setSlidingFrictionCoefficient(0.5);
        speciesS1L2->setSlidingStiffness(speciesS1L2->getStiffness()*2.0/7.0);
        
        
        speciesS2L1->setCollisionTimeAndRestitutionCoefficient(tc, COR ,mass_small2, mass_large);
        speciesS2L1->setSlidingDissipation(speciesS2L1->getDissipation());
        speciesS2L1->setSlidingFrictionCoefficient(0.5);
        speciesS2L1->setSlidingStiffness(speciesS2L1->getStiffness()*2.0/7.0);
        
        speciesS2L2->setCollisionTimeAndRestitutionCoefficient(tc, COR ,mass_small2, mass_large2);
        speciesS2L2->setSlidingDissipation(speciesS2L2->getDissipation());
        speciesS2L2->setSlidingFrictionCoefficient(0.5);
        speciesS2L2->setSlidingStiffness(speciesS2L2->getStiffness()*2.0/7.0);
        
        
        speciesL1L2->setCollisionTimeAndRestitutionCoefficient(tc, COR ,mass_large, mass_large2);
        speciesL1L2->setSlidingDissipation(speciesL1L2->getDissipation());
        speciesL1L2->setSlidingFrictionCoefficient(0.5);
        speciesL1L2->setSlidingStiffness(speciesL1L2->getStiffness()*2.0/7.0);
        
        
        
        
        //Set the min to zero. This should alread to true but good to check
        setXMin(0.0);
        setYMin(0.0);
        setZMin(0.0);
        
        //Look up the size of the box in x and y i.e. the cross
        double Xmax = getXMax();
        double Ymax = getYMax();
        
        //Calculate (and set) the height such that we have 50% solids fraction on average. This gives nice packing but does mean it fills in serveral steps.
        double Zmax = 2.0*4.0/3.0*constants::pi*(numLargeToBeInserted * pow(radius_l, 3.0) + numSmallToBeInserted * pow(radius_s, 3.0)) / Xmax / Ymax;
        setZMax(Zmax+radius_l);
        
        std::cout << "Zmax:" << Zmax << std::endl;
        
        //Set gravity to magitude 1 (from the used non-dim) and pointing down.
        setGravity(Vec3D(0.0,0.0,-1.0));
        
        
        ////////////////////////////////////
        //Next sections sets up the walls///
        ////////////////////////////////////
        
        // InfiniteWall w0;
        //Put in the front and back "y-walls"
        // w0.set(Vec3D(0.0, -1.0, 0.0), Vec3D(0.0,getYMin(),0.0));
        // w0.setSpecies(speciesAll);
        // wallHandler.copyAndAddObject(w0);
        
        //  w0.set(Vec3D(0.0, 1.0, 0.0), Vec3D(0.0,getYMax(),0.0));
        //  wallHandler.copyAndAddObject(w0);
        
        
        //Put in the front and back periodic "y-walls"
        //
        PeriodicBoundary b0;
        b0.set(Vec3D(0.0, 1.0, 0.0), getYMin(), getYMax());
        boundaryHandler.copyAndAddObject(b0);
        
        
        
        
        
        //Put in the side walls "x-walls"
        left = new InfiniteWall;
        left->setSpecies(speciesWall);
        left->set(Vec3D(-1.0, 0.0, 0.0), Vec3D(getXMin(),0.0,0.0));
        wallHandler.addObject(left);
        
        
        right=new InfiniteWall;
        right->setSpecies(speciesWall);
        right->set(Vec3D(1.0, 0.0, 0.0), Vec3D(getXMax(),0.0,0.0));
        wallHandler.addObject(right);
        
        //This is functions which creates the lambda for the moving walls
        setupMovingWalls();
        
        
        //Put in the base wall "z-wall"
        InfiniteWall wz;
        wz.set(Vec3D(0.0, 0.0, -1.0), Vec3D(0.0,0.0,getZMin()));
        wallHandler.copyAndAddObject(wz);
        
        
        //Move on to step 2 as the walls are created
        step=2;
        
        
        //Flag if you are saving data during the filling process
        if (showFillingFlag)
        {
            
            //If you are saving data save every 0.1 seconds. This is convient w.r.t to how the filling process works
            setSaveCount(0.1/getTimeStep());
        }
        else
        {
            //If you do not want saving set the counter to 2. So it will only save the very first step.
            //Note, you cannot set this to 0 or 1; this seems to be a MecuryDPM bug
            setSaveCount(helpers::getSaveCountFromNumberOfSavesAndTimeMaxAndTimeStep(2, getTimeMax(), getTimeStep()));
        }
        
        //This does one set of filling. Note, this is done everytime step as well so this is not really required.
        createParticles();
        
        
        
    }
    
    
    
    /*!
     * \brief Create and inserts in the drum new partciles random at non-overlapping locations in the drum. If also check the locations is not already filled by a particle
     * \details This has two different filling proceduces one for normal graded and one homogeneously mixed. This is determined by the flag segregatedFilling.
     *  Also this is called done at the end of the initial conditions setup and every time step (while the code is in step 2).
     */
    void createParticles()
    {
        
        
        // Report the current amount of particles that ar still outstanding to be inserted
        std::cout << "\n \n \n";
        std::cout << "STEP 2: Inserting particles " << std::endl;
        std::cout << "---------------------------" << std::endl;
        std::cout << "\n \n \n";
        
        std::cout << "Number of large particles to be inserted:" << numLargeToBeInserted << std::endl;
        std::cout << "Number of small particles to be inserted:" << numSmallToBeInserted << std::endl;
        
        // Create a  paticles to use as the template for creating new particles
        BaseParticle P0;
        
        if (segregatedFilling)
        {
            
            setZMax(getZMax()*0.05);
            //add smalls at the bottom
            
            std::cout << "Number of smalls to the bottom to be added " << numSmallAtBottom <<std::endl;
            
            while (numSmallAtBottom>0)
            {
                for (int n=numSmallAtBottom; n>0; n--)
                {
                    if (numSmallToBeInserted%2 == 0)
                    {
                        //It is a small partile
                        P0.setRadius(radius_s2);
                 
                        P0.setSpecies(speciesS2);
                    }
                    else
                    {
                        //It is a slightly smaller large partile
                        P0.setRadius(radius_s);
                   
                        P0.setSpecies(speciesS1);
                    }
                    

                    
                    if (!particleInsertable(P0)) break;
                    
                    
                    
                    BaseParticle* q = particleHandler.copyAndAddObject(P0);
                    for (BaseBoundary *it : boundaryHandler)
                        it->createPeriodicParticle(q, particleHandler);
                    
                    
                    
                    
                    
                    numSmallToBeInserted--;
                    numSmallAtBottom--;
                    
                    
                    //  std::cout << "numSmallAtBottom added " << numSmallAtBottom <<std::endl;
                }
                
            }
            // set back to correct Zmax afther inserting smalls at bottom
            setZMax(getZMax()*20);
            
            
            std::cout << "Added bottoms smalls"  <<std::endl;
           
     
            
            
            
            
            //If segreated filling first insert the large particles
            if (numLargeToBeInserted>0)
            {
                
                
                //Try to insert all the large particles
                while (numLargeToBeInserted > 0)
                {
                    if (numLargeToBeInserted%2 == 0)
                    {
                        // It is a slightly smaller large partile
                        
                        P0.setRadius(radius_l2);
                        P0.setSpecies(speciesL2);
                        
                    }
                    else
                    {
                        //It is a large partile
                        P0.setRadius(radius_l);
                        P0.setSpecies(speciesL1);
                    }
                    
                    
                    
                    

                    
                    
                    if (!particleInsertable(P0)) break;
                    
                    //  Vec3D position;
                    // position.X=getXMax()*0.5;
                    // position.Y=getYMax()*0.5;
                    //   position.Z=getZMax()*0.2+radius_l*2;
                    //  P0.setPosition(position);
                    
                    BaseParticle* q = particleHandler.copyAndAddObject(P0);
                    
                    for (BaseBoundary *it : boundaryHandler)
                        it->createPeriodicParticle(q, particleHandler);
                    
                    
                    
                    numLargeToBeInserted--;
                    
                }
                step=3;
                
                
            }
            //Once the large all inserted ask about the small
            else
            {
                
                //Try to insert all the small particles at once
                while (numSmallToBeInserted > 0)
                {
                    
                    if (numSmallToBeInserted%2 == 0)
                    {
                        //It is a small partile
                        P0.setRadius(radius_s);
                        P0.setSpecies(speciesS1);
                    }
                    else
                    {
                        //It is a slightly smaller small partile
                        P0.setRadius(radius_s2);
                        P0.setSpecies(speciesS2);
                    }
                    

                    if (!particleInsertable(P0)) break;
                    
                    
                    BaseParticle* q = particleHandler.copyAndAddObject(P0);
                    
                    for (BaseBoundary *it : boundaryHandler)
                        it->createPeriodicParticle(q, particleHandler);
                    
                    numSmallToBeInserted--;
                    
                }
                step=3;
                
                
            }
            
            
            
        }
        //This is the mixex initial conditions case
        else
        {

        }
        
        //end else insert mixed particles
        
        //remove periodic particles
        for (unsigned int i = particleHandler.getNumberOfObjects(); i >= 1; i--)
            if (particleHandler.getObject(i - 1)->getPeriodicFromParticle() != nullptr)
            {
                while (particleHandler.getObject(i - 1)->getInteractions().size() > 0)
                {
                    interactionHandler.removeObjectKeepingPeriodics(particleHandler.getObject(i - 1)->getInteractions().front()->getIndex());
                }
                particleHandler.removeObject(i - 1);
            }
        
        
        
        
        
        
        
        
        //Report how many particles are still left to be inserted after you have run out of space.
        std::cout << "Finished creating particles" << std::endl;
        std::cout << "Number of large particles still to be inserted:" << numLargeToBeInserted << std::endl;
        std::cout << "Number of small particles still to be inserted:" << numSmallToBeInserted << std::endl;
        std::cout << "Number of smalls to the bottom still to be added " << numSmallAtBottom <<std::endl;
        
        
        //If you managed to insert all target particles move on to step 4.
        if ((numSmallToBeInserted==0) && (numLargeToBeInserted==0))
        {
            // If you are here are partices are inserted and you are moving to step 4; relax to very low KE
            step=4;
            
            std::cout << "\n \n \n";
            std::cout << "STEP 4: Relaxing particles " << std::endl;
            std::cout << "---------------------------" << std::endl;
            std::cout << "\n \n \n";
            checkTime=getTime()+1.0;
            
            
            
            
            
        }
        //If you still have particles to insert move on to step 3.
        else
        {
            // If you are here, the drum is not full and you are settling the inserted particles. 10 times greater than step 4. Note, this is just to clear some space for the new particles.
            
            step=3;
            
            
            std::cout << "\n \n \n";
            std::cout << "STEP 3: Settling the inserted particles " << std::endl;
            std::cout << "---------------------------" << std::endl;
            std::cout << "\n \n \n";
            checkTime=getTime()+1.0;
            
        }
        
    }
    
    /*!
     * \brief This functions moves the particles to a place it can be inserted. If no place can be found it returns false.
     * \details The functions looks for a place to insert the particle and returns true if it has been found. It also moves the particle to the point the space was found so it can be directly inserted. Maybe merge with the insertion step itself later.
     * On entry the partilce must have the right species and size. This just finds a position for it to go.
     * Expended this function so that the newly inserted particles are above the old. This stops new small beinging inserted in the old packing which can break the homogenious nature and if the case of segregation initial condition even place particles to low.
     */
    bool particleInsertable(BaseParticle& P0)
    {
        //Set up some local need data; counter for the number of failed attempted and a vector for the position.
        int failCounter=0;
        Vec3D position;
        
        // Choose a random particle position and zero intial velocity
        do
        {
            position.X=random.getRandomNumber(getXMin(),getXMax());
            position.Y=random.getRandomNumber(getYMin(),getYMax());
            position.Z=random.getRandomNumber(zMinInsert,getZMax());
            P0.setPosition(position);
            P0.setVelocity(Vec3D(0.0, 0.0, 0.0));
            
            failCounter++;
            
            // If you failed to find a locations 1000 times. Give up and let the particles settle out of the way.
            if (failCounter==1000) return false;
            
            
        }
        // If there is no particle in the way we are good to go
        while (!checkParticleForInteraction(P0));
        
        return true;
        
    }
    
    /*!
     * \brief Set the lambda functions so that the left and right walls do move
     * \details This uses lamdba functions to move the wall but it is a bit of a trick and requires the passaing of point to the wall; because the normal has to be changes and I saw no other easy way to do that with the current interface. There must be one and this should be revisited later.
     */
    void setupMovingWalls()
    {
        left->setPrescribedPosition([this](double time)
                                    {
                                        //If you are in step 5 start to rotate the angle of the wall sinsodally. note took a while to get this formaular right.
                                        if (step==5)
                                        {
                                            double angle=shearBoxAngle*std::sin(2.0*constants::pi*time/shearBoxTimePeriod);
                                            
                                            left->setNormal(Vec3D(-1.0*std::cos(angle),0.0,std::sin(angle)));
                                        }
                                        return(Vec3D(getXMin(),0.0,0.0));
                                        
                                    });
        
        
        right->setPrescribedPosition([this](double time)
                                     {
                                         //As above the anlge is minus the other wall as the wall is the other way round.
                                         if (step==5)
                                         {
                                             double angle=shearBoxAngle*std::sin(2.0*constants::pi*time/shearBoxTimePeriod);
                                             
                                             right->setNormal(Vec3D(std::cos(angle),0.0,-1.0*std::sin(angle)));
                                         }
                                         return(Vec3D(getXMax(),0.0,0.0));
                                         
                                     });
        
       // right->set(Vec3D(1.0, 0.0, 0.0), Vec3D(getXMax(),0.0,0.0));
        
        
    }
    
    
    /*!
     * \brief This is the work horse of initially filling the shear box. It deals with step 2-4
     * \details This code is called after everytime step and deals with the following parts of setting up the problem
     * Step 2 : Fill particles: insert particles (with step 3)
     * Step 3 : Let the particles settle so more can be inserted not return to step 2
     * Note, if you are doing normally graded insertion step 2 and 3 will be done seperate for large and small particles.
     * Step 4 : This is a full settle until the particles are almost stationary
     */
    void actionsBeforeTimeStep() override
    {
        
        //If step 2; just call createParticles which does the insert check and will put us on either step 3 or step 4 depending if all the required particles have been inserted.
        if (step==2)
            createParticles();
        
        // We are settling the particles either to get new space for inserting (step 3) or to make the particles stationary step 4.
        if (step==4 || step==3)
        {
            //   std::cout << "check"  << std::endl;
            //If it has been a while seen the last time we checked the KE>
            if (getTime() > checkTime)
            {
                
                //Report out the current KE.
                std::cout << "Current KE " << getKineticEnergy() << std::endl;
                if (getKineticEnergy() < (particleHandler.getNumberOfObjects()/100.0))
                {
                    if (step==4)
                    {
                        checkTime=getTime()+1.0;
                        if (getKineticEnergy() < (particleHandler.getNumberOfObjects()/1000.0))
                        {
                            //Now everything is settled move on to step 5
                            step=5;
                            //Make the current time minus one time step; so everything start next time step at 0
                            setTime(-getTimeStep());
                            
                            //Force a save now
                            setLastSavedTimeStep(NEVER);
                            
                            //Set up the number of saves. The users sets the number of saves per time period
                            /// \bug This does seem to drift slightly. But the time shots in teh file seem right so it may be a rounding issues but some more investigation is required
                            int numberOfSaves=round(getTimeMax()/shearBoxTimePeriod*saveCountPerCycle);
                            setSaveCount(helpers::getSaveCountFromNumberOfSavesAndTimeMaxAndTimeStep(numberOfSaves+1, getTimeMax(), getTimeStep())-1);
                            
                            //Report we have moved on to step 5.
                            std::cout << "\n \n \n";
                            std::cout << "STEP 5: Move the wall " << std::endl;
                            std::cout << "---------------------------" << std::endl;
                            std::cout << "\n \n \n";
                            // Why?
                            checkTime=getTime()+1.0;
                            
                            
                            
                            
                        }
                    }
                    else
                    {
                        
                        //Calc the max height of the settled particles for the new insert step and then insert again. So no new particle can be inside the previous inserted particles.
                        for (int i = 0; i < particleHandler.getNumberOfObjects(); i++)
                        {
                            BaseParticle* P0 = particleHandler.getObject(i);
                            //if (P0->getIndSpecies() == 0)
                            // {
                                zMinInsert = std::max(zMinInsert,P0->getPosition().getComponent(2)+radius_l+1);
                            // }
                        }
                        
                        
                        //Go back to step 2 after step 3.
                        step=2;
                    }
                }
                else
                {
                    checkTime=getTime()+1.0;
                }
            }
        }
        
        
        
    }
    
    
    /*!
     * \brief sets the parameters which control the shear box
     * \param[in] angleDegree the maximum angle the shear box is tited to (in degrees)
     * \param[in] timePeriod the time period of the oscillations
     * \details This is a single set function. Note, internally the angle is stored in radians
     */
    void setShearParameters(double angleDegree, double timePeriod)
    {
        
        shearBoxAngle=angleDegree/360.0*2.0*constants::pi;
        shearBoxTimePeriod=timePeriod;
        
    }
    
    /*!
     * \brief set the size ratio between the large and small particles in the sytem
     * \param[in] sizeRatio a double which stores the ratio of the sizes of the particles
     * \details This is a set function for the size of the particles. Due to the implicit non-dim begin used the small particles have a radius of 0.5 and the large particles 0.5%sizeRatio.
     */
    void setSizeRatio(double sizeRatio)
    {
        radius_s=0.5;
        radius_l=radius_s*sizeRatio;
    }
    
    /*!
     * \brief Function to set the initial conditions to normally graded instead of the default homogenously mixted
     * \details This sets the flag segregatedFilled to true; which defaults to false. Onces set there is no why to unset this flad at the moment. So not calling this functions is the alternative.
     */
    void makeInitialConditionsNormallyGraded()
    {
        segregatedFilling=true;
    }
    
    /*!
     * \brief Functions to turn of the saving of data during the filling process
     * \details This sets the flag ShowFillingFlag to true; which defailts to false. Onces set there is no why to unset this flag. If this function is called saves to the data file will be written every 0.1 dimensionless sections while the drum is filling. If it is not set; only the first time step is writen and then it does not start to write data again until the shear box starts to shear i.e. step 5.
     */
    void showFilling()
    {
        showFillingFlag=true;
        
    }
    
    /*!
     * \brief Set the number of small and large particles in the simulation
     * \param[in] small an integer indicating the number of small particles to be inserted into the shear box
     * \param[out] large an integer indicating the number of large particles to be inserted into the shear box
     * \details Set, how many particles of each size and to be used in the simulations. Note, you do not have to worry if they fit as the height of the box will be automatically resized so that to accomodate all the required particles
     * \todo function has no error checking
     */
    void setNumberOfParticles(int small, int large)
    {
        numSmallToBeInserted=small;
        numLargeToBeInserted=large;
        
        volfrac =  (small*(4/3)*constants::pi*pow(radius_s,3)) / ( (large*(4/3)*constants::pi*pow(radius_l,3)) +         (small*(4/3)*constants::pi*pow(radius_s,3)) );
        
        std::cout << "Small partile volume fraction is " << volfrac <<std::endl;
        
        
    }
    
    /*!
     * \brief Set the number of saves (snapshots) taken per cycle of the shear box
     * \param[in] count unsigned int which is the number of saves for each cycle of the box must be one or greated
     * \details This is the number of saves (snapshots) per oscillation. So if you set it to one you should always see the shear box unright.
     * \bug This does seem to drift slightly due to a rounding error this should be looked into
     * \todo This current has no error checking that the number is set greater than one.
     */
    void setSaveCountPerCycle(unsigned int count)
    {
        saveCountPerCycle=count;
        
    }
    
    /*!
     * \brief set the horzontial dimensions to the box
     * \param[in] length double which is the size of the box in the x-direction
     * \param[in] width double which is the size of the box in the y-direction
     * \details Set the horzontial dimension of the box to be from (0,length) and (0,width) in the x and y directions respectivly.
     * Note the size of the box in the z-direction (height, agaist gravity) is automatically calculated based on this sizes and the number of partilces trying to be inserted.
     * \todo Error checking needs adding
     */
    void setBoxDimensions(double length, double width)
    {
        setXMin(0.0);
        setXMax(length);
        
        setYMin(0.0);
        setYMax(width);
        
    }
    
    void setNumSmallAtBottom(int _numSmallAtBottom)
    {
        numSmallAtBottom   =_numSmallAtBottom;
        numSmallToBeInserted=numSmallToBeInserted+_numSmallAtBottom;
    }
    
    /// The radius of the small particles
    double radius_s;
    /// The radius of the smaller small particles
    double radius_s2;
    /// The radius of the large particles
    double radius_l;
    /// The radius of the smaller large particle
    double radius_l2;
    /// small particle volume fraction
    double volfrac;
    
    /// Store the current step of the problem
    unsigned int step;
    
    /// Store the current step of the problem
    unsigned int numSmallAtBottom;
    
    /// Internal varible used in set 3 and 4 to see what to check if the KE is lower enough again
    double checkTime;
    /// Internal varibles used to check when the top of the flow is after one set of inserting. Used to make sure the later sets in inserted particles are above the previous
    double zMinInsert;
    
    /// The maximum angle of the shear box in radians
    double shearBoxAngle;
    /// The time period of oscilation of the shear box
    double shearBoxTimePeriod;
    
    /// The number of small particles the user requested to be inserted
    int numSmall;
    /// The number of large particles the uesr requested to be inserted
    int numLarge;
    /// The number of small particles which still have to be inserted. Note, starts at numSmall and falls to zero during steps 2 and 3.
    int numSmallToBeInserted;
    /// The number of large particles which still have to be inserted. Note, starts at numLarge and falls to zero during steps 2 and 3.
    int numLargeToBeInserted;
    
    /// The number of save (snapshots) taken per oscillation of the box
    unsigned int saveCountPerCycle;
    /// If true the box is normally graded initially; if true it is homogenious.
    bool segregatedFilling;
    /// If true the filling steps and saved to disk; if false only the first frame is
    bool showFillingFlag;
    
    /// Pointer to the left (moving) wall
    InfiniteWall* left;
    ///Pointer to the right (moving) wall
    InfiniteWall* right;
    
    /// Pointer to the one species in the problem; hence the name,
    LinearViscoelasticFrictionSpecies* speciesS1;
   LinearViscoelasticFrictionSpecies* speciesS2;
    LinearViscoelasticFrictionSpecies* speciesL1;
    LinearViscoelasticFrictionSpecies* speciesL2;
   LinearViscoelasticFrictionSpecies* speciesWall;
    
};

int main(int argc, char** argv)
{
    //Create an instance of the problem.
    ShearBox problem;
    
    
    //If you are restarted you have an extra flag and in this case append not overwrite the files.
    if (argc>1)
    {
        problem.setAppend(true);
    }
    
    //Set Xballs to colour on size and in the solid mode with no velocity vectors.
    problem.setXBallsAdditionalArguments("-cmode 7");
    
    //Sim parameters all should be clear from the access functions.
    problem.setTimeMax(1000000.0);
    problem.setBoxDimensions(13,6.7);
    problem.setShearParameters(15.0,100);  //angle, period
    problem.setSizeRatio(2.1);
    problem.makeInitialConditionsNormallyGraded();
    problem.setNumberOfParticles(950,105);
    problem.setNumSmallAtBottom(40);

    
    problem.setSaveCountPerCycle(1);
//    problem.showFilling();
    
    //problem.autoNumber();
    problem.setName("ShearBox");
    
    
    //Now, start the simulation
    problem.solve(argc, argv);
    return 0;
}
