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

//updated force calculation to the changes in MD.cc (further, a few variables needed to be privatized)
#include <sstream>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <cmath>

#include "Chute.h"
#include "Boundaries/PeriodicBoundary.h"
#include "sys/stat.h"
///Particles of a single Species
class ChuteWithPeriodicInflow : public Chute
{
public:
    ChuteWithPeriodicInflow(std::string restart_file)
    {
        loadPeriodicBox(restart_file);
    }
    
    ChuteWithPeriodicInflow(std::string restart_file, int numRepetitions)
    {
        loadPeriodicBox(restart_file);
        AddContinuingBottom(numRepetitions);
    }
    
    ChuteWithPeriodicInflow(std::string restart_file, int numRepetitions, int numRepetitionsInWidth)
    {
        loadPeriodicBox(restart_file);
        AddContinuingBottom(numRepetitions);
        ExtendInWidth(numRepetitionsInWidth);
    }
    
    double getInfo(const BaseParticle& P) const
            {
        return P.getIndSpecies();
    }
    
    ///loads periodic chute data from restart file
    void loadPeriodicBox(std::string const restart_file)
    {
        FillChute = false;

        std::cout << "constructor" << std::endl;

        //load particles and chute setup into periodic chute
        setName(restart_file.c_str());
        readRestartFile();

		set_PeriodicBoxLength(getXMax());
		set_PeriodicBoxNSpecies(speciesHandler.getNumberOfObjects());
		
		//we don't want to treat the data as restarted, i.e. we start the program at time=0 and create new output files
		setRestarted(false);

        //keep file name but create files in the local directory, i.e. remove folder
        size_t found = restart_file.find_last_of("/\\");
        setName(restart_file.substr(found + 1).c_str());

        // adds new species for the flow particles (0 for the periodic inflow particles)
        for (int i = 0; i < get_PeriodicBoxNSpecies(); i++)
            addSpecies(Species[i]);

        setName("ChuteWithPeriodicInflow");
    }
    
    void AddContinuingBottom(int numRepetitions)
    {
        // creates bottom outside periodic chute of species 1
        double lengthPeriodicChute = getXMax();
        int N = particleHandler.getNumberOfObjects();
        particleHandler.setStorageCapacity(N * (2. + numRepetitions));
        //allows for a gap between the infow and flow regimes
        double Gap = 0.;

        for (int j = 1; j < numRepetitions; j++)
        {
            for (int i = 0; i < N; i++)
            {
                if (particleHandler.getObject(i)->isFixed() | FillChute)
                {
                    particleHandler.addObject(particleHandler.getObject(i)->copy());
                    particleHandler.getLastObject()->setSpecies(speciesHandler.getObject(1));
                    particleHandler.getLastObject()->move(Vec3D(Gap + lengthPeriodicChute * j, 0., 0.));
                }
            }
        }

        setXMax(numRepetitions * lengthPeriodicChute);
        //setHGridNumberOfBucketsToPower(particleHandler.getStorageCapacity());
    }

    void ExtendInWidth(int numRepetitionsInWidth)
    {
        PeriodicBoundary* perw = static_cast<PeriodicBoundary*>(boundaryHandler.getLastObject());
        // creates bottom outside periodic chute of species 1
        double widthPeriodicChute = getYMax();
        int N = particleHandler.getNumberOfObjects();
        particleHandler.setStorageCapacity(N * numRepetitionsInWidth);
        for (int j = 1; j < numRepetitionsInWidth; j++)
        {
            for (int i = 0; i < N; i++)
            {
                //if (particleHandler.getObject(i)->isFixed()) {
                particleHandler.addObject(particleHandler.getObject(i)->copy());
                //Particles.back()->get_IndSpecies()=1;
                particleHandler.getLastObject()->move(Vec3D(0.0, widthPeriodicChute * j, 0.0));
                //}
            }
        }
        setYMax(numRepetitionsInWidth * widthPeriodicChute);
        perw->set(Vec3D(0, 1, 0), getYMin(), getYMax());
        //setHGridNumberOfBucketsToPower(particleHandler.getNumberOfObjects());
    }

    ///Do not add, only remove particles
    void actionsBeforeTimeStep()
    {
        cleanChute();
    }

    ///Remove particles if they fall below a certain height (allows them to become supercritical)
    void cleanChute()
    {
        //clean outflow every 100 time steps
        static int count = 0, maxcount = 100;
        if (count > maxcount)
        {
            count = 0;
            // delete all outflowing particles
            for (unsigned int i = 0; i < particleHandler.getNumberOfObjects();)
            {
                //~ if (particleHandler.getObject(i)->getPosition().Z<getZMin()-10*getInflowParticleRadius()){
                if (particleHandler.getObject(i)->getPosition().X > getXMax())
                {
                    //~ cout << "Remove particle" << endl;
                    particleHandler.removeObject(i);
                }
                else
                    i++;
            }
        }
        else
            count++;
    }
    
    ///Do not add bottom
    void setupInitialConditions()
    {
    }
    
    void integrateBeforeForceComputation()
    {
        for (std::vector<BaseParticle*>::iterator it = particleHandler.begin(); it != particleHandler.end(); ++it)
        {
            (*it)->integrateBeforeForceComputation(getTimeStep());

            // This shifts particles that moved through periodic walls
            PeriodicBoundary* perw = static_cast<PeriodicBoundary*>(boundaryHandler.getObject(0));
            PeriodicBoundary* perw1 = static_cast<PeriodicBoundary*>(boundaryHandler.getObject(1));
            if ((*it)->getIndSpecies() == 0 && perw->getDistance(**it) < 0)
            {
                //if (iP->getIndSpecies()==0 && static_cast<PeriodicWall*>(boundaryHandler.getObject(0))->getDistance(*iP)<0) {
                if (!perw->isClosestToLeftBoundary())
                {
                    if (particleHandler.getStorageCapacity() <= particleHandler.getNumberOfObjects())
                        particleHandler.setStorageCapacity(particleHandler.getStorageCapacity() + 10000);
                    particleHandler.addObject((*it)->copy());
                    const ParticleSpecies* species = speciesHandler.getObject(particleHandler.getLastObject()->getIndSpecies() + get_PeriodicBoxNSpecies());
                    particleHandler.getLastObject()->setSpecies(species);
                    perw->shiftPosition((*it));
                    if (!getHGridUpdateEachTimeStep())
                    {
                        hGridRemoveParticle((*it));
                        hGridUpdateParticle((*it));
                    }
                }
                else
                {
                    static int count = -1;
                    count++;
                    if (!(count % dataFile.getSaveCount()))

                        std::cout << "Warning: Particle " << (*it)->getIndex() << " is left of xmin: x="
                                << (*it)->getPosition().X << ", v_x=" << (*it)->getVelocity().X << std::endl;
                }
            }
            if (getIsPeriodic() > 1 && perw1->getDistance(**it) < 0)
            {
                perw1->shiftPosition(*it);
                if (!getHGridUpdateEachTimeStep())
                {
                    hGridRemoveParticle(*it);
                    hGridUpdateParticle(*it);
                }
            }
        }

    }
    
    ///add some particular output
    void printTime() const
    {
        int counter = 0;
        
        double speed1 = 0;
        double speed0 = 0;
        double n1 = 0;
        double n0 = 0;
        for (int i = 0; i < particleHandler.getNumberOfObjects(); i++)
        {
            if (particleHandler.getObject(i)->getIndSpecies() == 0)
                counter++;
            if (!particleHandler.getObject(i)->isFixed())
            {
                if (particleHandler.getObject(i)->getIndSpecies() == 0)
                {
                    speed0 += particleHandler.getObject(i)->getVelocity().X;
                    n0++;
                }
                else
                {
                    speed1 += particleHandler.getObject(i)->getVelocity().X;
                    n1++;
                }
            }
        }
        speed0 /= n0;
        speed1 /= n1;
        
        std::cout << "t=" << std::setprecision(5) << std::left << std::setw(6) << getTime()
                << ", n=" << n1 << "(inflow:" << n0 << ")"
                << ", speed= " << speed1 << "(" << speed0 << ")"
                << std::endl;
        //~ cout
        //~ << " A " << PeriodicBoxLength
        //~ << " A " << PeriodicBoxNSpecies
        //~ << " A " << FillChute << endl;
    }
    //%%%%%%%%%%%%%%%%% Copied from MD.cc and then edited
    void computeInternalForces(BaseParticle* P1, BaseParticle* P2)
    {
<<<<<<< .mine
        // For particles of the same species, set species vector to Species(PI);
        // for particles of different species, set species vector to MixedSpecies(PI,PJ)
        BaseSpecies* pSpecies = (PI->getIndSpecies() == PJ->getIndSpecies()) ? &Species[PI->getIndSpecies()] : speciesHandler.getMixedObject(PI->getIndSpecies(), PJ->getIndSpecies());
=======
        //this is because the rough bottom allows overlapping fixed particles
        if (P1->isFixed() && P2->isFixed())
            return;
>>>>>>> .r824

        ///Tangential spring information is always store in the real particle with highest index
        ///When a Periodic contact is encountered it is always encoutered twice (for normal periodic walls), but the force is only applied to the real particle
        ///The tangential spring information for periodic particles is stored in the normal particle (and thus twice for normal periodic interactions)
        ///When a Particle is removed the tangential spring information has to be moved
        
        //Dificult cases for tangential springs in comination with periodic walls are:
        //
        // A new contact over a periodic wall:
        // Starting situation: There are no tangential springs stored at all
        // Creating periodic particles: There are no tangential springs so nothing happens
        // Contact detection: There are two contacts detected, one (CA) for P1 with P2per and one (CB) for P2 with P1per
        // Switch to the 4 cases:
        //   CA: PI=P1, PJ=P2per PJreal=P2
        //   CB: PI=P2, PJ=P1per PJreal=P1
        // Reconnecting springs:
        //   CA: PI=P1 and PJ=P2per do not have a spring, so a new spring is created in either PI or PJ (could be in periodic particle)
        //   CB: PI=P2 and PJ=P1Per do not have a spring, so a new spring is created in either PI or PJ (could be in periodic particle)
        // Removing periodic particles: One of the springs will be stored in a periodic particle and thus removed, the other spring is kept and used (this is the real particle with the kighest index)
        
        // Reconnect to a contact over a periodic wall:
        // Starting situation: There is a tangential spring stored in the particle with the highest index, lets assume this is P1
        // Creating periodic particles: P1 has a tangential spring, so P1per has a reversed copy of this spring
        // Contact detection: There are two contacts detected, one (CA) for P1 with P2per and one (CB) for P2 with P1per
        // Switch to the 4 cases:
        //   CA: PI=P1, PJ=P2per PJreal=P2
        //   CB: PI=P2, PJ=P1per PJreal=P1
        // Reconnecting springs:
        //   CA: PI=P1 and PJ=P2per have a spring in P1, so we reconnect to this spring in the normal way and integrate it.
        //   CB: PI=P2 and PJ=P1Per have a spring in P1per, however this spring has to be a reversed spring since it is in PJ.
        // Removing periodic particles: The spring in P1per is removed
        
        //Cases (start from P1 and P2 and go to PI and PJ
        //1 Normal-Normal		->PI=max(P1,P2), PJ=min(P1,P2)
        //2 Periodic-Normal		->(PI=P2 PJ=real(P1))
        //3 Normal-Periodic		->(PI=P1 PJ=real(P2))
        //4 Periodic-Periodic	->do nothing
        
        //Just some statements to handle the 4 cases
        BaseParticle *PI, *PJ, *PJreal;
        bool isPeriodic;
        BaseParticle *P1Per = P1->getPeriodicFromParticle();
        BaseParticle *P2Per = P2->getPeriodicFromParticle();
        if (P1Per == NULL)
        {
            if (P2Per == NULL) //N-N
                if (P1->getIndex() < P2->getIndex())
                {
                    PI = P2;
                    PJ = P1;
                    PJreal = P1;
                    isPeriodic = false;
                }
                else
                {
                    PI = P1;
                    PJ = P2;
                    PJreal = P2;
                    isPeriodic = false;
                }
            else //N-P
            {
                PI = P1;
                PJ = P2;
                PJreal = P2Per;
                isPeriodic = true;
            }
        }
        else
        {
            if (P2Per == NULL) //P-N
            {
                PI = P2;
                PJ = P1;
                PJreal = P1Per;
                isPeriodic = true;
            }
            else
                return;
        }

#ifdef DEBUG_OUTPUT
        std::cerr << "In computing internal forces between particle "<<PI->getPosition()<<" and "<<PJ->getPosition()<<std::endl;
#endif

        //Get the square of the distance between particle i and particle j
        Mdouble dist_squared = Vec3D::getDistanceSquared(PI->getPosition(), PJ->getPosition());
        Mdouble interactionRadii_sum = PI->getInteractionRadius() + PJ->getInteractionRadius();

#ifdef DEBUG_OUTPUT_FULL
        std::cerr << "Square of distance between " << dist_squared << " square sum of radii " << radii_sum*radii_sum <<std::endl;
#endif

        // True if the particles are in contact
        if (dist_squared < (interactionRadii_sum * interactionRadii_sum))
        {
            // For particles of the same species, set species vector to Species(PI);
            // for particles of different species, set species vector to MixedSpecies(PI,PJ)
            CSpecies* pSpecies = (PI->getIndSpecies() == PJ->getIndSpecies()) ? &Species[PI->getIndSpecies()] : getMixedSpecies(PI->getIndSpecies(), PJ->getIndSpecies());
            
            // Calculate distance between the particles
            Mdouble dist = sqrt(dist_squared);

            // Compute normal vector

            Vec3D normal = (PI->getPosition() - PJ->getPosition()) / dist;

            // Compute the overlap between the particles
            Mdouble radii_sum = PI->getRadius() + PJ->getRadius();
            Mdouble deltan = std::max(0.0, radii_sum - dist);
            
            Vec3D force = Vec3D(0, 0, 0);
            ;
            Vec3D forcet;
            forcet.setZero();
            Vec3D forcerolling;
            forcerolling.setZero();
            Vec3D forcetorsion;
            forcetorsion.setZero();
            Mdouble fdotn = 0;
            CTangentialSpring* TangentialSpring = NULL;
            
            //evaluate shortrange non-contact forces
            if (pSpecies->getAdhesionForceType() != AdhesionForceType::NONE)
                fdotn += computeShortRangeForceWithParticle(PI, PJ, PJreal, pSpecies, dist);
            
            if (deltan > 0) //if contact forces
            {
                
                // Compute the relative velocity vector v_ij
                Vec3D vrel;
                if (!pSpecies->getSlidingFrictionCoefficient())
                {
                    vrel = (PI->getVelocity() - PJ->getVelocity());
                }
                else
                {
                    vrel = (PI->getVelocity() - PJ->getVelocity()) + Vec3D::cross(normal, PI->getAngularVelocity() * (PI->getRadius() - .5 * deltan) + PJ->getAngularVelocity() * (PJ->getRadius() - .5 * deltan));
                }

                // Compute the projection of vrel onto the normal (can be negative)
                Mdouble vdotn = -Vec3D::dot(vrel, normal);

                //update restitution coeff
                if (pSpecies->getIsRestitutionCoefficientConstant())
                    pSpecies->updateDissipation(PI->getMass(), PJ->getMass());
                Mdouble a = 0, R = 0;

                // Compute normal force on particle i due to contact
                if (pSpecies->getForceType() == ForceType::HERTZ_MINDLIN || pSpecies->getForceType() == ForceType::HERTZ_MINDLIN_DERESIEWICZ)
                {
                    //R is twice the effective radius
                    R = 2.0 * PI->getRadius() * PJ->getRadius() / (PI->getRadius() + PJ->getRadius());
                    a = sqrt(R * deltan);
                    //pSpecies->getStiffness() stores the elastic modulus
                    Mdouble kn = 4. / 3. * pSpecies->getStiffness() * a;
                    fdotn += kn * deltan + pSpecies->getDissipation() * vdotn;
                }
                else
                {
                    fdotn += pSpecies->getStiffness() * deltan + pSpecies->getDissipation() * vdotn;
                }
                force += normal * fdotn;
                
                //If tangential forces are present
                if (pSpecies->getSlidingFrictionCoefficient() || pSpecies->getRollingFrictionCoefficient() || pSpecies->getTorsionFrictionCoefficient())
                {
                    //call tangential spring
                    if (pSpecies->getSlidingStiffness() || pSpecies->getRollingStiffness() || pSpecies->getTorsionStiffness())
                        TangentialSpring = getTangentialSpring(PI, PJ, PJreal);

                    //Compute norm of normal force
                    Mdouble norm_fn = fabs(fdotn);
                    
                    //calculate sliding friction
                    if (pSpecies->getSlidingFrictionCoefficient())
                    {
                        //Compute the tangential component of vrel
                        Vec3D vrelt = vrel + normal * vdotn;
                        //Compute norm of vrelt
                        Mdouble vdott = vrelt.getLength();
                        
                        if (pSpecies->getSlidingStiffness())
                        {
                            Vec3D* delta = &(TangentialSpring->delta);

                            //Integrate the spring
                            ///Both options are up to first order the same (the first one is nicer becaus it always keeps the spring tangential, whereas the second one is in a nicer intergration form)
                            //(*delta) += vrelt * dt - Vec3D::Dot(*delta,normal)*normal;
                            Vec3D ddelta = (vrelt - Vec3D::dot(*delta, PI->getVelocity() - PJ->getVelocity()) * normal / dist) * getTimeStep();
                            (*delta) += ddelta;

                            //Calculate test force including viscous force
                            if (pSpecies->getForceType() == ForceType::HERTZ_MINDLIN)
                            {
                                //pSpecies->getSlidingStiffness() stores the elastic shear modulus
                                Mdouble kt = 8. * pSpecies->getSlidingStiffness() * a;
                                TangentialSpring->SlidingForce += -kt * ddelta;
                                forcet = TangentialSpring->SlidingForce - pSpecies->getSlidingDissipation() * vrelt;
                            }
                            else if (pSpecies->getForceType() == ForceType::HERTZ_MINDLIN_DERESIEWICZ)
                            {
                                //pSpecies->getSlidingStiffness() stores the elastic shear modulus
                                forcet = TangentialSpring->SlidingForce - pSpecies->getSlidingDissipation() * vrelt;
                                Mdouble kt = 8. * pSpecies->getSlidingStiffness() * a * std::pow(1 - Vec3D::getLength(forcet) / pSpecies->getSlidingFrictionCoefficient() / fdotn, 0.33);
                                TangentialSpring->SlidingForce += -kt * ddelta;
                                forcet = TangentialSpring->SlidingForce - pSpecies->getSlidingDissipation() * vrelt;
                            }
                            else
                            {
                                forcet = (-pSpecies->getSlidingDissipation()) * vrelt - pSpecies->getSlidingStiffness() * (*delta);
                            }
                            Mdouble forcet2 = forcet.getLengthSquared();

                            //tangential forces are modelled by a spring-damper of elastisity kt and viscosity getSlidingDissipation() (sticking),
                            //but the force is limited by Coulomb friction (sliding):
                            //f_t = -getSlidingDissipation()*vrelt, if getSlidingDissipation()*vrelt<=mu_s*fdotn, f_t=mu+s*fdotn*t, else
                            double muact = (TangentialSpring->sliding) ? (pSpecies->getSlidingFrictionCoefficient()) : (pSpecies->getSlidingFrictionCoefficientStatic()); // select mu from previous sliding mode
                            if (forcet2 <= mathsFunc::square(muact * norm_fn))
                            {
                                //sticking++;
                                TangentialSpring->sliding = false;
                            }
                            else
                            {
                                //sliding++;
                                TangentialSpring->sliding = true;
                                ///\todo{The spring should be cut back such that fdott=mu*fdotn. This is simple for getSlidingDissipation()=0; we have to think about what happens in the sliding case with tang. dissipation; same for walls; Update Dec-2-2011: fixed}
                                Mdouble norm_forcet = sqrt(forcet2);
                                forcet *= pSpecies->getSlidingFrictionCoefficient() * norm_fn / norm_forcet;
                                TangentialSpring->SlidingForce = forcet + pSpecies->getSlidingDissipation() * vrelt;
                                (*delta) = TangentialSpring->SlidingForce / (-pSpecies->getSlidingStiffness());
                            }

                            //Add tangential force to total force
                            force += forcet;

                        }
                        else
                        { //if no tangential spring
                          //tangential forces are modelled by a damper of viscosity getSlidingDissipation() (sticking),
                          //but the force is limited by Coulomb friction (sliding):
                          //f_t = -getSlidingDissipation()*vrelt, if getSlidingDissipation()*vrelt<=mu_s*fdotn, f_t=mu+s*fdotn*t, else
                            if (vdott * pSpecies->getSlidingDissipation() <= pSpecies->getSlidingFrictionCoefficientStatic() * norm_fn)
                            { //sticking++;
                                forcet = -pSpecies->getSlidingDissipation() * vrelt;
                            }
                            else
                            { //sliding++;
                              //set force to Coulomb limit
                                forcet = -(pSpecies->getSlidingFrictionCoefficient() * norm_fn / vdott) * vrelt;
                            }
                            //Add tangential force to total force
                            force += forcet;
                        }
                    }

                    //calculate rolling friction
                    if (pSpecies->getRollingFrictionCoefficient())
                    {
                        //From Luding2008, objective rolling velocity (eq 15) w/o 2.0!
                        Mdouble reducedRadiusI = PI->getRadius() - .5 * deltan;
                        Mdouble reducedRadiusJ = PJ->getRadius() - .5 * deltan;
                        Mdouble reducedRadiusIJ = 2.0 * reducedRadiusI * reducedRadiusJ / (reducedRadiusI + reducedRadiusJ);
                        Vec3D vrolling = -reducedRadiusIJ * Vec3D::cross(normal, PI->getAngularVelocity() - PJ->getAngularVelocity());
                        if (pSpecies->getRollingStiffness())
                        {
                            Vec3D* RollingSpring = &(TangentialSpring->RollingSpring);

                            //Integrate the spring
                            (*RollingSpring) += vrolling * getTimeStep(); // - Vec3D::Dot(*RollingSpring,normal)*normal;

                                    //Calculate test force including viscous force
                            forcerolling = (-pSpecies->getRollingDissipation()) * vrolling - pSpecies->getRollingStiffness() * (*RollingSpring);
                            Mdouble forcerolling2 = forcerolling.getLengthSquared();

                            //tangential forces are modelled by a spring-damper of elastisity kt and viscosity getSlidingDissipation() (sticking),
                            //but the force is limited by Coulomb friction (sliding):
                            //f_t = -getSlidingDissipation()*vrelt, if getSlidingDissipation()*vrelt<=mu_s*fdotn, f_t=mu+s*fdotn*t, else
                            double muact = (TangentialSpring->slidingRolling) ? (pSpecies->getRollingFrictionCoefficient()) : (pSpecies->getRollingFrictionCoefficientStatic()); // select mu from previous sliding mode
                            if (forcerolling2 <= mathsFunc::square(muact * norm_fn))
                            {
                                //sticking++;
                                TangentialSpring->slidingRolling = false;
                            }
                            else
                            {
                                //sliding++;
                                TangentialSpring->slidingRolling = true;
                                forcerolling *= pSpecies->getRollingFrictionCoefficient() * norm_fn / sqrt(forcerolling2);
                                (*RollingSpring) = (forcerolling + pSpecies->getRollingDissipation() * vrolling) / (-pSpecies->getRollingStiffness());
                            }

                            //Add tangential force to torque
                            Vec3D Torque = reducedRadiusIJ * Vec3D::cross(normal, forcerolling);
                            PI->addTorque(Torque);
                            PJreal->addTorque(-Torque);
                        }
                    }

                    //calculate torsive friction
                    if (pSpecies->getTorsionFrictionCoefficient())
                    {
                        //From Luding2008, spin velocity (eq 16) w/o 2.0!
                        Mdouble RadiusIJ = 2.0 * PI->getRadius() * PJ->getRadius() / (PI->getRadius() + PJ->getRadius());
                        Vec3D vtorsion = RadiusIJ * Vec3D::dot(normal, PI->getAngularVelocity() - PJ->getAngularVelocity()) * normal;
                        if (pSpecies->getTorsionStiffness())
                        {
                            //~ std::cout << "Error; not yet implemented" << std::endl;
                            //~ exit(-1);
                            Vec3D* TorsionSpring = &(TangentialSpring->TorsionSpring);

                            //Integrate the spring
                            (*TorsionSpring) = Vec3D::dot((*TorsionSpring) + vtorsion * getTimeStep(), normal) * normal;

                            //Calculate test force including viscous force
                            forcetorsion = (-pSpecies->getTorsionDissipation()) * vtorsion - pSpecies->getTorsionStiffness() * (*TorsionSpring);
                            Mdouble forcetorsion2 = forcetorsion.getLengthSquared();

                            //tangential forces are modelled by a spring-damper of elastisity kt and viscosity getSlidingDissipation() (sticking),
                            //but the force is limited by Coulomb friction (sliding):
                            //f_t = -getSlidingDissipation()*vrelt, if getSlidingDissipation()*vrelt<=mu_s*fdotn, f_t=mu+s*fdotn*t, else
                            double muact = (TangentialSpring->slidingTorsion) ? (pSpecies->getTorsionFrictionCoefficient()) : (pSpecies->getTorsionFrictionCoefficientStatic()); // select mu from previous sliding mode
                            if (forcetorsion2 <= mathsFunc::square(muact * norm_fn))
                            {
                                //sticking++;
                                TangentialSpring->slidingTorsion = false;
                            }
                            else
                            {
                                //sliding++;
                                TangentialSpring->slidingTorsion = true;
                                //~ std::cout << "sliding " << std::endl;
                                forcetorsion *= pSpecies->getTorsionFrictionCoefficient() * norm_fn / sqrt(forcetorsion2);
                                (*TorsionSpring) = (forcetorsion + pSpecies->getTorsionDissipation() * vtorsion) / (-pSpecies->getTorsionStiffness());
                            }

                            //Add tangential force to torque
                            Vec3D torque = RadiusIJ * forcetorsion;
                            PI->addTorque(torque);
                            PJreal->addTorque(-torque);

                        }
                    }
                } //end if tangential forces

                ///\todo TW: the following 13 lines concern only the sliding spring and could be moved into the if statement above
                //Make force Hertzian (note: deltan is normalized by the normal distanct of two particles in contact, as in Silbert)
                if (pSpecies->getForceType() == ForceType::HERTZ)
                    force *= sqrt(deltan / (PI->getRadius() + PJ->getRadius()));

            }
            else
            { //end if contact forces
                force += fdotn * normal;
            }

            //Add forces to total force
            //PI->addForce(force);
            //if (!isPeriodic)
            //    PJreal->addForce(-force);

            // Add torque due to tangential forces: t = Vec3D::Cross(l,f), l=dist*Wall.normal
            //if (pSpecies->getSlidingFrictionCoefficient())
            //{
            //    Vec3D Vec3D::Cross = Vec3D::Cross(normal, force);
            //    PI->addTorque(-Vec3D::Cross * (PI->getRadius() - .5 * deltan));
            //    if (!isPeriodic)
            //        PJreal->addTorque(-Vec3D::Cross * (PJ->getRadius() - .5 * deltan));
            //}
            //BEGIN: CHANGE
            if (InPeriodicBox(PI) != InPeriodicBox(PJreal))
            {
                if (InPeriodicBox(PI))
                {
                    if (PI->getPosition().X > get_PeriodicBoxLength() - 4.0 * particleHandler.getLargestParticle()->getRadius())
                    {
                        //~ if (PJreal->isFixed()) cout << "PJreal" << endl;
                        PJreal->addForce(-force);
                    }
                }
                else
                {
                    if (PJreal->getPosition().X > get_PeriodicBoxLength() - 4.0 * particleHandler.getLargestParticle()->getRadius())
                    {
                        //~ if (Particles[PI].isFixed()) cout << "PI" << endl;
                        PI->addForce(force);
                    }
                }
            }
            else
            {
                PI->addForce(force);
                PJreal->addForce(-force);
            }
            //END: CHANGE

            // Add torque due to tangential forces: t = Vec3D::Cross(l,f), l=dist*Wall.normal
            //if (pSpecies->getSlidingFrictionCoefficient()) {
            //	Vec3D Vec3D::Cross = Vec3D::Cross(normal, force);
            //	PI    ->add_Torque(-Vec3D::Cross * (PI->getRadius() - .5 * deltan));
            //	PJreal->add_Torque(-Vec3D::Cross * (PJ->getRadius() - .5 * deltan));
            //}
            //BEGIN: CHANGE
            if (pSpecies->getSlidingFrictionCoefficient())
            {
                Vec3D cross = Vec3D::cross(normal, force);
                if (InPeriodicBox(PI) != InPeriodicBox(PJreal))
                {
                    if (InPeriodicBox(PI))
                    {
                        if (PI->getPosition().X > get_PeriodicBoxLength() - 4.0 * particleHandler.getLargestParticle()->getRadius())
                        {
                            PJreal->addTorque(-cross * (PJreal->getRadius() - .5 * deltan));
                        }
                    }
                    else
                    {
                        if (PJreal->getPosition().X > get_PeriodicBoxLength() - 4.0 * particleHandler.getLargestParticle()->getRadius())
                        {
                            PI->addTorque(-cross * (PI->getRadius() - .5 * deltan));
                        }
                    }
                }
                else
                {
                    PI->addTorque(-cross * (PI->getRadius() - .5 * deltan));
                    PJreal->addTorque(-cross * (PJreal->getRadius() - .5 * deltan));
                }
            }
            //END:CHANGE
            
            // output for ene and stat files:
            if (eneFile.getSaveCurrentTimeStep())
            {
                if (!isPeriodic)
                    addElasticEnergy(0.5 * (pSpecies->getStiffness() * mathsFunc::square(deltan) + (TangentialSpring ? (pSpecies->getSlidingStiffness() * TangentialSpring->delta.getLengthSquared() + pSpecies->getRollingStiffness() * TangentialSpring->RollingSpring.getLengthSquared() + pSpecies->getTorsionStiffness() * TangentialSpring->TorsionSpring.getLengthSquared()) : 0.0)));
                else
                    addElasticEnergy(0.25 * (pSpecies->getStiffness() * mathsFunc::square(deltan) + (TangentialSpring ? (pSpecies->getSlidingStiffness() * TangentialSpring->delta.getLengthSquared() + pSpecies->getRollingStiffness() * TangentialSpring->RollingSpring.getLengthSquared() + pSpecies->getTorsionStiffness() * TangentialSpring->TorsionSpring.getLengthSquared()) : 0.0)));
            }
            if (fStatFile.getSaveCurrentTimeStep() || statFile.getSaveCurrentTimeStep() || getDoCGAlways())
            {

                Mdouble fdott = forcet.getLength();
                Mdouble deltat_norm = TangentialSpring ? (-TangentialSpring->delta.getLength()) : 0.0;

                ///\todo Define the center this way or are radii involved? Or maybe just use middle of overlap region?
                ///Thomas: Yes, we should correct that for polydispersed problems; however, then we also have to correct it in StatisticsVector::gatherContactStatistics.
                Vec3D centre = 0.5 * (PI->getPosition() + PJ->getPosition());

                ///The second particle (i.e. the particle the force acts on)
                ///is always a flow particle
                ///\todo{Is it the first particle the force acts on?}

                if (!PI->isFixed())
                {
                    if (statFile.getSaveCurrentTimeStep() || getDoCGAlways())
                        gatherContactStatistics(PJreal->getIndex(), PI->getIndex(), centre, deltan, deltat_norm, fdotn, fdott, -normal, -(fdott ? forcet / fdott : forcet));
                    if (fStatFile.getSaveCurrentTimeStep())
                        fStatFile.getFstream() << getTime() << " " << PJreal->getIndex() << " " << PI->getIndex() << " " << centre << " " << radii_sum - dist << " " << deltat_norm << " " << fdotn << " " << fdott << " " << -normal << " " << -(fdott ? forcet / fdott : forcet) << std::endl;
                }
                if (!PJreal->isFixed() && !isPeriodic)
                {
                    if (statFile.getSaveCurrentTimeStep() || getDoCGAlways())
                        gatherContactStatistics(PI->getIndex(), PJreal->getIndex(), centre, deltan, deltat_norm, fdotn, fdott, normal, (fdott ? forcet / fdott : forcet));
                    if (fStatFile.getSaveCurrentTimeStep())
                        fStatFile.getFstream() << getTime() << " " << PI->getIndex() << " " << PJreal->getIndex() << " " << centre << " " << radii_sum - dist << " " << deltat_norm << " " << fdotn << " " << fdott << " " << normal << " " << (fdott ? forcet / fdott : forcet) << std::endl;
                }
            }
        } // end if particle i and j are overlapping
    }

    int Check_and_Duplicate_Periodic_Particle(BaseParticle* i, int nWallPeriodic)
    {
        int C = 0; //Number of particles created
        for (int k = nWallPeriodic; k < getIsPeriodic(); k++)
        { //Loop over all still posible walls
            ///assumes first periodic wall is in x direction
            PeriodicBoundary* perw = static_cast<PeriodicBoundary*>(boundaryHandler.getObject(k));
            if ((k ? true : InPeriodicBox(i)) && (perw->getDistance(*i) < i->getRadius() + getMaxInflowParticleRadius()))
            {
                SphericalParticle F0 = *i;
                perw->shiftPosition(i);

                //If Particle is Mdouble shifted, get correct original particle
                BaseParticle* From = i;
                while (From->getPeriodicFromParticle() != NULL)
                    From = From->getPeriodicFromParticle();
                F0.setPeriodicFromParticle(From);

                particleHandler.copyAndAddObject(F0);
                C++;

                //Check for Mdouble shifted particles
                C += Check_and_Duplicate_Periodic_Particle(particleHandler.getLastObject(), k + 1);
            }
        }
        return (C);
    }

    void writeXBallsScript() const
    {

        std::stringstream file_name;
        std::ofstream script_file;
        file_name << getName() << ".disp";
        script_file.open((file_name.str()).c_str());

        ///First put in all the script lines. All these lines do is move you to the correct directory from any location
        script_file << "#!/bin/bash" << std::endl;
        script_file << "x=$(echo $0 | cut -c2-)" << std::endl;
        script_file << "file=$PWD$x" << std::endl;
        script_file << "dirname=`dirname \"$file\"`" << std::endl;
        script_file << "cd $dirname" << std::endl;
        
        double scale;
        int format;
        
        int width = 1570, height = 860;
        double ratio = getSystemDimensions() < 3 ? (getXMax() - getXMin()) / (getYMax() - getYMin()) : (getXMax() - getXMin()) / (getZMax() - getZMin());
        if (ratio > width / height)
            height = width / ratio;
        else
            width = height * ratio;
        //~ cout << ratio << "r" << endl;
        if (getSystemDimensions() < 3)
        { // dim = 1 or 2
            format = 8;
            if (getXBallsScale() < 0)
            {
                scale = 1.0 / std::max(getYMax() - getYMin(), getXMax() - getXMin());
            }
            else
            {
                scale = getXBallsScale();
            }
        }
        else
        { //dim==3
            format = 14;
            if (getXBallsScale() < 0)
            {
                scale = width / 480 / (getXMax() - getXMin());
            }

            else
            {
                scale = getXBallsScale();
            }

        }

        script_file << "../xballs -format " << format
                << " -f " << dataFile.getFullName()
                << " -s " << scale
                << " -w " << width + 140
                << " -h " << height + 140
                << " -cmode " << getXBallsColourMode()
                << " -cmax -scala 4 -sort -oh 50 "
                << getXBallsAdditionalArguments()
                << " $*";
        if (getXBallsVectorScale() > -1)
        {
            script_file << " -vscale " << getXBallsVectorScale();
        }
        script_file.close();
        
        //This line changes the file permision and gives the owner (i.e. you) read, write and execute permission to the file
#ifdef UNIX
        chmod((file_name.str().c_str()),S_IRWXU);
#endif
    }

    void outputXBallsDataParticlee(const unsigned int i, const unsigned int format, std::ostream& os) const
            {
        os
        << particleHandler.getObject(i)->getPosition().X << " "
                << particleHandler.getObject(i)->getPosition().Y << " "
                << particleHandler.getObject(i)->getPosition().Z << " "
                << particleHandler.getObject(i)->getVelocity().X << " "
                << particleHandler.getObject(i)->getVelocity().Y << " "
                << particleHandler.getObject(i)->getVelocity().Z << " "
                //~ << (particleHandler.getObject(i)->isFixed()?0.0:particleHandler.getObject(i)->Force.X) << " "
                //~ << (particleHandler.getObject(i)->isFixed()?0.0:particleHandler.getObject(i)->Force.Y) << " "
                //~ << (particleHandler.getObject(i)->isFixed()?0.0:particleHandler.getObject(i)->Force.Z) << " "
                << particleHandler.getObject(i)->getRadius() << " "
                << particleHandler.getObject(i)->getIndSpecies() << " "
                << particleHandler.getObject(i)->getOrientation().Y << " "
                << particleHandler.getObject(i)->getOrientation().Z << " "
                << particleHandler.getObject(i)->getAngularVelocity().X << " "
                << particleHandler.getObject(i)->getAngularVelocity().Y << " "
                << particleHandler.getObject(i)->getAngularVelocity().Z << " "
                << particleHandler.getObject(i)->getIndSpecies() << std::endl;
    }

    double get_PeriodicBoxLength()
    {
        return PeriodicBoxLength;
    }
    void set_PeriodicBoxLength(double new_)
    {
        PeriodicBoxLength = new_;
    }
    int get_PeriodicBoxNSpecies()
    {
        return PeriodicBoxNSpecies;
    }
    void set_PeriodicBoxNSpecies(int new_)
    {
        PeriodicBoxNSpecies = new_;
    }

    bool InPeriodicBox(BaseParticle* P)
    {
        return P->getIndSpecies() < PeriodicBoxNSpecies;
    }

    SphericalParticle inflowParticle_;

private:
    ///stores the length of the periodic box
    double PeriodicBoxLength;
    ///stores the number of species in the periodic box
    int PeriodicBoxNSpecies;

    bool FillChute;
};
