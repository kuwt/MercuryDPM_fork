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

#include<iostream>
#include <Species/LinearViscoelasticSlidingFrictionSpecies.h>
#include "Mercury3D.h"
#include "Walls/InfiniteWall.h"

//#define DEBUG_OUTPUT

class Binary : public Mercury3D
{

public:

    void setupInitialConditions() override
    {
        // defines species
        //0 : side walls
        //1 : particle species 1
        //2 : particle species 2
        //3 : base wall
        auto species0 = speciesHandler.copyAndAddObject(LinearViscoelasticSlidingFrictionSpecies());
        auto species1 = speciesHandler.copyAndAddObject(species0);
        auto species2 = speciesHandler.copyAndAddObject(species0);
        auto species3 = speciesHandler.copyAndAddObject(species0);
        auto species10 = speciesHandler.getMixedObject(species1, species0);
        auto species20 = speciesHandler.getMixedObject(species2, species0);
        auto species21 = speciesHandler.getMixedObject(species2, species1);
        auto species31 = speciesHandler.getMixedObject(species3, species1);
        auto species32 = speciesHandler.getMixedObject(species3, species2);

        //set densities of particle species 1 and 2
        species1->setDensity(particle_density1);
        species2->setDensity(particle_density2);

        //computes mass of particles of species 1 and 2, based on radius
        double particle_mass1 = 4.0 / 3.0 * particle_density1 * constants::pi * pow(particle_radius1, 3);
        double particle_mass2 = 4.0 / 3.0 * particle_density2 * constants::pi * pow(particle_radius2, 3);

        //Set properties of particle - particle collisions (both species 1)
        species1->setCollisionTimeAndRestitutionCoefficient(tc, r_particle11, particle_mass1);
        //Set properties of particle - particle collisions (both species 2)
        species2->setCollisionTimeAndRestitutionCoefficient(tc, r_particle22, particle_mass2);
        //Set properties of particle - particle collisions (species 1 - species 2)
        species21->setCollisionTimeAndNormalAndTangentialRestitutionCoefficient(tc, r_particle12, 1, 2 * particle_mass1 * particle_mass2 / (particle_mass1 + particle_mass2));

        //Set properties of particle - side wall collisions (species 1)
        //note: here, tangential restitution is set to 1 (i.e. no tangential dissipation). Can alternatively set to r_wall, for instance.
        species10->setCollisionTimeAndNormalAndTangentialRestitutionCoefficient(tc, r_wall, 1, 2 * particle_mass1);
        species10->setSlidingFrictionCoefficient(mu_wall);
        //Set properties of particle - side wall collisions (species 1)
        species20->setCollisionTimeAndNormalAndTangentialRestitutionCoefficient(tc, r_wall, 1, 2 * particle_mass2);
        species20->setSlidingFrictionCoefficient(mu_wall);

        //Set properties of particle - base wall collisions (species 1)
        species31->setCollisionTimeAndNormalAndTangentialRestitutionCoefficient(tc, r_base, 1, 2 * particle_mass1);
        species31->setSlidingFrictionCoefficient(mu_wall);
        //Set properties of particle - base wall collisions (species 2)
        species32->setCollisionTimeAndNormalAndTangentialRestitutionCoefficient(tc, r_base, 1, 2 * particle_mass2);
        species32->setSlidingFrictionCoefficient(mu_wall);

        //general wall properties
        InfiniteWall w0;
        w0.setSpecies(species3);

        // set left side wall in x-direction
        w0.set(Vec3D(-1.0, 0.0, 0.0),getMin());
        w0.setPrescribedPosition([this] (double time) {
            return Vec3D(getXMin(),getYMin(),getZMin()+shaker_amp *
            std::sin(getTime() * 2.0 * shaker_freq * constants::pi));
        });
        wallHandler.copyAndAddObject(w0);

        // set right side wall in x-direction
        w0.set(Vec3D(1.0, 0.0, 0.0),getMax());
        w0.setPrescribedPosition([this] (double time) {
            return Vec3D(getXMax(),getYMax(),getZMax()+shaker_amp *
            std::sin(getTime() * 2.0 * shaker_freq * constants::pi));
        });
        wallHandler.copyAndAddObject(w0);

        // set left side wall in y-direction
        w0.set(Vec3D(0.0, -1.0, 0.0), getMin());
        w0.setPrescribedPosition([this] (double time) {
            return Vec3D(getXMin(),getYMin(),getZMin()+shaker_amp *
            std::sin(getTime() * 2.0 * shaker_freq * constants::pi));
        });
        wallHandler.copyAndAddObject(w0);

        // set right side wall in y-direction
        w0.set(Vec3D(0.0, 1.0, 0.0),getMax());
        w0.setPrescribedPosition([this] (double time) {
            return Vec3D(getXMax(),getYMax(),getZMax()+shaker_amp *
            std::sin(getTime() * 2.0 * shaker_freq * constants::pi));
        });
        wallHandler.copyAndAddObject(w0);

        // set base wall
        w0.set(Vec3D(0.0, 0.0, -1.0),getMin());
        w0.setPrescribedPosition([this] (double time) {
            return Vec3D(getXMin(),getYMin(),getZMin()+shaker_amp * std::sin(getTime() * 2.0 * shaker_freq * constants::pi));
        });
        wallHandler.copyAndAddObject(w0);

        //Put the particles on a grid with small random velocities
        SphericalParticle p0;
        double max_radius = std::max(particle_radius1, particle_radius2);
        unsigned int N = numberOfParticles;
        double x = max_radius;
        double y = max_radius;
        double z = getZMin() + max_radius;
        for (int i = 0; i < N; i++)
        {
            x += max_radius * 2.;
            if (x > getXMax() - max_radius)
            {
                x = max_radius;
                y += max_radius * 2.;
            }
            if (y > getYMax() - max_radius)
            {
                x = max_radius;
                y = max_radius;
                z += max_radius * 2.;
            }

            p0.setPosition(Vec3D(x, y, z));
            p0.setVelocity(Vec3D(random.getRandomNumber(-0.01, 0.01), random.getRandomNumber(-0.01, 0.01), random.getRandomNumber(-0.01, 0.01)));

            //species one for even numbers, species 2 for odd --> should give an initially well-mixed system
            if (i % 2 == 0) {
                p0.setRadius(particle_radius1);
                p0.setSpecies(species1);
            } else {
                p0.setRadius(particle_radius2);
                p0.setSpecies(species2);
            }
            particleHandler.copyAndAddObject(p0);
        }
        logger(INFO,"Inserted % particles",particleHandler.getNumberOfObjects());
        // gravity
        setGravity(Vec3D(0.0, 0.0, -9.81));
        // write out some information about the simulation
        logger(INFO,"Time scale oscillation: %s", 1./shaker_freq);
        logger(INFO,"Time scale gravity: %s", 2.0*std::min(particle_radius1,particle_radius2)/getGravity().getLength());
        logger(INFO,"Time scales collision: %s", species1->getCollisionTime(particle_mass1));
    }


    void setParticleRadius(double pr1, double pr2)
    {
        particle_radius1 = pr1;
        particle_radius2 = pr2;
    }

    void setWallRestitution(double Restitution)
    {
        r_wall = Restitution;
    }

    void setWallFriction(double mu)
    {
        mu_wall = mu;
    }

    void setParticleRestitution(double c11, double c22, double c12)
    {
        r_particle11 = c11;
        r_particle22 = c22;
        r_particle12 = c12;
    }

    void setParticleRestitution(double Restitution)
    {
        r_particle11 = Restitution;
        r_particle12 = Restitution;
        r_particle22 = Restitution;
    }

    void setParticleDensity(double density1, double density2)
    {
        particle_density1 = density1;
        particle_density2 = density2;
    }

    void setCollisionTime(double tc_in)
    {
        tc = tc_in;
    }

    void setNumberOfParticles(int np)
    {
        numberOfParticles = np;
    }

    void setFrequency(double f)
    {
        shaker_freq = f;
    }

    void setInitialAmplitude(double a)
    {
        shaker_amp = a;
    }

    void setBaseRestitution(double Restitution)
    {
        r_base = Restitution;
    }

    void setSwitchPlateAmplitude(double time_, double amp_)
    {
        switch_time = time_;
        switch_amp = amp_;
    }
    
protected:

    void actionsBeforeTimeStep() override
    {
        if (getTime() > switch_time) {
            shaker_amp = switch_amp;
        }
    }

private:
    //parameters set in the main
    double particle_radius1;
    double particle_radius2;
    double particle_density1;
    double particle_density2;
    double tc;

    double r_wall;
    double r_base;
    double mu_wall;
    double r_particle12;
    double r_particle22;
    double r_particle11;

    unsigned int numberOfParticles;

    //when to start shaker
    double shaker_amp;
    double shaker_freq;

    //second switch
    double switch_time=constants::inf;
    double switch_amp;

};

int main() {

    // Set up a problem of type binary
    Binary dpm;
    // Set name of output files
    dpm.setName("Binary");
    // Set simulation time
    dpm.setTimeMax(21);
    // Set container size
    dpm.setMax({0.1,0.025,0.4});
    // Set particle and wall properties (5mm glass particles)
    // number of particles
    dpm.setNumberOfParticles(250);
    // radius of particles species 1/2
    dpm.setParticleRadius(2.50002e-3, 2.5e-3);
    // density of particles species 1/2
    dpm.setParticleDensity(2500, 2500); //2500 & 7850 --> glass & steel.
    // set stiffness and dissipation such that you get a constant collision time and coefficient of restitution
    double tc = 1e-4;
    dpm.setCollisionTime(tc);
    // coefficient of restitution or particle - particle collisions (for species 1-1, 1-2, and 2-2)
    dpm.setParticleRestitution(0.91, 0.91, 0.91); //(11, 22, 12) OR just leave as a single value for equal Restitution.
    // coefficient of restitution or particle - base collisions (for all species)
    dpm.setBaseRestitution(0.6);
    // coefficient of restitution or particle - side wall collisions (for all species)
    dpm.setWallRestitution(0.7);
    // friction coefficient of walls
    dpm.setWallFriction(0.0); //
    // frequency of base
    dpm.setFrequency(70);
    // amplitude of base
    dpm.setInitialAmplitude(0.862e-3); //normally 0.862
    // at a certain time, switch amplitude
    //dpm.setSwitchPlateAmplitude(3.0, 0.862e-3); //(time, new amplitude)
    // uncomment to number your output files (prevents old simulations from being overwritten by new simulations)
    //dpm.autoNumber();
    // set time step
    dpm.setTimeStep(tc / 50);
    // set frequency of writing output
    dpm.setSaveCount(21000);
    // randomise seed //TURN THIS OFF IF I WANT REPRODUCIBLE DATA!!
    //random.randomise();
    // start the solver
    dpm.solve();
}
