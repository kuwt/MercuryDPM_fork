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

#include<iostream>
#include <Species/LinearViscoelasticSlidingFrictionSpecies.h>
#include "Mercury3D.h"
#include "Particles/BaseParticle.h"
#include "Walls/InfiniteWall.h"

//#define DEBUG_OUTPUT

class my_problem : public Mercury3D
{

public:

    void setupInitialConditions()
    {

        double particle_mass1 = 4.0 / 3.0 * particle_density1 * constants::pi * pow(particle_radius1, 3);

        double particle_mass2 = 4.0 / 3.0 * particle_density2 * constants::pi * pow(particle_radius2, 3);

        //0 : Walls
        //1,2 Particles
        //3 Base Wall
        auto species0 = speciesHandler.copyAndAddObject(LinearViscoelasticSlidingFrictionSpecies());
        auto species1 = speciesHandler.copyAndAddObject(species0);
        auto species2 = speciesHandler.copyAndAddObject(species0);
        auto species3 = speciesHandler.copyAndAddObject(species0);
        auto species10 = speciesHandler.getMixedObject(species1, species0);
        auto species20 = speciesHandler.getMixedObject(species2, species0);
        auto species21 = speciesHandler.getMixedObject(species2, species1);
        auto species31 = speciesHandler.getMixedObject(species3, species1);
        auto species32 = speciesHandler.getMixedObject(species3, species2);

        //Set up partilce 1 properties - Please note walls are species 0.
        species1->setCollisionTimeAndRestitutionCoefficient(tc, r_particle11, particle_mass1);
        species10->setCollisionTimeAndNormalAndTangentialRestitutionCoefficient(tc, r_wall, 1, 2 * particle_mass1); //note: here, tangential cor is set to 1 (i.e. no tangential friction). Can alternatively set to r_wall, for instance.
        species10->setSlidingFrictionCoefficient(mu_wall);
        species1->setDensity(particle_density1);

        //Set up particle 2 properties
        species2->setCollisionTimeAndRestitutionCoefficient(tc, r_particle22, particle_mass2);
        species20->setCollisionTimeAndNormalAndTangentialRestitutionCoefficient(tc, r_wall, 1, 2 * particle_mass2); //note: here, tangential cor is set to 1 (i.e. no tangential friction). Can alternatively set to r_wall, for instance.
        species20->setSlidingFrictionCoefficient(mu_wall);
        species21->setCollisionTimeAndNormalAndTangentialRestitutionCoefficient(tc, r_particle12, 1, 2 * particle_mass1 * particle_mass2 / (particle_mass1 + particle_mass2));
        species2->setDensity(particle_density2);

        //Set up particle 3 properties
        species31->setCollisionTimeAndNormalAndTangentialRestitutionCoefficient(tc, r_base, 1, 2 * particle_mass1);
        species32->setCollisionTimeAndNormalAndTangentialRestitutionCoefficient(tc, r_base, 1, 2 * particle_mass2);
        species31->setSlidingFrictionCoefficient(mu_wall);
        species32->setSlidingFrictionCoefficient(mu_wall);

        //Make it a 3D problem
        setParticleDimensions(3);
        setSystemDimensions(3);

        //Five solid walls - open top
        InfiniteWall w0;
        w0.set(Vec3D(-1.0, 0.0, 0.0),Vec3D(getXMin(),0.0,0.0));
        w0.setPrescribedPosition([this] (double time)
        {
            double t = getTime()-1.0;
            if (t > 0.0)
            {
                return Vec3D(getXMin(),0.0,shaker_amp * std::sin(t * 2.0 * shaker_freq * constants::pi));
            }
            else
            {
                return Vec3D(getXMin(),0.0,0.0);
            }
        });
        wallHandler.copyAndAddObject(w0);

        w0.set(Vec3D(+1.0, 0.0, 0.0),Vec3D(getXMax(),0.0,0.0));
        w0.setPrescribedPosition([this] (double time)
        {
            double t = getTime()-1.0;
            if (t > 0.0)
            {
                return Vec3D(getXMax(),0.0,shaker_amp * std::sin(t * 2.0 * shaker_freq * constants::pi));
            }
            else
            {
                return Vec3D(getXMax(),0.0,0.0);
            }
        });
        wallHandler.copyAndAddObject(w0);

        w0.set(Vec3D(0.0, -1.0, 0.0), Vec3D(0.0,getYMin(),0.0));
        w0.setPrescribedPosition([this] (double time)
        {
            double t = getTime()-1.0;
            if (t > 0.0)
            {
                return Vec3D(0.0,getYMin(),shaker_amp * std::sin(t * 2.0 * shaker_freq * constants::pi));
            }
            else
            {
                return Vec3D(0.0,getYMin(),0.0);
            }
        });
        wallHandler.copyAndAddObject(w0);

        w0.set(Vec3D(0.0, +1.0, 0.0), Vec3D(0.0,getYMax(),0.0));
        w0.setPrescribedPosition([this] (double time)
        {
            double t = getTime()-1.0;
            if (t > 0.0)
            {
                return Vec3D(0.0,getYMax(),shaker_amp * std::sin(t * 2.0 * shaker_freq * constants::pi));
            }
            else
            {
                return Vec3D(0.0,getYMax(),0.0);
            }
        });
        wallHandler.copyAndAddObject(w0);

        w0.set(Vec3D(0.0, 0.0, -1.0), Vec3D(0.0,0.0,getZMin()));
        w0.setIndSpecies(3);
        w0.setPrescribedPosition([this] (double time)
        {
            double t = getTime()-1.0;
            if (t > 0.0)
            {
                return Vec3D(0.0,0.0,getZMin() + shaker_amp * std::sin(t * 2.0 * shaker_freq * constants::pi));
            }
            else
            {
                return Vec3D(0.0,0.0,getZMin());
            }
        });
        wallHandler.copyAndAddObject(w0);
        
        //Put the particles on a grid with small random velocities
        BaseParticle p0;
        double max_radius = std::max(particle_radius1, particle_radius2);
        unsigned int N = numberOfParticles;
        double x = max_radius * 1.01;
        double y = max_radius * 1.01;
        double z = getZMin() + max_radius * 1.01;
        for (int i = 0; i < N; i++)
        {
            x += max_radius * 2.01;
            if (x > getXMax() - max_radius)
            {
                x = max_radius * 1.01;
                y += max_radius * 2.01;
            }
            if (y > getYMax() - max_radius)
            {
                x = max_radius * 1.01;
                y = max_radius * 1.01;
                z += max_radius * 2.01;
            }

            p0.setPosition(Vec3D(x, y, z));
            random.randomise(); //TURN THIS OFF IF I WANT REPRODUCIBLE DATA!!
            p0.setVelocity(Vec3D(random.getRandomNumber(-0.01, 0.01), random.getRandomNumber(-0.01, 0.01), random.getRandomNumber(-0.01, 0.01)));

            //species one for even numbers, species 2 for odd --> should give me an initially well-mixed system!
            if (i % 2 == 0)
            {
                p0.setRadius(particle_radius1);
                p0.setSpecies(speciesHandler.getObject(1));
            }
            else
            {
                p0.setRadius(particle_radius2);
                p0.setIndSpecies(2);
            }
            particleHandler.copyAndAddObject(p0);
        }
        setGravity(Vec3D(0.0, 0.0, -9.81));

        particleHandler.computeAllMasses();

    }

    void set_ParticleRadius(double pr1, double pr2)
    {
        particle_radius1 = pr1;
        particle_radius2 = pr2;
    }

    void set_WallCOR(double cor)
    {
        r_wall = cor;
    }

    void set_WallFriction(double mu)
    {
        mu_wall = mu;
    }

    void set_ParticleCOR(double c11, double c22, double c12)
    {
        r_particle11 = c11;
        r_particle22 = c22;
        r_particle12 = c12;
    }

    void set_ParticleCOR(double cor)
    {
        r_particle11 = cor;
        r_particle12 = cor;
        r_particle22 = cor;
    }

    void set_ParticleDensity(double density1, double density2)
    {
        particle_density1 = density1;
        particle_density2 = density2;
    }

    void set_CollisionTime(double tc_in)
    {
        tc = tc_in;
    }

    void set_NumberOfParticles(int np)
    {
        numberOfParticles = np;
    }

    void set_Frequency(double f)
    {
        shaker_freq = f;
    }

    void set_Amplitude(double a)
    {
        shaker_amp = a;
    }

    void set_BaseCOR(double cor)
    {
        r_base = cor;
    }

    void set_switch_plate_amplitude(double time_, double amp_)
    {
        switch_time = time_;
        switch_amp = amp_;
    }
    
protected:

    void actionsBeforeTimeStep()
    {
        //After t=1.0 start to move the bottom wall
        double t = getTime();
		///todo{DK: This is the old moving wall implementation, the new one has to be tested}
        /*if (t > 1.0)
        {
            wallHandler.getObject(4)->move(getZMin() - shaker_amp * sin(t * 2.0 * shaker_freq * constants::pi));

            double vel = shaker_amp * 2.0 * shaker_freq * constants::pi * cos(t * 2.0 * shaker_freq * constants::pi);

            wallHandler.getObject(0)->move(Vec3D(0, 0, vel), getTimeStep());

            wallHandler.getObject(1)->move(Vec3D(0, 0, vel), getTimeStep());

            wallHandler.getObject(2)->move(Vec3D(0, 0, vel), getTimeStep());

            wallHandler.getObject(3)->move(Vec3D(0, 0, vel), getTimeStep());
        }*/

        if (t > switch_time)
            set_Amplitude(switch_amp);

    }

private:

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

    double shaker_amp;
    double shaker_freq;

    double switch_time;
    double switch_amp;

};

int main(int argc UNUSED, char *argv[] UNUSED)
{

    //Set problem up
    my_problem problem;
    problem.setName("f07030");
    problem.setTimeMax(21);

    //Set Container Geometry
    problem.setXMax(100. / 1000.);
    //vector<int> study_num=problem.get_numbers(1,12);
    //problem.setXMax(40.0+5.0*study_num[1]);

    problem.setYMax(25. / 1000.);
    problem.setZMax(400.0 / 1000.);

    //Set Partilce and Wall properties
    //
    //214 5mm glass particles
    problem.set_NumberOfParticles(250);
    problem.set_ParticleRadius(2.50002 / 1000, 2.5 / 1000);
    problem.set_ParticleDensity(2500, 2500); //2500 & 7850 --> glass & steel.
    double tc = 1e-4;
    problem.set_CollisionTime(tc);
    problem.set_WallCOR(0.7);
    problem.set_WallFriction(0.0); //
    problem.set_ParticleCOR(0.91, 0.91, 0.91); //(11, 22, 12) OR just leave as a single value for equal COR.
    problem.set_BaseCOR(0.6);
    
    problem.set_switch_plate_amplitude(3.0, 0.862 / 1000); //(time, new amplitude)

    //if (study_num[0] > 0)
    //		{
    //			std::cout << "Whole study has started" << std::endl;
    //			exit(0);
    //		}
    //	else
    //If the study is not complete save the data to disk and move on
    //		{
    //std::cout << "Going to launch a second/third/... code" <<std::endl;
    //problem.launch_new("Param");
    //		}
    
    //Shaker prop
//This is fake lowering the frequence and raising amplitude by a factor of 5.

    problem.set_Frequency(70);
    problem.set_Amplitude(0.862 / 1000); //normally 0.862
    //problem.set_Amplitude(2.5/1000);
    
    //Now run the code and solve  - time is set here because I am using the autodetect
    problem.autoNumber();
    problem.setTimeStep(tc / 50);
    problem.setSaveCount(1000 * 21);
    
    problem.solve();
}
