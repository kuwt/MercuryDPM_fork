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

#include <sstream>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <sys/types.h>
#include <sys/stat.h>
#include "scr/Chute.h"
#include "scr/Time.h"

using namespace std;
///quick=true uses very soft particles, false uses harder particles
bool quick = true;

class VariableBottom : public Chute
{
public:

    ///sets parameters of particles and time stepping to the L3 type used in Silbert's papers
    void set_silbert_parameters()
    {
        setInflowParticleRadius(.5);
        setFixedParticleRadius(.5);
        setDensity(6 / pi);
        //time stepping
        if (quick)
        {
            setTimeStep(2e-4);
            setCollisionTimeAndRestitutionCoefficient(0.01, .88);
        }
        else
        {
            setTimeStep(1e-4);
            setCollisionTimeAndRestitutionCoefficient(0.005, .88);
        }
        //particle properties
        setSlidingStiffness(2.0 / 7.0 * getStiffness());
        setSlidingDissipation(0.0);
        //setSlidingDissipation(getDissipation());
        setSlidingFrictionCoefficient(0.5);
    }

    ///sets parameters of the chute
    void set_chute_parameters()
    {
        //chute properties
        setChuteAngle(24, 1.0);
        setChuteWidth(10);
        setChuteLength(80);
        setZMax(30);
        setTimeMax(1e20);
        setSaveCount(1e3);
        //create a wall below the chute
        set_NWall(2);
        Walls[0].set(Vec3D(0.0, 0.0, -1.0), 3.4 * MaxInflowParticleRadius);
        Walls[1].set(Vec3D(-1.0, 0.0, 0.0), MaxInflowParticleRadius);
        //create periodic walls
        set_NWallPeriodic(1);
        WallsPeriodic[0].set(Vec3D(0.0, 1.0, 0.0), getYMin(), getYMax());
    }

    void setupInitialConditions()
    {
        createBottom();
    }

    void createBottom()
    {
        unsigned int NFirstBottom = max(1., round(getXMax() / 20 * 0.0));
        unsigned int NSecondBottom = max(1., round(getXMax() / 20 * 1.0) - 1);
        DPMBase FirstBottom, SecondBottom;
        if (!FirstBottom.readDataFile(quick ? "../ini_files/H20A24L0.5M0.5Quick.ini" : "../ini_files/H20A22L0.5M0.5.ini", 14))
        {
            cerr << "1st input data not found exiting " << endl;
            exit(-1);
        }
        if (!SecondBottom.readDataFile("../ini_files/bottom0.5.ini", 14))
        {
            cerr << "2nd input data not found exiting " << endl;
            exit(-1);
        }
        set_Nmax(FirstBottom.get_N() * NFirstBottom + SecondBottom.get_N() * NSecondBottom
                + getChuteLength() * getChuteWidth() * getZMax());
        set_N(0);
        for (int j = 0; j < FirstBottom.get_N(); j++)
        {
            if (FirstBottom.getObjects()[j].Velocity.GetLength2() == 0.0)
                FirstBottom.getObjects()[j].fixParticle();
        }
        for (int j = 0; j < SecondBottom.get_N(); j++)
        {
            if (SecondBottom.getObjects()[j].Velocity.GetLength2() == 0.0)
                SecondBottom.getObjects()[j].fixParticle();
        }
        for (unsigned int i = 0; i < NFirstBottom; i++)
            for (int j = 0; j < FirstBottom.get_N(); j++)
            {
                Particles.push_back(FirstBottom.getObjects()[j]);
                Particles.back().Position.X += i * 20;
            }
        for (unsigned int i = NFirstBottom; i < NFirstBottom + NSecondBottom; i++)
            for (int j = 0; j < SecondBottom.get_N(); j++)
            {
                Particles.push_back(SecondBottom.getObjects()[j]);
                Particles.back().Position.X += i * 20;
            }
    }

    void add_particles()
    {
        static DPMBase Inflow;
        static bool MDcreated = false;
        if (!MDcreated)
        {
            if (Inflow.readDataFile(quick ? "../ini_files/H20A24L0.5M0.5Quick.ini" : "../ini_files/H20A22L0.5M0.5.ini", 14))
            {
                for (int j = 0; j < Inflow.get_N(); j++)
                {
                    if (Inflow.getObjects()[j].Position.X > 2)
                    {
                        Inflow.getObjects()[j] = Inflow.getObjects().back();
                        Inflow.getObjects().pop_back();
                        j--;
                    }
                }
                MDcreated = true;
            }
            else
            {
                cerr << "Input data not found exiting " << endl;
                exit(-1);
            }
        }
        //try to find new insertable particle
        for (int j = 0; j < Inflow.get_N(); j++)
        {
            if (IsInsertable(Inflow.getObjects()[j]))
                num_created++;
        };
    }

    void actionsBeforeTimeStep()
    {
        static int counter = 0;
        if ( ++counter > 10)
        {
            counter = 0;
            add_particles();
        }
        cleanChute();
    }

    void cleanChute()
    {
        //clean outflow every 100 time steps
        static int count = 0, maxcount = 100;
        if (count > maxcount)
        {
            count = 0;
            // delete all outflowing particles
            for (unsigned int i = 0; i < Particles.size();)
            {
                if (Particles[i].Position.X > getXMax() && Particles[i].Position.Z < -2.0 * MaxInflowParticleRadius)
                {
#ifdef DEBUG_OUTPUT_FULL
                    cout << "erased:" << Particles[i] << endl;
#endif
                    particleHandler.removeObject(i);
                }
                else
                    i++;
            }
        }
        else
            count++;
    }

    void printTime() const
    {
        static int Nold = get_N();
        static double told = getTime();
        cout << "t=" << setprecision(3) << left << setw(6) << getTime()
                << ", tmax=" << setprecision(3) << left << setw(6) << getTimeMax()
                << ", N=" << setprecision(3) << left << setw(6) << Particles.size()
                //~ << ", Nmax=" << setprecision(3) << left << setw(6) << Particles.capacity()
                << ", dN/dt=" << setprecision(3) << left << setw(6) << (get_N() - Nold) / (getTime() - told)
                << endl;
        Nold = get_N();
        told = getTime();
        //~ static unsigned int counter=0;
        //~ if (++counter>10) {counter=0; cout.flush();}
    }

};

int main(int argc, char *argv[])
{
    VariableBottom problem;
    if (quick)
        problem.setName("VariableBottomQuick");
    else
        problem.setName("VariableBottom");

    problem.set_silbert_parameters();
    problem.set_chute_parameters();
    problem.readArguments(argc, argv);
    problem.solve();
}
