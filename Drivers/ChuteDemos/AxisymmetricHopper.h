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

#include <string.h>

#include "Chute.h"
#include "Walls/AxisymmetricIntersectionOfWalls.h"
#include "Math/ExtendedMath.h"
#include "Species/LinearViscoelasticFrictionSpecies.h"
using namespace std;

///\image html AxisymmetricHopper.png "AxisymmetricHopper"
class AxisymmetricHopper : public Chute
{
public:

    AxisymmetricHopper()
    {
        setXMin(0);
        setYMin(0);
        setZMin(0);
        setXMax(20);
        setYMax(20);
        setZMax(10);
        setTimeMax(0);
        FunnelMinRadius = -1; //radius at funnel outlet
        FunnelMaxRadius = -1; //radius at upper funnel surface
        FunnelHeight = -1; //height/thickness of funnel
        FunnelInflowHeight = -1; //height/thickness of region where particles are included
        num_created = 0;
        max_failed = 100;
        FunnelPointOnAxis.setZero();
    }

    void setupInitialConditions()
    {
        //do not write first time step, as there are zero particles in it
        setLastSavedTimeStep(1);
    
        //Rudi's chute: rmin=H=12.5, rmax=34, dndt=95 500/s= 747 /sqrt(d/g)
        //Rudi's large chute: rmin=14, H=25, rmax=55.5, dndt=215 000/s = 1682 /sqrt(d/g)
        //tmax=3s=380 sqrt(d/g)
        //real chute: dndt=700 000/s = 5477 /sqrt(d/g)

        if (Vec3D::getLength(FunnelPointOnAxis) == 0.0)
            FunnelPointOnAxis = Vec3D(0.5 * (getXMax() + getXMin()), 0.5 * (getYMax() + getYMin()), 0);
        Mdouble PossibleRadius = 0.5 * (getXMax() - getXMin());
        if (FunnelMaxRadius == -1)
            FunnelMaxRadius = PossibleRadius;
        if (FunnelMinRadius == -1)
            FunnelMinRadius = FunnelMaxRadius / 3.;
        if (FunnelHeight == -1)
            FunnelHeight = (FunnelMaxRadius - FunnelMinRadius) / mathsFunc::sin(45. * constants::pi / 180.); //60 deg (2./cos(60.*constants::pi/180.))
        setZMax(getZMax() + FunnelHeight);
        if (FunnelInflowHeight == -1)
            FunnelInflowHeight = min(FunnelHeight / 4., 25. * getInflowParticleRadius());
        
        //create particle properties
        Mdouble tc = .1;

        auto species = speciesHandler.copyAndAddObject(LinearViscoelasticFrictionSpecies());        
        species->setDensity(6.0 / constants::pi);
        //species->setCollisionTimeAndRestitutionCoefficient(tc, 2.0 / 3.0, 0.5 * species->getMassFromRadius(0.5 * (getMinInflowParticleRadius() + getMaxInflowParticleRadius()))));
        species->setCollisionTimeAndRestitutionCoefficient(tc, 2.0 / 3.0, species->getMassFromRadius(0.5 * (getMinInflowParticleRadius() + getMaxInflowParticleRadius())));
        setTimeStep(tc / 10.);
        setSaveCount(helpers::getSaveCountFromNumberOfSavesAndTimeMaxAndTimeStep(50, getTimeMax(), getTimeStep()));
        species->setSlidingStiffness(2. / 7. * species->getStiffness());
        species->setSlidingDissipation(2. / 7. * species->getDissipation());
        species->setSlidingFrictionCoefficient(0.4);
        species->setRollingStiffness(2. / 7. * species->getStiffness());
        species->setRollingDissipation(2. / 7. * species->getDissipation());
        species->setRollingFrictionCoefficient(0.1);

        //set walls
        AxisymmetricIntersectionOfWalls w0;
        w0.setSpecies(speciesHandler.getObject(0));
        //for a prism wall, define points ABC (or more points)
        //as the corners of the prism base (ordered in clockwise
        //direction) and D as the 'periodic direction' of the prism
        //(upwards: Vec3D::Cross(A-B,D) has to point into the wall)
        Vec3D A(FunnelMinRadius, 0, getZMax() - FunnelHeight);
        Vec3D B(FunnelMaxRadius, 0, getZMax());
        Vec3D C(FunnelMaxRadius, 0, getZMax() - FunnelHeight);
        Vec3D D(0, 1, 0); //Periodic direction of the prism
        w0.addObject(Vec3D::cross(A - B, D), A);
        w0.addObject(Vec3D::cross(B - C, D), B);
        w0.addObject(Vec3D::cross(C - A, D), C);
        w0.setPosition(FunnelPointOnAxis);
        w0.setAxis(Vec3D(0, 0, 1));
        wallHandler.copyAndAddObject(w0);

        write(std::cout, false);
    }

    void actionsBeforeTimeStep()
    {
        add_particles();
        cleanChute();
    }

    void create_inflow_particle()
    {
        //the following formula yields polydispersed particle radii:
        P0.setRadius(random.getRandomNumber(getMinInflowParticleRadius(), getMaxInflowParticleRadius()));
        //P0.computeMass(speciesHandler);
        Mdouble R = random.getRandomNumber(0, FunnelMaxRadius - P0.getRadius());
        Mdouble A = random.getRandomNumber(0, 2 * constants::pi);
        double xhat = R * mathsFunc::cos(A);
        double zhat = random.getRandomNumber(getZMax() - FunnelInflowHeight + P0.getRadius(), getZMax() - P0.getRadius());
        double c = mathsFunc::cos(getChuteAngle());
        double s = mathsFunc::sin(getChuteAngle());
        P0.setPosition(FunnelPointOnAxis + Vec3D(xhat * c - zhat * s, R * mathsFunc::sin(A), zhat * c + xhat * s));
        P0.setVelocity(Vec3D(0, 0, 0));
        P0.setSpecies(speciesHandler.getObject(0));
    }

    void add_particles()
    {
        unsigned int failed = 0;

        //try max_failed times to find new insertable particle
        while (failed <= max_failed)
        {
            create_inflow_particle();
            if (checkParticleForInteraction (P0))
            {
                particleHandler.copyAndAddObject(P0);
                failed = 0;
                num_created++;
            }
            else
            {
                failed++;
            }
        };
    }

    virtual void cleanChute()
    {
        //clean outflow every 100 time steps
        static int count = 0, maxcount = 100;
        if (count > maxcount)
        {
            count = 0;
            // delete all outflowing particles
            for (unsigned int i = 0; i < particleHandler.getNumberOfObjects();)
            {
                if (particleHandler.getObject(i)->getPosition().Z < getZMin())
                {
                    particleHandler.removeObject(i);
                }
                else
                    i++;
            }
        }
        else
            count++;
    }

    bool readNextArgument(int& i, int argc, char *argv[])
    {
        if (!strcmp(argv[i], "-inflowHeight"))
        {
            setInflowHeight(atof(argv[i + 1]));
            setZMax(atof(argv[i + 1]));
        }
        else if (!strcmp(argv[i], "-FunnelMinRadius"))
        {
            set_FunnelMinRadius(atof(argv[i + 1]));
        }
        else if (!strcmp(argv[i], "-FunnelMaxRadius"))
        {
            set_FunnelMaxRadius(atof(argv[i + 1]));
        }
        else if (!strcmp(argv[i], "-FunnelHeight"))
        {
            set_FunnelHeight(atof(argv[i + 1]));
        }
        else if (!strcmp(argv[i], "-FunnelInflowHeight"))
        {
            set_FunnelInflowHeight(atof(argv[i + 1]));
        }
        else
            return Chute::readNextArgument(i, argc, argv); //if argv[i] is not found, check the commands in Mercury3D
        return true; //returns true if argv[i] is found
    }

    void set_FunnelMinRadius(Mdouble new_)
    {
        FunnelMinRadius = new_;
    }

    void set_FunnelMaxRadius(Mdouble new_)
    {
        FunnelMaxRadius = new_;
    }

    void set_FunnelHeight(Mdouble new_)
    {
        FunnelHeight = new_;
    }

    void set_FunnelInflowHeight(Mdouble new_)
    {
        if (FunnelHeight < new_)
        {
            cout << "error: FunnelInflowHeight=" << new_ << " < FunnelHeight=" << FunnelHeight << endl;
            exit(-1);
        }
        FunnelInflowHeight = new_;
    }


protected:
    Vec3D FunnelPointOnAxis;
    Mdouble FunnelMinRadius;
    Mdouble FunnelMaxRadius;
    Mdouble FunnelHeight;
    Mdouble FunnelInflowHeight;
    SphericalParticle P0;
    unsigned int max_failed;
    unsigned int num_created;
};
