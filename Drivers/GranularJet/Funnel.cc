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

#include "Funnel.h"
#include "Boundaries/PeriodicBoundary.h"

///This is the actually constructor it is called do both constructors above.
void Funnel::constructor()
{
    set_funH(0);
    set_funnz(0);
    set_funa(0);
    set_funfr(0.33);

}

////////////////////////////////////////////////////////////////////////////////////////////////////
///Writes the restart data for the funnel ///
////////////////////////////////////////////////////////////////////////////////////////////////////

void Funnel::write(std::ostream& os, bool writeAllParticles) const
        {
    Chute::write(os, writeAllParticles);
    os
    << " funa " << get_funa()
            << " funD " << get_funD()
            << " funHf " << get_funHf()
            << " funnz " << get_funnz()
            << " funOx " << get_funOx()
            << " funOy " << get_funOy()
            << " funfr " << get_funfr() << std::endl;

}
//
//////////////////////////////////////////////////////////////////////////////////////////////////////
/////Prints the restart data for the funnel ///
//////////////////////////////////////////////////////////////////////////////////////////////////////
//
//void Funnel::write(std::ostream& os, bool print_all) const
//{
//    Chute::write(os);
//    os
//    << " Funnel angle:  " << get_funa() / constants::pi * 180
//            << " Funnel Diameter: " << get_funD()
//            << " Funnel falling Height: " << get_funHf()
//            << " Funnel number of particles along the funnel:  " << get_funnz()
//            << " Funnel origin:  (" << get_funOx() << " , " << get_funOy() << " )"
//            << " Funnel filling ratio: " << get_funfr();
//
//}

////////////////////////////////////////////////////////////////////////////////////////////////////
///Reads the restart data for the funnel ///
////////////////////////////////////////////////////////////////////////////////////////////////////

void Funnel::read(std::istream& is, ReadOptions opt)
{
    Chute::read(is, opt);
    std::cout << "Funnel read: " << wallHandler.getNumberOfObjects() << std::endl;
    std::string dummy;
    double funOx, funOy;
    is >> dummy >> funa >> dummy >> funD >> dummy >> funHf >> dummy >> funnz >> dummy >> funOx >> dummy >> funOy >> dummy >> funfr;
    set_funO(funOx, funOy);

    //Calculate additional information
    std::cout << "In Funnel::read, calling update_funnel()" << std::endl;
    update_funnel();
}

void Funnel::create_inflow_particle()
{
    //cout <<" void Funnel::create_inflow_particle()"<<endl;
 inflowParticle_.setRadius(random.getRandomNumber(getMinInflowParticleRadius(), getMaxInflowParticleRadius()));
    /// \bug This can no longer be called. Some one should check with Ruud's thesis what the density should be
    //    inflowParticle_.computeMass();

    Vec3D position, velocity;
    position.Z = random.getRandomNumber(Chute::getInflowHeight() - cos(get_funa()) * get_fundiag(), Chute::getInflowHeight());
    position.Y = random.getRandomNumber(get_funOy() - get_funrmax(), get_funOy() + get_funrmax());

    //P0.getPosition().Z=random.get_RN(Chute::getInflowHeight()-cos(get_funa())*get_fundiag(),Chute::getInflowHeight()); // Generate an Z value from the lower bound up to zmax
//	P0.getPosition().Y=random.get_RN(get_funOy()-get_funrmax(),get_funOy()+get_funrmax()); // Generate an Y value in the diameter range of the funnel.
    double funb = acos((inflowParticle_.getPosition().Y - get_funOy()) / get_funrmax()); // Determine the angle beta which the Y value makes with respect to the origin in the XY plane
    position.X = random.getRandomNumber(get_funOx() - sin(funb) * get_funrmax(), get_funOx() + sin(funb) * get_funrmax());
    //	P0.getPosition().X=random.get_RN(get_funOx()-sin(funb)*get_funrmax(),get_funOx()+sin(funb)*get_funrmax()); // Generate an X value on a distance smaller then Y from the origin.
    inflowParticle_.setPosition(position);
    //Apply ChuteAngle rotation:
    //1) Determine the X distance from the current particle to the origin of the funnel
    double Xo = inflowParticle_.getPosition().X - get_funOx();
    //2) Recalculate X: Origin  + Lifting of a particle      + Moving of the origin
    //P0.getPosition().X = get_funOx() + Xo * cos(getChuteAngle()) + (getInflowHeight()-P0.getPosition().Z)*sin(getChuteAngle());
    position.X = random.getRandomNumber(get_funOx() - sin(funb) * get_funrmax(), get_funOx() + sin(funb) * get_funrmax());
    //3) Recalculate Z: Original Z+ Lifting of a particle
    //P0.getPosition().Z = P0.getPosition().Z + Xo * sin(getChuteAngle());
    position.Z = inflowParticle_.getPosition().Z + Xo * sin(getChuteAngle());

    inflowParticle_.setPosition(position);

    velocity.X = getInflowVelocity() * random.getRandomNumber(-getInflowVelocityVariance(), getInflowVelocityVariance()) + getInflowVelocity();
    velocity.Y = getInflowVelocity() * random.getRandomNumber(-getInflowVelocityVariance(), getInflowVelocityVariance());
    velocity.Z = getInflowVelocity() * random.getRandomNumber(-getInflowVelocityVariance(), getInflowVelocityVariance());

    inflowParticle_.setVelocity(velocity);
    //P0.getVelocity().X = InflowVelocity * random.get_RN(-InflowVelocityVariance,InflowVelocityVariance) + InflowVelocity;
    //P0.getVelocity().Y = InflowVelocity * random.get_RN(-InflowVelocityVariance,InflowVelocityVariance);
    //P0.getVelocity().Z = InflowVelocity * random.get_RN(-InflowVelocityVariance,InflowVelocityVariance);

}

////////////////////////////////////////////////////////////////////////////////////////////////////
///This initially set up the particles///
////////////////////////////////////////////////////////////////////////////////////////////////////
void Funnel::setupInitialConditions()
{
    std::cout << " in Funnel::setupInitialConditions " << std::endl;

    //check funnel
    check_funnel();

    //create walls
    create_walls();

    //creates the funnel.
    create_funnel();
    //creates the bottom of the chute
    createBottom();

    write(std::cout);
}

////////////////////////////////////////////////////////////////////////////////////////////////////
///Creates the funnel; based on particles ///
////////////////////////////////////////////////////////////////////////////////////////////////////

void Funnel::create_funnel()
{

    // Define standard fixed particle
    SphericalParticle F0;
    F0.setRadius(getFixedParticleRadius());
    F0.setPosition(Vec3D(0., 0., 0.));
    const double* funh_ = get_funO();
    double zc = Chute::getInflowHeight() - F0.getRadius();
    bool isdone;
    int layercount = 0;
    if (get_funH() <= getFixedParticleRadius())
    {
        isdone = true;
    }
    else
    {
        isdone = false;
    }

    while (!isdone)
    {
        layercount++;
#ifdef DEBUG_OUTPUT
        cout<< "Current layer in Funnel creation = " << layercount << endl;
#endif
        //Current radius of the funnel
        //double funr_=get_funr()-(Chute::getInflowHeight()-zc+F0.Radius)*sin(get_funa())+F0.Radius;
        double funr_ = 1.0 / 2.0 * funD + ((funnz - layercount) * sin(funa)) * 2 * getFixedParticleRadius() + getFixedParticleRadius();
        //Calculate the # of particles in current circle:
        int nc = ceil((constants::pi * funr_) / (F0.getRadius()));
        double dth = 2 * constants::pi / nc; // angle increment.
        double Xo;
        Vec3D position;
        //Loop over those particles
        for (int i = 0; i < nc; i++)
        {
            //Create particles in a circle at certain height zc:
            position.X = funh_[0] + funr_ * cos(i * dth);
            position.Y = funh_[1] + funr_ * sin(i * dth);
            position.Z = zc;
            F0.setPosition(position);

            //Apply ChuteAngle rotation:
            //1) Determine the X distance from the current particle to the origin of the funnel
            Xo = F0.getPosition().X - get_funOx();
            //2) Recalculate X: Origin  + Lifting of a particle      + Moving of the origin
            //F0.getPosition().X = get_funOx() + Xo * cos(getChuteAngle()) + (getInflowHeight()-F0.getPosition().Z)*sin(getChuteAngle());
            position.X = get_funOx() + Xo * cos(getChuteAngle()) + (getInflowHeight() - F0.getPosition().Z) * sin(getChuteAngle());
            //3) Recalculate Z: Original Z+ Lifting of a particle
            //F0.getPosition().Z = F0.getPosition().Z + Xo * sin(getChuteAngle());
            position.Z = F0.getPosition().Z + Xo * sin(getChuteAngle());
            F0.setPosition(position);
            particleHandler.copyAndAddObject(F0);
            //cc++;
        }

        if (zc < (Chute::getInflowHeight() - get_funH()) + getFixedParticleRadius())
        {
            std::cout << "Lst funr_ = " << funr_ << " with zc = " << zc << std::endl;
            isdone = true;
        }

        if (funr_ <= (1.0 / 2.0 * get_funD()))
        {
            std::cout << "Lst funr_ = " << funr_ << std::endl;
            funr_ = 1.0 / 2.0 * get_funD() + 2.0 * F0.getRadius();
            std::cerr << "Warning: Funnel is closing! \n";
            isdone = true;
        }

        zc -= 2 * F0.getRadius() * cos(get_funa());
    }
    //set_Nmax(particleHandler.getNumberOfObjects());

    //finally, fix particles to the floor
    for (std::vector<BaseParticle*>::iterator it = particleHandler.begin(); it != particleHandler.end(); ++it)
        (*it)->fixParticle();

#ifdef DEBUG_OUTPUT
    //Print funnel info to console:
    std::cout << "Current zmax = " << getZMax()
    << "\nCurrent inflowheight = " << Chute::getInflowHeight()
    << "\nFunnel D = " << get_funD()
    << "\nFunnel R = " << get_funr()
    << "\nFunnel A = " << get_funa()/constants::pi*180
    << "\nFunnel HF = " << get_funHf()
    << "\nFunnel H = " << get_funH()
    << "\nFunnel rmax = " << get_funrmax()
    << "\nFunnel diag = " << get_fundiag()
    << "\nFunnel fr = " << get_funfr()
    << "\nFunnel Oy = " << get_funOy()
    << "\nFunnel Ox = " << get_funOx()
    << "\nCurrent number of particles = " << get_NmaxR() << std::endl;
#endif
}

////////////////////////////////////////////////////////////////////////////////////////////////////
///Update the funnel if it did not fit into the simulation boundaries (xmin,xmax, ymin or ymax)s ///
////////////////////////////////////////////////////////////////////////////////////////////////////

void Funnel::update_funnel()
{
    std::cout << " in Funnel::update_funnel()" << std::endl;
    //set_funnz(get_funnz()*(1.0+get_funfr()));
    set_funH(funnz * cos(funa) * 2 * getFixedParticleRadius());
    Chute::setInflowHeight(funHf + funH);
    set_funr(1.0 / 2.0 * funD + (funnz * sin(funa)) * 2 * getFixedParticleRadius() + getFixedParticleRadius());
    //set_funH(Chute::getInflowHeight()-funHf);
    set_fundiag(get_funfr() * get_funnz() * getFixedParticleRadius()); // Length (in m) of the filling length over the cone
    set_funrmax(get_funr() - get_fundiag() * sin(get_funa()) + getFixedParticleRadius());
    //create_walls();
}

////////////////////////////////////////////////////////////////////////////////////////////////////
///Checks the funnel with the boundaries (xmin,xmax, ymin or ymax)s ///
////////////////////////////////////////////////////////////////////////////////////////////////////

void Funnel::check_funnel()
{
    //Call update_funnel to initiate funr
    update_funnel();
    //Check parameters and reinitialise if needed
    if ((get_funr() + get_funOx()) > getXMax())
    {
        double oldxmax = getXMax();
        setXMax(get_funr() + get_funOx());
        std::cerr << "Funnel did not fit in boundaries, changed xmax from " << oldxmax << " to " << getXMax() << std::endl;
    }
    if ((get_funOx() - get_funr()) < getXMin())
    {
        double oldxmin = getXMin();
        setXMin(get_funOx() - get_funr());
        std::cerr << "Funnel did not fit in boundaries, changed xmin from " << oldxmin << " to " << getXMin() << std::endl;
    }

    if ((get_funOy() + get_funr()) > getYMax())
    {
        double oldymax = getYMax();
        setYMax(get_funOy() + get_funr());
        std::cerr << "Funnel did not fit in boundaries, changed ymax from " << oldymax << " to " << getYMax() << std::endl;
    }
    if ((get_funOy() - get_funr()) < getYMin())
    {
        double oldymin = getYMin();
        setYMin(get_funOy() - get_funr());
        std::cerr << "Funnel did not fit in boundaries, changed ymin from " << oldymin << " to " << getYMin() << std::endl;
    }

    update_funnel();
}

////////////////////////////////////////////////////////////////////////////////////////////////////
///Creates or updates the wall according to the boundaries (xmin,xmax, ymin or ymax)s ///
////////////////////////////////////////////////////////////////////////////////////////////////////

void Funnel::create_walls()
{
    //set side walls - solid if not a periodic
    if (Chute::getIsPeriodic())
    {
        PeriodicBoundary b0;
        b0.set(Vec3D(0.0, 1.0, 0.0), getYMin(), getYMax());
        boundaryHandler.copyAndAddObject(b0);
        //set_NWallPeriodic(1);
        wallHandler.clear(); //set_NWall(0);
        //WallsPeriodic[0].set(Vec3D( 0.0, 1.0, 0.0), getYMin(), getYMax());

    }

    //SIDE WALLS: REMOVED, clean chute should be updated to clean sidewalls to (nov 28)
    else
    {

        //set_NWallPeriodic(0);
        wallHandler.clear(); //set_NWall(0);
        //wallHandler.getObject(0)->set(Vec3D( 0.0,-1.0, 0.0), -getYMin());
        //Walls[1].set(Vec3D( 0.0, 1.0, 0.0),  getYMax());

    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////
///Removes particles which left the x-y domain///
////////////////////////////////////////////////////////////////////////////////////////////////////
void Funnel::cleanChute()
{
    //clean outflow every 100 time steps
    static int count = 0, maxcount = 100;
    if (count > maxcount)
    {
        count = 0;
        // delete all outflowing particles
        for (unsigned int i = 0; i < particleHandler.getNumberOfObjects();)
        {
            //if (Particles[i].get_Position().X>getXMax()||Particles[i].get_Position().X<getXMin()||Particles[i].get_Position().Y<getYMin()||Particles[i].get_Position().Y>getYMax()||Particles[i].get_Position().Z+Particles[i].Radius<-getZMax())
            if (particleHandler.getObject(i)->getPosition().Z + 10 * particleHandler.getObject(i)->getRadius() < getZMin())
            {
#ifdef DEBUG_OUTPUT_FULL
                std::cout << "erased:" << particleHandler.getObject(i) << std::endl;
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

bool Funnel::readNextArgument(int& i, int argc, char *argv[])
{
    if (!strcmp(argv[i], "-funnz"))
    {
        set_funnz(atof(argv[i + 1]));
    }
    else if (!strcmp(argv[i], "-funfr"))
    {
        set_funfr(atof(argv[i + 1]));
    }
    else if (!strcmp(argv[i], "-funr"))
    {
        set_funr(atof(argv[i + 1]));
    }
    else if (!strcmp(argv[i], "-funO"))
    {
        set_funO(atof(argv[i + 1]), atof(argv[i + 2]));
        i++; //two arguments
    }
    else if (!strcmp(argv[i], "-fundiag"))
    {
        set_fundiag(atof(argv[i + 1]));
    }
    else if (!strcmp(argv[i], "-funrmax"))
    {
        set_funrmax(atof(argv[i + 1]));
    }
    else if (!strcmp(argv[i], "-funa"))
    {
        set_funa(atof(argv[i + 1]));
    }
    else if (!strcmp(argv[i], "-funH"))
    {
        set_funH(atof(argv[i + 1]));
    }
    else if (!strcmp(argv[i], "-funHf"))
    {
        set_funHf(atof(argv[i + 1]));
    }
    else if (!strcmp(argv[i], "-funD"))
    {
        set_funD(atof(argv[i + 1]));
    }
    else
    {
        //if argv[i] is not found, check the commands in Chute
        return Chute::readNextArgument(i, argc, argv);
    }
    return true; //returns true if argv[i] is found
}

void Funnel::setName_()
{
    std::stringstream name;
    name << "A" << getChuteAngleDegrees()
            << "X" << getChuteLength()
            << "Y" << getChuteWidth()
            << "L" << getFixedParticleRadius() / getInflowParticleRadius()
            << "nz" << get_funnz()
            //~ << "r" << get_funr()
            //~ << "Ox" << get_funOx()
            //~ << "Oy" << get_funOy()
            //~ << "diag" << get_fundiag()
            //~ << "rmax" << get_funrmax()
            //~ << "H" << get_funH()
            //~ << "a" << get_funa()
            //~ << "Hf" << get_funHf()
            << "D" << get_funD();
    setName(name.str().c_str());
}
