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

#include "Mercury3D.h"
#include "Species/SinterFrictionReversibleAdhesiveSpecies.h"
#include "Species/SinterSpecies.h"
#include "Species/HertzianSinterSpecies.h"
#include "Walls/InfiniteWall.h"
#include "Boundaries/PeriodicBoundary.h"
using constants::pi;

/// One particle, sintering slowly to a wall
template <class SpeciesType>
class SingleParticle : public Mercury3D
{
    typedef typename SpeciesType::InteractionType InteractionType;


public:
    /**
     * Default constructor; creates (empty) species, particle, and geometry sets output details
     * @return
     */
    explicit SingleParticle(Mdouble radius=1e-6) {
        //create empty species
        species = speciesHandler.copyAndAddObject(SpeciesType());

        //output settings
        setName("SingleParticleConstantTemperature");
        //dataFile.setFileType(FileType::NO_FILE);
        //restartFile.setFileType(FileType::NO_FILE);
        setXBallsAdditionalArguments(" -v0 -solidf -noborder 4 -p 1");
        setTimeMax(0);
        restartFile.setSaveCount(constants::intMax);

        //set gravity and density
        //setGravity(Vec3D(0, 0, -1e-6*238.732414637843e12));
        setGravity(Vec3D(0, 0, -9.8));
        setGravity(getGravity()*5e5); //any higher, and gravity has an effect
        species->setDensity(1000); //overly heavy

        // Inserting particle
        SphericalParticle P;
        P.setSpecies(species);
        P.setRadius(radius);
        P.setPosition(Vec3D(0, 0, 0));
        particleHandler.copyAndAddObject(P);

        // Inserting particle
        P.setPosition(Vec3D(0, 0, -(2-1e-15)*radius));
        P.fixParticle();
        particleHandler.copyAndAddObject(P);
//        // Inserting wall
//        InfiniteWall W;
//        W.setSpecies(species);
//        W.set(Vec3D(0, 0, -1),Vec3D(0, 0, -(1-1e-15)*radius));
//        wallHandler.copyAndAddObject(W);

        //define domain
        setXMin(-radius);
        setYMin(-radius);
        setZMin(-radius);
        setXMax(radius);
        setYMax(radius);
        setZMax(radius);
    }

    /** set time step */
    void setupInitialConditions() override
    {
        //set time step
        Mdouble mass = species->getSmallestParticleMass();
        setTimeStep(species->computeTimeStep(mass));
        logger(INFO,"Collision time: %, time step: %", 50.*getTimeStep(), getTimeStep());
        if (getTimeMax()==0) {
            setTimeMax(1e3*getTimeStep()*fStatFile.getSaveCount());
            logger(INFO,"Time max: %", getTimeMax());
        }
        logger(INFO, "%", *species);

        helpers::writeToFile(getName()+".gnu",
                             "set logscale x\n"
                              "set logscale y\n"
                              "set xlabel 't'\n"
                              "set ylabel 'x/a'\n"
                              "p '"+getName()+".fstat' u 1:(sqrt($7/"+helpers::to_string(2.0*getXMax())+")) title 'x/d', 0.3*x**(1./"+(species->getSinterType()==SINTERTYPE::PARHAMI_MCKEEPING?'6':'2')+")\n"
        );

    }

    /** creates custom console output */
    void printTime() const override
    {
        static int counter = 0;
        if (!counter) writeEneHeader(std::cout);
        counter++;
        writeEneTimeStep(std::cout);
    }

    /** creates custom ene header */
    void writeEneHeader(std::ostream& os) const override
    {
        os << "time\tdisplacement\trelContactRadius\trelContactRadiusMax\tforce\n";
    }

    /** creates custom ene output */
    void writeEneTimeStep(std::ostream& os) const override
    {
        if (interactionHandler.getNumberOfObjects() == 0)
        {
            if (getTime()>0) logger(INFO, "writeEneTimeStep: No interaction found at time %",getTime());
            return;
        }

        Mdouble mo = dynamic_cast<const InteractionType*>(interactionHandler.getLastObject())->getPlasticOverlap();
        Mdouble o = interactionHandler.getObject(0)->getOverlap();
        Mdouble f = interactionHandler.getObject(0)->getForce().getLength();
        Mdouble r = particleHandler.getMeanRadius();
        Mdouble z = particleHandler.getObject(0)->getPosition().Z;

        os << getTime() //time
           << "\t" << -z //displacement
           << "\t" << sqrt(o/r) //relContactRadius x/d
           << "\t" << sqrt(mo/r) //relPlasticContactRadius
           << "\t" << f/(species->getSinterAdhesion()*particleHandler.getLastObject()->getRadius()) //force
           << std::endl;

        std::cout << getTime() //time
           << "\t" << -z //displacement
         << "\t" << sqrt(o/r) //relContactRadius
         << "\t" << sqrt(mo/r) //relPlasticContactRadius
         << "\t" << f/(species->getSinterAdhesion()*particleHandler.getLastObject()->getRadius()) //force
         << std::endl;
    }

    void scaleMass(Mdouble scale)
    {
        logger(INFO,"Scaling density, inverse gravity, square of timestep, and square of diss. coeff. up by a factor %",scale);
        species->setDensity(scale * species->getDensity());
        setGravity(getGravity()/scale);
        //setTimeStep(sqrt(scale) * getTimeStep());
        species->setDissipation(sqrt(scale) * species->getDissipation());
    }

    SpeciesType* species;
};

int main(int argc UNUSED, char *argv[] UNUSED)
{
    SingleParticle<SinterSpecies> sp;
    sp.scaleMass(1e12);
    //set species properties
    Mdouble kMax=20000;
    sp.species->setPlasticParameters(0.1*kMax, kMax, 0.1*kMax, 0.16);
    Mdouble radius = sp.particleHandler.getSmallestInteractionRadius();
    Mdouble mass = sp.species->getSmallestParticleMass();
    sp.species->setStiffnessAndRestitutionCoefficient(sp.species->getLoadingStiffness(),0.8,mass);
    //flip between volume and surface sintering
    if (true)
    {
        sp.species->setSinterType(SINTERTYPE::CONSTANT_RATE);
        sp.species->setSinterRate(std::pow(0.3,2.0));
    } else {
        sp.species->setSinterType(SINTERTYPE::PARHAMI_MCKEEPING);
        sp.species->setSinterRate(std::pow(0.3,6.0));
    }
    sp.setTimeMax(1);
    sp.setSaveCount(10);
    sp.species->setSinterAdhesion(1200e-9/radius);
    sp.solve();

    logger(INFO,"Execute 'gnuplot %.gnu' to view output",sp.getName());
}
