// This file is part of MercuryDPM.
// 
// MercuryDPM is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
// 
// MercuryDPM is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
// 
// You should have received a copy of the GNU General Public License
// along with MercuryDPM.  If not, see <http://www.gnu.org/licenses/>.
// 
// Copyright 2013 The Mercury Developers Team
// For the list of developers, see <http://www.MercuryDPM.org/Team>

#include "Mercury3D.h"
#include <Species/SinterFrictionSpecies.h>
#include <Walls/InfiniteWall.h>
#include <Walls/AxisymmetricIntersectionOfWalls.h>
#include <iomanip>
#include <Interactions/NormalForceInteractions/LinearPlasticViscoelasticInteraction.h>

/** This code loads a simulation with two SinterFrictionSpecies and sinters the particles.
*/
class Sintering : public Mercury3D {
public:

    //set default values
    Sintering(int argc, char *argv[])
    {
        std::stringstream addToName;
        for (unsigned i=1; i<argc; i++)
        {
            if (!strcmp(argv[i], "-H"))
            {
                addToName << 'H' << atof(argv[i + 1]);
            }
            else if (!strcmp(argv[i], "-W"))
            {
                addToName << 'W' << atof(argv[i + 1]);
            }
            else if (!strcmp(argv[i], "-E"))
            {
                addToName << 'E' << atof(argv[i + 1]);
            }
            else if (!strcmp(argv[i], "-D"))
            {
                addToName << 'D' << atof(argv[i + 1]);
            }
            else if (!strcmp(argv[i], "-R"))
            {
                addToName << 'R' << atof(argv[i + 1]);
            }
        }
        setName("SinterBed0"+addToName.str());
        readRestartFile();
        setRestarted(false);
        setName("SinterBed1"+addToName.str());
        logger(INFO,"Name: %",getName());

        particleSpecies = dynamic_cast<SinterFrictionSpecies*>(speciesHandler.getObject(0));
        logger.assert(particleSpecies,"particle species has to be of sinter type");
//        wallSpecies = dynamic_cast<SinterFrictionSpecies*>(speciesHandler.getObject(1));
//        logger.assert(wallSpecies,"wall species has to be of sinter type");
//        wallParticleSpecies = speciesHandler.getMixedObject(wallSpecies, particleSpecies);
    }

    /** creates custom console output */
    void printTime() const override
    {
        static int counter = 0;
        if (++counter==1) writeEneHeader(std::cout);
        writeEneTimestep(std::cout);
    }

    /** creates custom ene header */
    void writeEneHeader(std::ostream& os) const override
    {
    }

    /** creates custom ene output */
    void writeEneTimestep(std::ostream &os) const override {
        Mdouble overlap = 0.0;
        Mdouble plasticOverlap = 0.0;
        for (auto i:interactionHandler) {
            overlap += i->getOverlap();
            auto j = dynamic_cast<const SinterInteraction *>(i);
            if (j) {
                plasticOverlap += j->getPlasticOverlap();
            } else {
                auto j = dynamic_cast<const LinearPlasticViscoelasticInteraction *>(i);
                logger.assert(true, "Species is neither Sinter nor LinearPlasticViscoelastic");
                plasticOverlap += j->getMaxOverlap();
            }
        }
        unsigned N = interactionHandler.getNumberOfObjects();
        overlap /= N;
        plasticOverlap /= N;
        Mdouble radius = particleHandler.getObject(0)->getRadius();

        os
        << "t/tMax " << std::left << std::setw(12) << getTime()/getTimeMax()
        << " eneRatio " << std::left << std::setw(12) << getKineticEnergy() / getElasticEnergy()
        << " delta/r " << std::left << std::setw(12) << overlap/radius
        << " delta0/r " << std::left << std::setw(12) << plasticOverlap/radius
        << std::endl;
    }
    
    SinterFrictionSpecies* particleSpecies; //pointer to the particle properties
//    SinterFrictionSpecies* wallSpecies; //pointer to the particle properties
//    SinterFrictionMixedSpecies* wallParticleSpecies; //pointer to the wall-particle collision properties
};

int main(int argc, char *argv[])
{
	Sintering s(argc,argv);
    Mdouble radius = s.particleHandler.getObject(0)->getRadius();
    s.particleSpecies->setSinterForceAndTime(0.1*3.88772 * 1e-15 * radius/*N*/, 86.8e-3/*s*/, radius); //100 times too quick
    //s.wallParticleSpecies->setSinterForceAndTime(3.88772 * 1e-15 * radius/*N*/, 86.8/*s*/, radius);
    s.setTimeMax(6e3*s.getTimeStep());
    s.setSaveCount(200);
    Mdouble mass = s.particleSpecies->getMassFromRadius(radius);
    logger(INFO, "  ts=r^4*nu/ad=%", radius * radius * radius * radius /s.particleSpecies->getInverseSinterViscosity() / s.particleSpecies->getSinterAdhesion());
    logger(INFO, "  tg=sqrt(d/g)=%", sqrt(2.0*radius/s.getGravity().getLength()));
    logger(INFO, "  tc=sqrt(m/k)=%", sqrt(mass/s.particleSpecies->getLoadingStiffness()));
    std::cout << "tMax=" << s.getTimeMax() << std::endl;
    s.solve();
    
    helpers::writeToFile("SinterBed1.gnu",
                         "set xlabel 'time [s]'\n"
                         "set ylabel ''\n"
                         "r = "+std::to_string(radius)+"\n"
                         "p 'SinterBed1.ene' u 2:($6/r) w lp title 'delta/r', '' u 2:($8/r) w lp title 'delta_0/r'\n"
                         );
    std::cout << "Execute 'gnuplot SinterBed1.gnu' to view output" << std::endl;
//./SinterBed0.xballs -noborder 4 -cmode 7 -cmaxset 1 -s 50000 -w 1205 -moh 540 -mo 50
}
