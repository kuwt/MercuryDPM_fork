//Copyright (c) 2013-2023, The MercuryDPM Developers Team. All rights reserved.
//For the list of developers, see <http://www.MercuryDPM.org/Team>.
#include "HorizontalMixer.h"

int main(int argc, char* argv[])
{
    Mdouble revolutionsPerSecond = 0.25;
    Mdouble rotationSpeed = 0.1*constants::pi*revolutionsPerSecond;
    Mdouble particleRadius = 0.02;
    Mdouble timeMin = 0.0;
    Mdouble fillHeight = 2.0;
    HorizontalMixer mixer(particleRadius, rotationSpeed, timeMin, fillHeight);

    //set name
    mixer.setName("HorizontalMixer");

    //remove old files (note this is a bit dangerous)
//    mixer.removeOldFiles();

    //set species
    auto s = mixer.speciesHandler.copyAndAddObject(LinearViscoelasticSlidingFrictionSpecies());
    s->setDensity(2000);
    Mdouble mass = s->getMassFromRadius(particleRadius);
    s->setCollisionTimeAndRestitutionCoefficient(0.01, 0.5, mass);
    s->setSlidingFrictionCoefficient(0.5);
    s->setSlidingStiffness(2.0 / 7.0 * s->getStiffness());
    s->setSlidingDissipation(2.0 / 7.0 * s->getDissipation());
    
    //set timestep
    mixer.setTimeStep(0.2 * s->getCollisionTime(mass));
    //save every 2 collision times for smooth viewable output
    mixer.setSaveCount((unsigned) (200.0 * s->getCollisionTime(mass) / mixer.getTimeStep()));
    logger(INFO, "Savecount: %", mixer.dataFile.getSaveCount());
    
    mixer.setTimeMax(100);
    mixer.writeScript();
    
    mixer.solve(argc, argv);
    
    return 0;
}
