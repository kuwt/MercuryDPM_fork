//Copyright (c) 2013-2014, The MercuryDPM Developers Team. All rights reserved.
//For the list of developers, see <http://www.MercuryDPM.org/Team>.
#include "TriolietDietFeeder.h"

int main(int argc UNUSED, char *argv[] UNUSED)
{
    Mdouble revolutionsPerSecond = 0.25;
    Mdouble rotationSpeed = 1.0*constants::pi*revolutionsPerSecond;
    Mdouble particleRadius = 0.1;
    Mdouble timeMin = 2.0;

    TriolietDietFeeder tdf(particleRadius, rotationSpeed, timeMin);
    tdf.setName("TriolietDietFeederDemo");
    tdf.setGravity(Vec3D(0,0,-9.8));
    auto s = tdf.speciesHandler.copyAndAddObject(LinearViscoelasticSlidingFrictionSpecies());
    s->setDensity(2000);
    Mdouble mass = s->getMassFromRadius(particleRadius);
    s->setCollisionTimeAndRestitutionCoefficient(0.025, 0.5, mass);
    s->setSlidingFrictionCoefficient(0.5);
    s->setSlidingStiffness(2.0/7.0*s->getStiffness());
    s->setSlidingDissipation(2.0/7.0*s->getDissipation());
    tdf.setTimeStep(0.2 * s->getCollisionTime(mass));
    //save every 2 collision times for smooth viewable output
    tdf.setSaveCount(20.0*s->getCollisionTime(mass)/tdf.getTimeStep());
    std::cout << "Savecount: " << tdf.dataFile.getSaveCount() << std::endl;
    tdf.setTimeMax(4);
    tdf.setBlades({0,3,6,9,11});
//    tdf.fStatFile.setFileType(FileType::MULTIPLE_FILES_PADDED);
//    tdf.dataFile.setFileType(FileType::MULTIPLE_FILES_PADDED);
//    tdf.restartFile.setFileType(FileType::MULTIPLE_FILES_PADDED);
    tdf.setWallsWriteVTK(FileType::ONE_FILE);
    tdf.setParticlesWriteVTK(true);
    tdf.solve();
}
