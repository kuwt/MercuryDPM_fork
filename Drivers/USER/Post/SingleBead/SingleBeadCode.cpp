#include <sstream>
#include <iostream>
#include <iomanip>
#include <cmath>
#include "Mercury3D.h"
#include "Species/LinearViscoelasticFrictionSpecies.h"
#include "Walls/AxisymmetricIntersectionOfWalls.h"
#include "Walls/InfiniteWall.h"
#include "Boundaries/PeriodicBoundary.h"
#include <chrono>

class SingleIntruderParticle : public Mercury3D
{
public:
    SingleIntruderParticle() = default;

    void setupInitialConditions() override
    {
        setGravity({0., 0., -9.81});
        setParticleVolumes();
        setSmallParticleSpecies(0.4, 0.1, 0.0, 0.8);
        setLargeParticleSpecies(0.4, 0.1, 0.0, 0.8);
        setDrumSpecies(0.4, 0.1, 0.0, 0.8);
        createMixedSpecies();
        insertParticles();

        // The smallest timestep is coming from the smallest particle species
        setTimeStep(getSmallParticleSpecies()->getCollisionTime(getSmallParticleVolume()*getSmallParticleDensity())/50.0);

        createDrum();
        if (makeEndWalls)
        {
            createDrumSides();
        }
        else
        {
            createDrumPeriodic();
        }

        checkTime = 1.0;
        logger(INFO,"=================================================================================");
        logger(INFO,"Starting to settle particles until getKineticEnergy() < 1.0e-2> and getTime > 1.0");
        logger(INFO,"=================================================================================\n");
    }

    void actionsAfterTimeStep() override
    {
        if ((getTime() >= checkTime) && (getKineticEnergy() <= 1.0e-2))
        {
            logger(INFO,"=================================================================================");
            logger(INFO,"Start drum rotation as we are later than checkTime and getKineticEnergy() < 1.0e-2");
            logger(INFO,"Current KE: %", getKineticEnergy());
            logger(INFO,"=================================================================================\n");

            for (unsigned int i = 0; i < wallHandler.getNumberOfObjects(); i++)
            {
                // Angular velocity is the set velocity transformed to rad/s, hence user sets rotations/sec
                logger(INFO,"Rotation of wallHandler.getObject(%) is set",i);
                wallHandler.getObject(i)->setAngularVelocity(getRotationalVelocity() * 2.0*constants::pi);
            }

            checkTime = getTimeMax();
        }

        if ((getTime() >= checkTime) && (fmod(getTime(),1.0) < getTimeStep()))
        {
            logger(INFO,"Current KE: %", getKineticEnergy());
        }
    }

    /// Particle-parameter set and get functions
    void setLargeParticleRadius(double lpr_) {largeParticleRadius = lpr_;}
    void setLargeParticleDensity(double lpd_) {largeParticleDensity = lpd_;}
    void setLargeParticleStiffness(double lps_) {largeParticleStiffness = lps_;}
    void setSmallParticleRadius(double spr_) {smallParticleRadius = spr_;}
    void setSmallParticleDensity(double spd_) {smallParticleDensity = spd_;}
    void setSmallParticleStiffness(double sps_) {smallParticleStiffness = sps_;}

    double getLargeParticleRadius() {return largeParticleRadius;}
    double getSmallParticleRadius() {return smallParticleRadius;}
    double getLargeParticleDensity() {return largeParticleDensity;}
    double getSmallParticleDensity() {return smallParticleDensity;}
    double getLargeParticleStiffness() {return largeParticleStiffness;}
    double getSmallParticleStiffness() {return smallParticleStiffness;}

    double getSmallParticleVolume() {return smallParticleVolume;}
    double getLargeParticleVolume() {return largeParticleVolume;}
    void setParticleVolumes()
    {
        setLargeParticleVolume();
        setSmallParticleVolume();
    }
    void setSmallParticleVolume()
    {
        smallParticleVolume = 4.0/3.0 * constants::pi * getSmallParticleRadius() * getSmallParticleRadius() * getSmallParticleRadius();
    }
    void setLargeParticleVolume()
    {
        largeParticleVolume = 4.0/3.0 * constants::pi * getLargeParticleRadius() * getLargeParticleRadius() * getLargeParticleRadius();
    }

    /// Species set and get functions
    void setParticleWallMuS(double pwmus) {particleWallMuS = pwmus;}
    void setParticleWallMuR(double pwmur) {particleWallMuR = pwmur;}
    void setParticleWallMuT(double pwmut) {particleWallMuT = pwmut;}
    void setParticleParticleMuS(double ppmus) {particleParticleMuS = ppmus;}
    void setParticleParticleMuR(double ppmur) {particleParticleMuR = ppmur;}
    void setParticleParticleMuT(double ppmut) {particleParticleMuT = ppmut;}
    double getParticleWallMuS() {return particleWallMuS;}
    double getParticleWallMuR() {return particleWallMuR;}
    double getParticleWallMuT() {return particleWallMuT;}
    double getParticleParticleMuS() {return particleParticleMuS;}
    double getParticleParticleMuR() {return particleParticleMuR;}
    double getParticleParticleMuT() {return particleParticleMuT;}

    void setParticleWallMu(double mus, double mur, double mut)
    {
        setParticleWallMuS(mus);
        setParticleWallMuR(mur);
        setParticleWallMuT(mut);
    }
    void setParticleParticleMu(double mus, double mur, double mut)
    {
        setParticleParticleMuS(mus);
        setParticleParticleMuR(mur);
        setParticleParticleMuT(mut);
    }

    void setLargeParticleSpecies(double mus, double mur, double mut, double cor)
    {
        largeParticleSpecies = speciesHandler.copyAndAddObject(LinearViscoelasticFrictionSpecies());

        logger.assert_always(largeParticleDensity > 0, "The largeParticleDensity is < 0, please set a value");
        largeParticleSpecies->setDensity(getLargeParticleDensity());
        largeParticleSpecies->setStiffness(getLargeParticleStiffness());

        double massLargeParticle = 4. / 3. * constants::pi * getLargeParticleRadius() * getLargeParticleRadius() * getLargeParticleRadius() * getLargeParticleDensity();
        double tc = largeParticleSpecies->getCollisionTime(massLargeParticle);
        largeParticleSpecies->setCollisionTimeAndRestitutionCoefficient(tc, cor, massLargeParticle);

        largeParticleSpecies->setSlidingFrictionCoefficient(mus);
        largeParticleSpecies->setSlidingStiffness(largeParticleSpecies->getStiffness()*2./7.);
        largeParticleSpecies->setSlidingDissipation(largeParticleSpecies->getDissipation()*2./7.);

        largeParticleSpecies->setRollingFrictionCoefficient(mur);
        largeParticleSpecies->setRollingStiffness(largeParticleSpecies->getStiffness()*2.0/7.0);
        largeParticleSpecies->setRollingDissipation(largeParticleSpecies->getDissipation()*2./7.);

        largeParticleSpecies->setTorsionFrictionCoefficient(mut);
        largeParticleSpecies->setTorsionStiffness(largeParticleSpecies->getStiffness()*2.0/7.0);
        largeParticleSpecies->setTorsionDissipation(largeParticleSpecies->getDissipation()*2./7.);
    }
    LinearViscoelasticFrictionSpecies* getLargeParticleSpecies() {return largeParticleSpecies;}

    void setSmallParticleSpecies(double mus, double mur, double mut, double cor)
    {
        smallParticleSpecies = speciesHandler.copyAndAddObject(LinearViscoelasticFrictionSpecies());
        logger.assert_always(smallParticleDensity > 0, "The smallParticleDensity is < 0, please set a value");
        smallParticleSpecies->setDensity(getSmallParticleDensity());
        smallParticleSpecies->setStiffness(getSmallParticleStiffness());

        double massSmallParticle = 4. / 3. * constants::pi * getSmallParticleRadius() * getSmallParticleRadius() * getSmallParticleRadius() * getSmallParticleDensity();
        double tc = smallParticleSpecies->getCollisionTime(massSmallParticle);
        smallParticleSpecies->setCollisionTimeAndRestitutionCoefficient(tc, cor, massSmallParticle);

        smallParticleSpecies->setSlidingFrictionCoefficient(mus);
        smallParticleSpecies->setSlidingStiffness(smallParticleSpecies->getStiffness()*2./7.);
        smallParticleSpecies->setSlidingDissipation(smallParticleSpecies->getDissipation()*2./7.);

        smallParticleSpecies->setRollingFrictionCoefficient(mur);
        smallParticleSpecies->setRollingStiffness(smallParticleSpecies->getStiffness()*2.0/7.0);
        smallParticleSpecies->setRollingDissipation(smallParticleSpecies->getDissipation()*2./7.);

        smallParticleSpecies->setTorsionFrictionCoefficient(mut);
        smallParticleSpecies->setTorsionStiffness(smallParticleSpecies->getStiffness()*2.0/7.0);
        smallParticleSpecies->setTorsionDissipation(smallParticleSpecies->getDissipation()*2./7.);
    }
    LinearViscoelasticFrictionSpecies* getSmallParticleSpecies() {return smallParticleSpecies;}

    void setDrumSpecies(double mus, double mur, double mut, double cor)
    {
        drumWallSpecies = speciesHandler.copyAndAddObject(LinearViscoelasticFrictionSpecies());
        logger.assert_always(smallParticleDensity > 0, "The smallParticleDensity is < 0, please set a value");
        drumWallSpecies->setDensity(getSmallParticleDensity());
        drumWallSpecies->setStiffness(getSmallParticleStiffness());

        double massSmallParticle = 4. / 3. * constants::pi * getSmallParticleRadius() * getSmallParticleRadius() * getSmallParticleRadius() * getSmallParticleDensity();

        double tc = drumWallSpecies->getCollisionTime(massSmallParticle);
        drumWallSpecies->setCollisionTimeAndRestitutionCoefficient(tc, cor, massSmallParticle);

        drumWallSpecies->setSlidingFrictionCoefficient(mus);
        drumWallSpecies->setSlidingStiffness(drumWallSpecies->getStiffness()*2./7.);
        drumWallSpecies->setSlidingDissipation(drumWallSpecies->getDissipation()*2./7.);

        drumWallSpecies->setRollingFrictionCoefficient(mur);
        drumWallSpecies->setRollingStiffness(drumWallSpecies->getStiffness()*2.0/7.0);
        drumWallSpecies->setRollingDissipation(drumWallSpecies->getDissipation()*2./7.);

        drumWallSpecies->setTorsionFrictionCoefficient(mut);
        drumWallSpecies->setTorsionStiffness(drumWallSpecies->getStiffness()*2.0/7.0);
        drumWallSpecies->setTorsionDissipation(drumWallSpecies->getDissipation()*2./7.);
    }
    LinearViscoelasticFrictionSpecies* getDrumWallSpecies() {return drumWallSpecies;}

    void createMixedSpecies()
    {
        logger.assert_always(smallParticleDensity > 0, "The smallParticleDensity is < 0, please set a value");
        double massSmallParticle = 4. / 3. * constants::pi * getSmallParticleRadius() * getSmallParticleRadius() * getSmallParticleRadius() * getSmallParticleDensity();
        double massLargeParticle = 4. / 3. * constants::pi * getLargeParticleRadius() * getLargeParticleRadius() * getLargeParticleRadius() * getLargeParticleDensity();

        auto speciesDrumAndLarge = speciesHandler.getMixedObject(drumWallSpecies,largeParticleSpecies);
        speciesDrumAndLarge->setCollisionTimeAndRestitutionCoefficient(largeParticleSpecies->getCollisionTime(massLargeParticle),
                                                                    largeParticleSpecies->getRestitutionCoefficient(massLargeParticle),
                                                                    massLargeParticle, massLargeParticle);
        speciesDrumAndLarge->setSlidingFrictionCoefficient(getParticleWallMuS());
        speciesDrumAndLarge->setSlidingStiffness(speciesDrumAndLarge->getStiffness()*2.0/7.0);
        speciesDrumAndLarge->setSlidingDissipation(speciesDrumAndLarge->getDissipation()*2./7.);
        speciesDrumAndLarge->setRollingFrictionCoefficient(getParticleWallMuR());
        speciesDrumAndLarge->setRollingStiffness(speciesDrumAndLarge->getStiffness()*2.0/7.0);
        speciesDrumAndLarge->setRollingDissipation(speciesDrumAndLarge->getDissipation()*2./7.);
        speciesDrumAndLarge->setTorsionFrictionCoefficient(getParticleWallMuT());
        speciesDrumAndLarge->setTorsionStiffness(speciesDrumAndLarge->getStiffness()*2.0/7.0);
        speciesDrumAndLarge->setTorsionDissipation(speciesDrumAndLarge->getDissipation()*2./7.);

        auto speciesDrumAndSmall = speciesHandler.getMixedObject(drumWallSpecies,smallParticleSpecies);
        speciesDrumAndSmall->setCollisionTimeAndRestitutionCoefficient(smallParticleSpecies->getCollisionTime(massSmallParticle),
                                                                    smallParticleSpecies->getRestitutionCoefficient(massSmallParticle),
                                                                    massSmallParticle, massSmallParticle);
        speciesDrumAndSmall->setSlidingFrictionCoefficient(getParticleWallMuS());
        speciesDrumAndSmall->setSlidingDissipation(speciesDrumAndSmall->getDissipation()*2./7.);
        speciesDrumAndSmall->setSlidingStiffness(speciesDrumAndSmall->getStiffness()*2.0/7.0);
        speciesDrumAndSmall->setRollingFrictionCoefficient(getParticleWallMuR());
        speciesDrumAndSmall->setRollingStiffness(speciesDrumAndSmall->getStiffness()*2.0/7.0);
        speciesDrumAndSmall->setRollingDissipation(speciesDrumAndSmall->getDissipation()*2./7.);
        speciesDrumAndSmall->setTorsionFrictionCoefficient(getParticleWallMuT());
        speciesDrumAndSmall->setTorsionStiffness(speciesDrumAndSmall->getStiffness()*2.0/7.0);
        speciesDrumAndSmall->setTorsionDissipation(speciesDrumAndSmall->getDissipation()*2./7.);

        auto speciesLargeAndSmall = speciesHandler.getMixedObject(largeParticleSpecies,smallParticleSpecies);
        speciesLargeAndSmall->setCollisionTimeAndRestitutionCoefficient(smallParticleSpecies->getCollisionTime(massSmallParticle),
                                                                  smallParticleSpecies->getRestitutionCoefficient(massSmallParticle),
                                                                  massSmallParticle, massLargeParticle);
        speciesLargeAndSmall->setSlidingFrictionCoefficient(getParticleParticleMuS());
        speciesLargeAndSmall->setSlidingDissipation(speciesLargeAndSmall->getDissipation()*2./7.);
        speciesLargeAndSmall->setSlidingStiffness(speciesLargeAndSmall->getStiffness()*2.0/7.0);
        speciesLargeAndSmall->setRollingFrictionCoefficient(getParticleParticleMuR());
        speciesLargeAndSmall->setRollingStiffness(speciesLargeAndSmall->getStiffness()*2.0/7.0);
        speciesLargeAndSmall->setRollingDissipation(speciesLargeAndSmall->getDissipation()*2./7.);
        speciesLargeAndSmall->setTorsionFrictionCoefficient(getParticleParticleMuT());
        speciesLargeAndSmall->setTorsionStiffness(speciesLargeAndSmall->getStiffness()*2.0/7.0);
        speciesLargeAndSmall->setTorsionDissipation(speciesLargeAndSmall->getDissipation()*2./7.);
    };

    /// Geometry creation, set and get functions
    void createDrum()
    {
        AxisymmetricIntersectionOfWalls drumWall;
        logger.assert_always(drumWallSpecies != nullptr, "Please set the drumWallSpecies");
        drumWall.setSpecies(getDrumWallSpecies());
        drumWall.setAngularVelocity(getInitialRotationalVelocity());
        drumWall.setPosition(Vec3D(0.0, 0.0, 0.0));
        drumWall.setOrientation(Vec3D(0.0, 1.0, 0.0));
        drumWall.addObject({1., 0., 0.},{getDrumRadius(), 0.0, 0.0});
        wallHandler.copyAndAddObject(drumWall);
    }
    void createDrumSides()
    {
        InfiniteWall w;
        w.setSpecies(getDrumWallSpecies());

        w.set(Vec3D(0.,-1.,0.),Vec3D(0.,-0.5*getDrumDepth(),0.));
        w.setAngularVelocity(getInitialRotationalVelocity());
        wallHandler.copyAndAddObject(w);
        w.set(Vec3D(0.,1.,0.),Vec3D(0., 0.5*getDrumDepth(), 0.));
        w.setAngularVelocity(getInitialRotationalVelocity());
        wallHandler.copyAndAddObject(w);
    }
    void createDrumPeriodic()
    {
        PeriodicBoundary b_;
        b_.set({0., 1., 0.},{0., -0.5*getDrumDepth(), 0.},{0., 0.5*getDrumDepth(), 0.});
        boundaryHandler.copyAndAddObject(b_);
    }

    void setDrumDepth(double dd_) {drumDepth = dd_;}
    double getDrumDepth() {return drumDepth;}
    void setDrumRadius(double dr_) {drumRadius = dr_;}
    double getDrumRadius() {return drumRadius;}
    void setDrumVolumeFillFraction(double dvff_) {drumVolumeFillFraction = dvff_;}
    double getDrumVolumeFillFraction() {return drumVolumeFillFraction;}

    void setInitialRotationalVelocity(Vec3D irv_) {initialRotationalVelocity = irv_;}
    Vec3D getInitialRotationalVelocity() {return initialRotationalVelocity;}
    void setRotationalVelocity(Vec3D rv_) {rotationalVelocity = rv_;}
    Vec3D getRotationalVelocity() {return rotationalVelocity;}

    void setMakeEndWalls(bool val_){makeEndWalls = val_;}

    /// Operative functions
    void insertParticles()
    {
        SphericalParticle pSmall;
        SphericalParticle pLarge;
        pSmall.setSpecies(getSmallParticleSpecies());
        pLarge.setSpecies(getLargeParticleSpecies());
        pSmall.setRadius(getSmallParticleRadius());
        pLarge.setRadius(getLargeParticleRadius());

        double drumVolume = (getDrumDepth() * constants::pi * getDrumRadius() * getDrumRadius());
        double volumeDesired = getDrumVolumeFillFraction() * (drumVolume - getLargeParticleVolume());
        int nSmall = std::ceil(volumeDesired / getSmallParticleVolume());
        logger(INFO,"Adding 1 large particle and % small particles",nSmall);

        double r;
        double theta;
        double y;

        // Insert large particle
        r = random.getRandomNumber(0.0, getDrumRadius() - getLargeParticleRadius());
        theta = random.getRandomNumber(0.0,2.0*constants::pi);
        y = random.getRandomNumber(-0.5*getDrumDepth()+getLargeParticleRadius(), 0.5*getDrumDepth() - getLargeParticleRadius());
        pLarge.setPosition({r * std::cos(theta), y, r * std::sin(theta)});
        pLarge.setVelocity({0., 0., 0.});
        particleHandler.copyAndAddObject(pLarge);

        int failCounter = 0;
        while (nSmall > 0)
        {
            r = random.getRandomNumber(0.0, getDrumRadius() - getSmallParticleRadius());
            theta = random.getRandomNumber(0.0,2*constants::pi);
            y = random.getRandomNumber(-0.5*getDrumDepth()+getSmallParticleRadius(), 0.5*getDrumDepth() - getSmallParticleRadius());

            pSmall.setPosition({r * std::cos(theta), y, r * std::sin(theta)});
            pSmall.setVelocity({0., 0., 0.});

            if (MercuryBase::checkParticleForInteraction(pSmall))
            {
                particleHandler.copyAndAddObject(pSmall);
                nSmall--;
                failCounter = 0;
            }
            else
            {
                failCounter++;
            }
            if (failCounter==100)
            {
                logger(INFO,"Failed % amount of times to add a small particle consecutively, adding no more particles",failCounter);
                break;
            }
        }
        hGridRebuild();
        logger(INFO,"Added % particles in the drum",particleHandler.getNumberOfObjects());
    }

private:
    LinearViscoelasticFrictionSpecies* smallParticleSpecies;
    LinearViscoelasticFrictionSpecies* largeParticleSpecies;
    LinearViscoelasticFrictionSpecies* drumWallSpecies;

    double largeParticleRadius;
    double smallParticleRadius;
    double largeParticleDensity;
    double smallParticleDensity;
    double largeParticleVolume;
    double smallParticleVolume;
    double largeParticleStiffness;
    double smallParticleStiffness;

    double particleWallMuS;
    double particleWallMuR;
    double particleWallMuT;
    double particleParticleMuS;
    double particleParticleMuR;
    double particleParticleMuT;

    double drumDepth;
    double drumRadius;
    Vec3D initialRotationalVelocity;
    Vec3D rotationalVelocity;
    double drumVolumeFillFraction;
    bool makeEndWalls;

    double checkTime;
    unsigned int step;

};

int main(int argc, char *argv[])
{
    SingleIntruderParticle problem;

    /// Setting drum parameters
    problem.setMakeEndWalls(true);
    problem.setInitialRotationalVelocity({0.0, 0.0, 0.0}); // Set in rotation/s
    problem.setRotationalVelocity({0.0, 0.1, 0.0}); // Set in rotations/s
    problem.setDrumRadius(0.2);
    problem.setDrumDepth(0.1);
    problem.setDrumVolumeFillFraction(0.25);

    /// Setting particle parameters
    problem.setSmallParticleRadius(5.0e-3);
    problem.setSmallParticleDensity(2500.0);
    problem.setSmallParticleStiffness(1.0e5);

    problem.setLargeParticleRadius(10.0e-3);
    problem.setLargeParticleDensity(2500.0);
    problem.setLargeParticleStiffness(1.0e5);

    problem.setParticleParticleMu(0.4, 0.1, 0.0);
    problem.setParticleWallMu(0.4, 0.1, 0.0);

    /// Overwriting values from command line for easy starting on cluster
    if (argc>1)
    {
        double val_ = 0.0;
        auto readspr = helpers::readFromCommandLine<double>(argc,argv,"-spr",val_);
        if (readspr > 0.0)
        {
            problem.setSmallParticleRadius(readspr);
        }
        auto readlpr = helpers::readFromCommandLine<double>(argc,argv,"-lpr",val_);
        if (readlpr > 0.0)
        {
            problem.setLargeParticleRadius(readlpr);
        }
        auto readspd = helpers::readFromCommandLine<double>(argc,argv,"-spd",val_);
        if (readspd > 0.0)
        {
            problem.setSmallParticleDensity(readspd);
        }
        auto readlpd = helpers::readFromCommandLine<double>(argc,argv,"-lpd",val_);
        if (readlpd > 0.0)
        {
            problem.setLargeParticleDensity(readlpd);
        }
    }

    /// Setting filename
    std::stringstream nameStream;
    std::string nameBase = "SingleIntruder_";
    nameStream << nameBase << "sr_" << (problem.getLargeParticleRadius()/problem.getSmallParticleRadius())
               << "_dr_" << (problem.getLargeParticleDensity()/problem.getSmallParticleDensity());
    std::string fullName = nameStream.str();
    problem.setName(fullName);

	problem.setSaveCount(5000);
    problem.setTimeMax(5.0);
    problem.setMin({-problem.getDrumRadius(), -0.5*problem.getDrumDepth(), -problem.getDrumRadius()});
    problem.setMax({problem.getDrumRadius(), 0.5*problem.getDrumDepth(), problem.getDrumRadius()});

    problem.dataFile.setFileType(FileType::MULTIPLE_FILES);
    problem.restartFile.setFileType(FileType::NO_FILE);
    problem.fStatFile.setFileType(FileType::NO_FILE);
    problem.eneFile.setFileType(FileType::NO_FILE);
    problem.setParticlesWriteVTK(true);
    problem.setWallsWriteVTK(true);

    problem.removeOldFiles();
	problem.solve();
	return 0;
}
