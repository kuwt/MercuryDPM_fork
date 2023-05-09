#ifndef MERCURYDPM_CALIBRATION_H
#define MERCURYDPM_CALIBRATION_H

#include <vector>
#include <cstring>
#include "Species/LinearViscoelasticFrictionReversibleAdhesiveSpecies.h"
#include "Species/LinearPlasticViscoelasticFrictionSpecies.h"
#include "Math/PSD.h"
#include "Math/ExtendedMath.h"
#include "Mercury3D.h"
#include "DPMBase.h"
#include "cmath"
#include "Logger.h"

using constants::pi;
using mathsFunc::cubic;
using helpers::readFromCommandLine;

class Calibration : public Mercury3D {
protected:
    PSD psd;
    ParticleSpecies* particleSpecies;
    ParticleSpecies* frictionalWallSpecies;
    ParticleSpecies* frictionlessWallSpecies;
    std::string param;

public:
    Calibration (int argc, char *argv[]) : Mercury3D() {
        // string of parameter values added to file name
        param = readFromCommandLine(argc,argv,"-param",std::string(""));
        //define output file settings
        setOutput(helpers::readFromCommandLine(argc,argv,"-output"));
        // set psd
        setPSD(argc, argv);
        // set species
        setSpecies(argc, argv);
    }

    std::string getParam() {return param;}

    // sets default output file settings
    void setOutput(bool output)
    {
        setFileType(FileType::NO_FILE);
        setSaveCount(NEVER);
        restartFile.setFileType(FileType::ONE_FILE);
        eneFile.setFileType(FileType::ONE_FILE);
        if (output) {
            setSaveCount(1000);
            dataFile.setFileType(FileType::ONE_FILE);
            wallHandler.setWriteVTK(FileType::ONE_FILE);
            //setParticlesWriteVTK(true);
        }
        //default xballs arguments
        setXBallsAdditionalArguments("-v0 -solidf");
    }

    // set psd from command line (-psd r1 p1 .. rn pn)
    void setPSD(int argc, char *argv[]) {
        for (unsigned i=0; i<argc-1; ++i) {
            if (!strcmp(argv[i],"-psd")) {
                std::vector<PSD::RadiusAndProbability> psdVector;
                PSD::RadiusAndProbability value;
                // distribution type given
                if (!strcmp(argv[i+1],"logNormal")) {
                    logger.assert_debug(argc>i+5,"Error in logNormal");
                    // https://en.wikipedia.org/wiki/Log-normal_distribution
                    double meanX = std::atof(argv[i+4]);
                    double stdX = std::atof(argv[i+5]);
                    if (!strcmp(argv[i+3],"radius")) {
                    } else if (!strcmp(argv[i+3],"diameter")) {
                        meanX /= 2;
                        stdX /= 2;
                    } else {
                        logger(ERROR,"setPSD: distribution type % not known", argv[i+3]);
                    }
                    double meanX2 = meanX*meanX;
                    double stdX2 = stdX*stdX;
                    double mean = log(meanX2/sqrt(meanX2+stdX2));
                    double std = sqrt(log(1+stdX2/meanX2));
                    double lnRMin = mean-2.5*std;
                    double lnRMax = mean+2.5*std;
                    unsigned n = 21;
                    psdVector.push_back({0.5*exp(lnRMin), 0});
                    if (std>0) {
                        for (int j = 1; j < n; ++j) {
                            double lnR = lnRMin + j / (double) n * (lnRMax - lnRMin);
                            value.radius = exp(lnR);
                            value.probability = 0.5 * (1.0 + erf((lnR - mean) / (sqrt(2) * std)));
                            psdVector.push_back(value);
                            logger(INFO, "psd % %", value.radius, value.probability);
                        }
                    }
                    psdVector.push_back({0.5*exp(lnRMax), 1});
                    if (!strcmp(argv[i+2],"number")) {
                        psd.setPSDFromVector(psdVector, PSD::TYPE::CUMULATIVE_NUMBER_DISTRIBUTION);
                    } else if (!strcmp(argv[i+2],"volume")) {
                        psd.setPSDFromVector(psdVector, PSD::TYPE::CUMULATIVE_VOLUME_DISTRIBUTION);
                    } else {
                        logger(ERROR,"setPSD: distribution type % not known", argv[i+2]);
                    }
                    return;
                }
                // distribution given
                for (unsigned j=i+4; j<argc-2 && argv[j][0]!='-' && argv[j+1][0]!='-'; j+=2) {
                    if (!strcmp(argv[i+3],"radius")) {
                        value.radius = std::atof(argv[j]);
                    } else if (!strcmp(argv[i+3],"diameter")) {
                        value.radius = std::atof(argv[j]) / 2.;
                    } else {
                        logger(ERROR,"setPSD: distribution type % not known", argv[i+3]);
                    }
                    value.probability = std::atof(argv[j+1]);
                    psdVector.push_back(value);
                }
                if (!strcmp(argv[i+1],"cumulative")) {
                    if (!strcmp(argv[i+2],"number")) {
                        psd.setPSDFromVector(psdVector, PSD::TYPE::CUMULATIVE_NUMBER_DISTRIBUTION);
                    } else if (!strcmp(argv[i+2],"volume")) {
                        psd.setPSDFromVector(psdVector, PSD::TYPE::CUMULATIVE_VOLUME_DISTRIBUTION);
                    } else {
                        logger(ERROR,"setPSD: distribution type % not known", argv[i+1]);
                    }
                } else if (!strcmp(argv[i+1],"probability")) {
                    if (!strcmp(argv[i+2],"number")) {
                        psd.setPSDFromVector(psdVector, PSD::TYPE::PROBABILITYDENSITY_NUMBER_DISTRIBUTION);
                    } else if (!strcmp(argv[i+2],"volume")) {
                        psd.setPSDFromVector(psdVector, PSD::TYPE::PROBABILITYDENSITY_VOLUME_DISTRIBUTION);
                    } else {
                        logger(ERROR,"setPSD: distribution type % not known", argv[i+1]);
                    }
                }  else {
                    logger(ERROR,"setPSD: distribution type % not known", argv[i+2]);
                }
                return;
            }
        }
        //if the "-psd" is not found
        logger(ERROR, "-psd argument not found");
    }

    // set species from command line (-species, -density, etc)
    void setSpecies(int argc, char *argv[]) {
        std::string stringSpecies = readFromCommandLine(argc,argv,"-species",std::string(""));
        if (stringSpecies == "LinearViscoelasticFrictionReversibleAdhesiveSpecies") {
            //set species
            auto species = speciesHandler.copyAndAddObject(LinearViscoelasticFrictionReversibleAdhesiveSpecies());
            //set density
            species->setDensity(readFromCommandLine(argc,argv,"-density", constants::NaN));
            double massMin = species->getMassFromRadius(psd.getMinRadius());
            double massD50 = species->getMassFromRadius(psd.getVolumeDx(50));
            //set constantRestitution
            species->setConstantRestitution(readFromCommandLine(argc,argv,"-constantRestitution"));
            //set collisionTime and restitutionCoefficient
            double collisionTime = readFromCommandLine(argc,argv,"-collisionTime",0.0);
            double restitutionCoefficient = readFromCommandLine(argc,argv,"-restitutionCoefficient",1.0);
            species->setCollisionTimeAndRestitutionCoefficient(collisionTime, restitutionCoefficient, massMin);
            // set slidingFriction, rollingFriction, torsionFriction
            species->setSlidingFrictionCoefficient(readFromCommandLine(argc,argv,"-slidingFriction",0.0));
            species->setSlidingStiffness(2. / 7. * species->getStiffness());
            species->setSlidingDissipation(2. / 7. * species->getDissipation());
            species->setRollingFrictionCoefficient(readFromCommandLine(argc,argv,"-rollingFriction",0.0));
            species->setRollingStiffness(2. / 5. * species->getStiffness());
            species->setRollingDissipation(2. / 5. * species->getDissipation());
            species->setTorsionFrictionCoefficient(readFromCommandLine(argc,argv,"-torsionFriction",0.0));
            species->setTorsionStiffness(2. / 5. * species->getStiffness());
            species->setTorsionDissipation(2. / 5. * species->getDissipation());
            // set adhesion
            species->setAdhesionStiffness(species->getStiffness());
            species->setAdhesionForceMax(9.8*massD50*readFromCommandLine(argc,argv,"-bondNumber",0.0));
            // set timeStep
            setTimeStep(0.05 * species->getCollisionTime(massMin));
            // define side-wall species (no friction/cohesion)
            auto frictionlessWallSpecies_ = speciesHandler.copyAndAddObject(species);
            auto mixedSpecies = speciesHandler.getMixedObject(species, frictionlessWallSpecies_);
            mixedSpecies->setRollingFrictionCoefficient(0.0);
            mixedSpecies->setSlidingFrictionCoefficient(0.0);
            mixedSpecies->setAdhesionForceMax(0.0);
            // define drum-wall species (high friction/ no cohesion)
            auto frictionalWallSpecies_ = speciesHandler.copyAndAddObject(species);
            mixedSpecies = speciesHandler.getMixedObject(species, frictionalWallSpecies_);
            mixedSpecies->setRollingFrictionCoefficient(std::max(1.0,species->getSlidingFrictionCoefficient())); //\todo TW infinity does not work (anymore?)
            mixedSpecies->setSlidingFrictionCoefficient(std::max(1.0,species->getRollingFrictionCoefficient()));
            mixedSpecies->setAdhesionForceMax(0.0);
            // cast down
            particleSpecies = species;
            frictionalWallSpecies = frictionalWallSpecies_;
            frictionlessWallSpecies = frictionlessWallSpecies_;
        } else {
            logger(ERROR,"Species % unknown", stringSpecies);
        }
        logger(INFO, "Species: %", *particleSpecies);
        logger(INFO, "Frictionless wall species: %", *speciesHandler.getMixedObject(particleSpecies, frictionlessWallSpecies));
        logger(INFO, "Frictional wall species: %", *speciesHandler.getMixedObject(particleSpecies, frictionalWallSpecies));
        logger(INFO, "Time step %", getTimeStep());
    }

    void writePSDToFile() {
        // convert psd to volume cumulative diameter
        std::ofstream out(getName()+".psd");
        auto vpsd = psd;
//        vpsd.convertCumulativeToProbabilityDensity();
//        vpsd.convertProbabilityDensityNumberDistributionToProbabilityDensityVolumeDistribution();
//        vpsd.convertProbabilityDensityToCumulative();
        out << "Diameter volumeProbability\n";
        for (auto p : vpsd.getParticleSizeDistribution())
            out << p.radius*2.0 << ' ' << p.probability << '\n';
        // convert psd to volume cumulative diameter
        std::ofstream out2(getName()+".psd2");
        out2 << "Diameter volumeProbability\n";
        for (auto p : particleHandler)
            out2 << p->getRadius()*2.0 << '\n';
        //get real distribution
        helpers::writeToFile(getName()+"PSD.m",
                             "psd = importdata('" + getName() + ".psd',' ',1);\n"
                             "d=psd.data(:,1);\n"
                             "p=psd.data(:,2);\n"
                             "plot(d*1e6,p,'k.-',\"LineWidth\",2)\n"
                             "xlabel(\"diameter [um]\")\n"
                             "ylabel(\"number cumulative\")\n"
                             "axis padded\n"
                             "\n"
                             "hold on\n"
                             "psd2 = importdata('" + getName() + ".psd2',' ',1);\n"
                             "d=psd2.data(:,1);\n"
                             "histogram(d*1e6,'Normalization','cdf')\n"
                             "legend('exact','data')\n"
                             "set(legend,'Location','best')\n"
                             "hold off\n"
                             "\n"
                             "saveas(gcf,'" + getName() + ".psd.png')\n");
        logger(INFO,"Run % to plot PSD",getName()+".m");
    }
};

#endif //MERCURYDPM_CALIBRATION_H