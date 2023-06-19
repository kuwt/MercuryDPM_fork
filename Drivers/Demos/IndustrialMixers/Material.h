//Copyright (c) 2013-2023, The MercuryDPM Developers Team. All rights reserved.
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

#ifndef MERCURY_MATERIAL_H
#define MERCURY_MATERIAL_H
#include "Mercury3D.h"
#include "Math/PSD.h"
#include "Species/LinearViscoelasticFrictionReversibleAdhesiveSpecies.h"
using helpers::readFromCommandLine;

class Material : public Mercury3D {
protected:
    PSD psd;
    ParticleSpecies* particleSpecies;
    ParticleSpecies* frictionalWallSpecies;
    ParticleSpecies* frictionlessWallSpecies;

    Material (int argc, char *argv[]) {
        // set psd
        setPSD(argc, argv);
        // set species
        setSpecies(argc, argv);
    }

    // set psd from command line (-psd r1 p1 .. rn pn)
    void setPSD(int argc, char *argv[]) {
        for (unsigned i=0; i<argc-1; ++i) {
            if (!strcmp(argv[i],"-psd")) {
                std::vector<DistributionElements> psdVector;
                DistributionElements value;
                // distribution type given
                if (!strcmp(argv[i + 1], "logNormal"))
                {
                    logger.assert_debug(argc > i + 5, "Error in logNormal");
                    // https://en.wikipedia.org/wiki/Log-normal_distribution
                    double meanX = std::atof(argv[i + 4]);
                    double stdX = std::atof(argv[i + 5]);
                    if (!strcmp(argv[i + 3], "radius"))
                    {
                    }
                    else if (!strcmp(argv[i + 3], "diameter"))
                    {
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
                            value.internalVariable = exp(lnR);
                            value.probability = 0.5 * (1.0 + erf((lnR - mean) / (sqrt(2) * std)));
                            psdVector.push_back(value);
                            logger(INFO, "psd % %", value.internalVariable, value.probability);
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
                for (unsigned j=i+4; j<argc-1 && argv[j][0]!='-' && argv[j+1][0]!='-'; j+=2) {
                    if (!strcmp(argv[i+3],"radius")) {
                        value.internalVariable = std::atof(argv[j]);
                    } else if (!strcmp(argv[i+3],"diameter")) {
                        value.internalVariable = std::atof(argv[j]) / 2.;
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
        Mdouble particleRadius = 1.5e-3;
        psd = PSD::getDistributionNormal(particleRadius,0.025*particleRadius,50);
        logger(WARN, "-psd argument not found; using default psd");
    }

    // set species from command line (-species, -density, etc)
    void setSpecies(int argc, char *argv[]) {
        double density = readFromCommandLine(argc,argv,"-density",1000);
        double collisionTime = readFromCommandLine(argc,argv,"-collisionTime",0.001);
        double restitutionCoefficient = readFromCommandLine(argc,argv,"-restitutionCoefficient",0.5);
        bool constantRestitution = readFromCommandLine(argc,argv,"-constantRestitution");
        double slidingFriction = readFromCommandLine(argc,argv,"-slidingFriction",0.5);
        double rollingFriction = readFromCommandLine(argc,argv,"-rollingFriction",0.0);
        double torsionFriction = readFromCommandLine(argc,argv,"-torsionFriction",0.0);
        double bondNumber = readFromCommandLine(argc,argv,"-bondNumber",0.0);

        //set species
        auto species = speciesHandler.copyAndAddObject(LinearViscoelasticFrictionReversibleAdhesiveSpecies());
        //set density
        species->setDensity(density);
        //set constantRestitution
        species->setConstantRestitution(constantRestitution);
        //set collisionTime and restitutionCoefficient
        double massD50 = species->getMassFromRadius(psd.getVolumeDx(50));
        double massMin = species->getMassFromRadius(psd.getMinRadius());
        species->setCollisionTimeAndRestitutionCoefficient(collisionTime, restitutionCoefficient, massMin);
        // set slidingFriction, rollingFriction, torsionFriction
        species->setSlidingFrictionCoefficient(slidingFriction);
        species->setSlidingStiffness(2. / 7. * species->getStiffness());
        species->setSlidingDissipation(2. / 7. * species->getDissipation());
        species->setRollingFrictionCoefficient(rollingFriction);
        species->setRollingStiffness(2. / 5. * species->getStiffness());
        species->setRollingDissipation(2. / 5. * species->getDissipation());
        species->setTorsionFrictionCoefficient(torsionFriction);
        species->setTorsionStiffness(2. / 5. * species->getStiffness());
        species->setTorsionDissipation(2. / 5. * species->getDissipation());
        // set adhesion
        species->setAdhesionStiffness(species->getStiffness());
        species->setAdhesionForceMax(9.8*massD50*bondNumber);

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

        logger(INFO, "Frictionless wall species: %", *speciesHandler.getMixedObject(particleSpecies, frictionlessWallSpecies));
        logger(INFO, "Frictional wall species: %", *speciesHandler.getMixedObject(particleSpecies, frictionalWallSpecies));
        logger(INFO, "Time step %", getTimeStep());
    }

};


#endif //MERCURY_MATERIAL_H
