//
// Created by irana on 7/28/17.
//

#include <fstream>
#include <iomanip>

#include "StatisticsVector.h"
#include "CG/TimeAveragedCG.h"
#include "CG/CG.h"
#include "Logger.h"
#include "BidispersedChute.h"
#include "Species/LinearViscoelasticSlidingFrictionSpecies.h"

class SmallPecletPeriodicChute : public BidispersedChute
{
public:
    ///In the constructor, we make a file for the centre of mass (COM) of the bulk, small and large particles over time
    ///Here, the file is opened and the header is written
    ///Additionally, we have no polydispersity at all around each particle diameter.
    SmallPecletPeriodicChute(BidispersedChuteParameters parameters) : BidispersedChute(parameters)
    {
        
        // since the size-ratio is so small here, we do not want a distribution for the particle radii
        setStandardDeviation(0);
    }
    
    ///In the destructor, we have to close the file with centres of mass.
    ~SmallPecletPeriodicChute()
    {
        comFile.close();
    }
    
    ///Here, we set the property of the species. Note, that the bottom is made of small particles, and that the micro-
    ///scopic friction coefficient (slidingFrictionCoefficient) is different for large and small particles. For the
    ///interaction between small and large particles, we use the standard microscopic friction computed from these, here
    ///that is the harmonic mean.
    void setSpeciesProperties() override
    {
        const Mdouble density = 6.0 / constants::pi;
        
        const Mdouble massSmall = 4.0 / 3 * constants::pi * pow(parameters.getSmallParticleRadius(), 3.0) * density;
        const Mdouble massLarge = 4.0 / 3 * constants::pi * pow(parameters.getLargeParticleRadius(), 3.0) * density;
        const Mdouble massFixed = 4.0 / 3 * constants::pi * pow(parameters.getFixedParticleRadius(), 3.0) * density;
        logger(INFO, "mass large: %", massLarge);
        
        auto sReference = LinearViscoelasticSlidingFrictionSpecies();
        //for the reference particles (d=1, m=1, see silbert):
        sReference.setDensity(density);
        sReference.setDissipation(25); //gamma^n
        sReference.setSlidingDissipation(sReference.getDissipation()); //  gamma^t
        sReference.setStiffness(2e5); // k^n
        sReference.setSlidingStiffness(2.0 / 7 * sReference.getStiffness()); // k^t
        sReference.setSlidingFrictionCoefficient(0.5); //mu
        
        const Mdouble r_c = sReference.getRestitutionCoefficient(1);
        const Mdouble tc_1 = sReference.getCollisionTime(1);
        
        //for the large particles
        auto sLarge = sReference;
        sLarge.setCollisionTimeAndRestitutionCoefficient(tc_1, r_c, massLarge);
        sLarge.setSlidingDissipation(sLarge.getDissipation());
        sLarge.setSlidingStiffness(2.0 / 7 * sLarge.getStiffness());
        sLarge.setSlidingFrictionCoefficient(0.5);
        
        //for the small particles
        auto sSmall = sReference;
        sSmall.setDensity(density);
        sSmall.setCollisionTimeAndRestitutionCoefficient(tc_1, r_c, massSmall);
        sSmall.setSlidingDissipation(sSmall.getDissipation()); //  gamma^t
        sSmall.setSlidingStiffness(2.0 / 7 * sSmall.getStiffness()); // k^t
        sSmall.setSlidingFrictionCoefficient(0.25); //mu
    
        //for the small particles
        auto sFixed = sReference;
        sSmall.setDensity(density);
        sSmall.setCollisionTimeAndRestitutionCoefficient(tc_1, r_c, massFixed);
        sSmall.setSlidingDissipation(sSmall.getDissipation()); //  gamma^t
        sSmall.setSlidingStiffness(2.0 / 7 * sSmall.getStiffness()); // k^t
        sSmall.setSlidingFrictionCoefficient(0.25); //mu
        
        //add all species to handler. Note that this must be done after setting all properties, since otherwise the mixed
        //species don't make sense
        speciesHandler.copyAndAddObject(sReference);
        speciesHandler.copyAndAddObject(sLarge);
        speciesHandler.copyAndAddObject(sSmall);
        speciesHandler.copyAndAddObject(sFixed);
        logger(INFO, "Added % species", speciesHandler.getSize());
    }
    
    ///overwriting some of the default options, here we set the species of the bottom particles to no. 3.
    void setupInitialConditions() override
    {
        BidispersedChute::setupInitialConditions();
        for (BaseParticle* particle : particleHandler)
        {
            if (particle->isFixed())
            {
                particle->setSpecies(speciesHandler.getObject(3));
            }
        }
        std::string filename = "COM" + std::to_string(parameters.getSmallParticleRadius()) + "_"
                               + std::to_string(parameters.getAngleInDegrees());
        comFile.open(filename, std::ofstream::out);
        comFile << std::setw(12);
        comFile << "time" << std::setw(12) << "com_flow" << std::setw(12) << "com_large" << std::setw(12)
                << "com_small" << std::endl;
    }
    
    void actionsOnRestart() override
    {
        std::string filename = "restartedCOM" + std::to_string(parameters.getSmallParticleRadius()) + "_"
                               + std::to_string(parameters.getAngleInDegrees());
        comFile.open(filename, std::ofstream::out);
        comFile << std::setw(12);
        comFile << "time" << std::setw(12) << "com_flow" << std::setw(12) << "com_large" << std::setw(12)
                << "com_small" << std::endl;
        setName("restarted" + getName());
    }
    
    ///compute the centre of mass of a) all flow particles b) all large particles c) all small particles
    ///also check if the flow is arrested.
    ///note, that we do not check every time step, but only every 1t.
    void actionsAfterTimeStep() override
    {
        static Mdouble nextCheckedTime = 1;
        if (getTime() > nextCheckedTime)
        {
            //compute COM of species 1
            Mdouble com1 = 0;
            int n1 = 0;
            //compute COM of species 2
            Mdouble com2 = 0;
            int n2 = 0;
            //compute COM of species 1 and 2
            Mdouble com = 0;
            int n = 0;
            for (const BaseParticle* const p : particleHandler)
            {
                if (p->getSpecies()->getIndex() == 1)
                {
                    com1 += p->getPosition().Z;
                    n1++;
                    com += p->getPosition().Z;
                    n++;
                }
                else if (p->getSpecies()->getIndex() == 2)
                {
                    com2 += p->getPosition().Z;
                    n2++;
                    com += p->getPosition().Z;
                    n++;
                }
            }
            com /= n;
            com1 /= n1;
            if (n2 > 0)
                com2 /= n2;
            Mdouble differenceNew = com1 - com2;
            nextCheckedTime += 1;
            logger(INFO, "COM species 1: %, COM species 2: %", com1, com2);
            int width = 10;
            comFile << std::setw(12);
            comFile << getTime() << std::setw(12) << com << std::setw(12) << com1 << std::setw(12) << com2 << std::endl;
            if (getKineticEnergy() < 1e-5)
                logger(ERROR, "The flow has arrested");
        }
    }

private:
    std::ofstream comFile;
};

int main(int argc, char* argv[])
{
    //make a chute with angle 23, inflow-height "height"*d^L and 50%/50% small/large particles
    Mdouble height = 6;
    Mdouble angle = 23;
    auto parameters = BidispersedChuteParameters(height, angle, 0.5);
    
    //set the large particle diameter to 1, and the small particle diameter sizeRatio times smaller
    Mdouble sizeRatio = 1.025;
    parameters.setLargeParticleRadius(sizeRatio / (1 + sizeRatio));
    parameters.setSmallParticleRadius(1 / (1 + sizeRatio));
    parameters.setFixedParticleRadius(0.5);
    SmallPecletPeriodicChute problem(parameters);
    
    //set the chute length and chute-width, note that the chute is periodic in both directions
    problem.setChuteLength(20);
    problem.setChuteWidth(6);
    problem.setRoughBottomType(RoughBottomType::MULTILAYER);
    
    //set the name depending on the parameters used
    std::stringstream name;
    name << 'H' << parameters.getInflowHeight()
         << "SizeRatio" << sizeRatio
         << "Width" << problem.getChuteWidth();
    problem.setName(name.str());
    
    //simulate for a long time, write a file every 50t and write a new restart file every time
    // (so do not overwrite earlier restart file)
    problem.setTimeMax(2000);
    problem.setSaveCount(50. / problem.getTimeStep());
    problem.restartFile.setFileType(FileType::ONE_FILE);
    
    //Turn this on if you want coarse-graining
    /*auto cg0 = problem.cgHandler.copyAndAddObject(TimeAveragedCG<CGCoordinates::Z, CGFunctions::Lucy>());
    cg0->statFile.setSaveCount(0.1 / problem.getTimeStep());
    cg0->statFile.setName(problem.getName() + ".LucyZ.stat");
    cg0->setNZ(200);
    cg0->setTimeMin(1900);*/
    problem.solve(argc, argv);
    return 0;
}
