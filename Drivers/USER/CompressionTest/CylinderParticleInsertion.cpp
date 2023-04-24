#include <Mercury3D.h>
#include <Species/LinearPlasticViscoelasticFrictionSpecies.h>
#include <Boundaries/CubeInsertionBoundary.h>
#include <Walls/AxisymmetricIntersectionOfWalls.h>
#include <CG/Coordinates/O.h>
#include <CG/CG.h>
#include <CMakeDefinitions.h>
#include <Boundaries/CubeDeletionBoundary.h>


class CylinderParticleInsertion : public Mercury3D
{

public:
    CG<CGCoordinates::O>* cg;
    PSD psd;
    
    explicit CylinderParticleInsertion(Mdouble particleDensity, Mdouble cylScaling)
    {
        //-----------------
        //Global parameters
        // Density
        std::string pD = helpers::to_string(particleDensity);
        // Scale the cylinder to fit the maximum radius of the PSD cylScaling times into one cylinder diameter
        std::string cS = helpers::to_string(cylScaling);
    
        setName("CylinderParticleInsertion_" + pD + "_" + cS);
    
        //read .csv file line by line to get PSD
    
        // set the PSD from cumulative volumetric distribution functions (CSDLactose, CSDMCC, CSDPVP are exemplary
        // CVDFs from QICPIC image analysis measurements)
        //Make particles and walls visible
        /// \todo TP: add InputData automatically to the build folder (define in master cmake file)
        psd.setPSDFromCSV(getMercuryDPMSourceDir() + "/Drivers/USER/CompressionTest/InputData/CVDFLactose.csv",
                          PSD::TYPE::CUMULATIVE_VOLUME_DISTRIBUTION, false, 1000000.0);
        // cut off size distribution if size ratio is too high
        psd.cutHighSizeRatio();
    
        //Calculate important properties; scale cylinder by maximum PSD radius
        cylHeight_ = cylScaling * psd.getMaxRadius();
        cylRad_ = cylScaling * psd.getMaxRadius();
//        Mdouble particleBulkVolume = particleNumber_*cubic(2.0*particleRadius_);
        cylVol_ = constants::pi * cylRad_ * cylRad_ * cylHeight_;
    
        setParticlesWriteVTK(true);
        wallHandler.setWriteVTK(FileType::MULTIPLE_FILES);
        //Set files
        dataFile.setFileType(FileType::ONE_FILE);
        restartFile.setFileType(FileType::ONE_FILE);
        fStatFile.setFileType(FileType::NO_FILE);
        eneFile.setFileType(FileType::ONE_FILE);
        dataFile.setSaveCount(1000);
        //setting system size
        setMax(Vec3D(cylRad_, cylRad_, 2 * cylHeight_));
        setMin(Vec3D(-cylRad_, -cylRad_, 0.0));
        //gravity
        setGravity(Vec3D(0.0, 0.0, -9.81));
        //add Luding-model
        LinearPlasticViscoelasticFrictionSpecies* s1;
        s1 = speciesHandler.copyAndAddObject(LinearPlasticViscoelasticFrictionSpecies());
        // stiffness > pD * g * CylHeight_ * Particlediameter / 0.01 (0.01 = overlap)
        // stiffness = 2e5 * maxMass * g / maxDiameter (Weinhart 2012) (Luding 98)
        const Mdouble k1 =
                2e2 / 3 * particleDensity * constants::pi * mathsFunc::square(psd.getMaxRadius()) * getGravity()
                        .getLength();
//        const Mdouble k1 = 1 * minParticleRadius; //[Stiffness depends on particle radius]
        s1->setPlasticParameters(k1, k1, 0.0, 0.1);
        s1->setDissipation(0.00000001);
        s1->setDensity(particleDensity);
        s1->setSlidingStiffness(s1->getLoadingStiffness() * 2.0 / 7.0);
        s1->setSlidingDissipation(s1->getDissipation() * 2.0 / 5.0);
        s1->setSlidingFrictionCoefficient(0.5);
        s1->setRollingStiffness(s1->getLoadingStiffness() * 2.0 / 7.0);
        s1->setRollingDissipation(s1->getDissipation() * 2.0 / 5.0);
        s1->setRollingFrictionCoefficient(0.3);
        s1->setTorsionStiffness(s1->getLoadingStiffness() * 2.0 / 7.0);
        s1->setTorsionDissipation(s1->getDissipation() * 2.0 / 5.0);
        s1->setTorsionFrictionCoefficient(0.0);
//        LinearPlasticViscoelasticFrictionSpecies* s2;
//        s2 = dynamic_cast<LinearPlasticViscoelasticSpecies*>(speciesHandler.getObject(1));
//        s2->setPlasticParameters(258.5,8*258.5
//                ,0,0);
//        s2->setDissipation(0.00);
//        s2->setDensity(300);
//        LinearPlasticViscoelasticMixedSpecies* mix;
//        mix= dynamic_cast<LinearPlasticViscoelasticMixedSpecies*>(speciesHandler.getMixedObject(s1,s2));
//        mix->setPlasticParameters(258.5,8*258.5
//                ,0,0.001);
//        mix->setDissipation(0.0);
    
        //Set timestep
        Mdouble minParticleRadius = psd.getMinRadius();
        Mdouble maxParticleRadius = psd.getMaxRadius();
        Mdouble minMass = s1->getMassFromRadius(minParticleRadius);
        const Mdouble collisionTime = s1->getCollisionTime(minMass);
        setTimeStep(collisionTime / 10.0);
        
        //add walls
        AxisymmetricIntersectionOfWalls w;
        w.setSpecies(s1);
        w.setPosition(Vec3D(0, 0, 0));
        w.setOrientation(Vec3D(0, 0, 1));
        w.addObject(Vec3D(1, 0, 0), Vec3D(cylRad_, 0, 0));  //Cylindric wall
        wallHandler.copyAndAddObject(w);
        InfiniteWall w1;
        w1.setSpecies(s1);
        w1.set(Vec3D(0, 0, -1), Vec3D(0, 0, 0));    //Bottom wall
        wallHandler.copyAndAddObject(w1);
//        w1.set(Vec3D(0,0,1), Vec3D(0,0,1)); //Top wall
//        wallHandler.copyAndAddObject(w1);
//        /// Infinite wall ///
//        //z-axis bottom
//        InfiniteWall* w;
//        w = dynamic_cast<InfiniteWall*>(wallHandler.getObject(0));
//        w->set(Vec3D(0, 0, -1), Vec3D(0, 0, getZMin()));
//        //z-axis top
//        wTop->set(Vec3D(0, 0, 1), Vec3D(0, 0, getZMax()));
//        //const Mdouble Lz0 = (getZMax()-getZMin());
//        // final (maximum) top-wall position to satisfy the given strain [loading stage].
//        //const Mdouble finalPos = Lz0*(1-strain_);

        // define cubic insertion
        CubeInsertionBoundary insertionBoundary;
    
        // add particles with PSD
        SphericalParticle particle;
//        double insertionVolume = 0.774/s1.getDensity();
//        logger(INFO,"Insertion volume: %", insertionVolume);
        particle.setSpecies(speciesHandler.getObject(0));
//        insertionBoundary.set(particle, 1000, Vec3D(0.018,0.003,0), Vec3D(0.053,0.05,0.085),Vec3D(0,0,0),Vec3D(0,0,0),radius, radius);
//        insertionBoundary.setInitialVolume(insertionVolume/2);
    
        // insert particles in the whole domain (between getMin and getMax) with initial velocity 0 and and a radius between 1 and 2 (uniform number distribution)
        Vec3D posMin = {getXMin(), getYMin(), getZMax() / 2};
        Vec3D posMax = getMax();
        Vec3D velMin = {0, 0, 0.0};
        Vec3D velMax = {0, 0, 0.0};
        unsigned maxFail = 1000;
        double radMin = 1;
        double radMax = 2;
        insertionBoundary.set(particle, maxFail, posMin, posMax, velMin, velMax);
        //set PSD
        insertionBoundary.setPSD(psd);
        insertionBoundary.setManualInsertion(true);
        insertionBoundary.setInitialVolume(1.1 * cylVol_);
//        insertionBoundary.setInitialVolume(0.1*cylVol_);
        // fill the cylinder volume with a given flowrate
//        insertionBoundary.setVolumeFlowRate(1000 * cylRad_);
    
        //add the insertion boundary to the handler
        boundaryHandler.copyAndAddObject(insertionBoundary);

//        Mdouble Vp = 0;
//        double N = particleHandler.getNumberOfObjects();
//        for (int i = 0; i < N; i++)
//        {
//            //particleHandler.getObject(i)->setVelocity(Vec3D(0.0, 0.0, 0.0));
//            particleHandler.getObject(i)->setSpecies(s1);
//            Vp = Vp + particleHandler.getObject(i)->getVolume();
//        }
//    }
//        const Mdouble l0 = (getZMax()-getZMin());
//        const Mdouble finalPos = l0*(1-strain_);
//        // time required for the wall to cross finalPos distance
//        //double loadingTime = strain_ /strainRate_;
//        // time required for the wall to cross finalPos distance
//        loadingTime = -log(finalPos/l0)/strainRate_;
//        // reset Time Max
//        // check if the time max is sufficient for both loading and unloading stages
//        logger(INFO, "Maximum time for Loading and unloading stages is twice the loading time(approximately)!");
//        setTimeMax(2 * loadingTime);
//        std::cout << "Time max: "<< loadingTime <<std::endl;
//        std::cout << "Sample Volume: "<< Vp <<std::endl;
//        // this for setting numbering to start form zero
//        getVtkWriter()->setFileCounter(0);

//    int countDown = 1; // counter for loading
//    int countUp = 1;   // counter for un-loading
//    void actionsAfterTimeStep() override {
//        // initial length
//        Mdouble l0 = (getZMax() - getZMin());
//        Mdouble finalPos = l0*(1-strain_);
//        // ALWAYS PUT THIS COMMANDS FIRST
//        // write a file for stress strain on the top wall
//        /// check! saveCount should be the same here
//        if (getNumberOfTimeSteps()%2000==0) {
//            static std::ofstream myFile("strainStressTopWall.dat");
//            myFile << getTime()
//                   << '\t' << wallHandler.getObject(1)->getPosition().Z
//                   << '\t' << (l0-wallHandler.getObject(1)->getPosition().Z)/l0
//                   << '\t' << wallHandler.getObject(1)->getForce()/(getXMax()*getYMax())
//                   << '\n';
//        }
//        //double loadingTime = -log(finalPos/l0)/strainRate_;
//        //double loadingTime = strain_ /strainRate_;
//        // if the position of the wall now is bigger than finalPos do nothing
//        if (getTime() >= loadingTime) {// UNLOADING STAGE
//            // here l0 = finalPos
//            countUp++;
//            // Calculate the next position of the wall
//            Mdouble le = finalPos * exp(+strainRate_ * countUp * getTimeStep());
//            /*linear calculation*/
//            //Mdouble l = l0*(1-strainRate_*countUp*getTimeStep());
//            // set the new location of the wall
//            wallHandler.getObject(1)->setPosition(Vec3D(0, 0, le));
//            //termination criterion >> when the force in z-direction = 1e-6
//            if (std::abs(wallHandler.getObject(1)->getForce().Z)<=0){
//                setTimeMax(getTime());
//            }
//            return;
//        }
//            /// LOADING STAGE ///
//        else if (finalPos < wallHandler.getObject(1)->getPosition().Z){
//            // why we do not use countDown here ? because time here already start from 1
//            // Calculate the next position of the wall
//            /*logarithmic calculation*/
//            Mdouble l = l0*exp(-strainRate_*getTime());
//            /*linear calculation*/
//            //Mdouble l = l0*(1-strainRate_*countDown*getTimeStep());
//            // set the new location of the wall
//            wallHandler.getObject(1)->setPosition(Vec3D(0,0,l));
//            return;
//        }
//
//    }

//public:
//    // constructor
//    explicit CylinderParticleInsertion(std::string& fileName){
//        setName(fileName);
//        readRestartFile();
//        setRestarted(false);
//        wTop = dynamic_cast<InfiniteWall*>(wallHandler.getObject(1));
//    }
//    // setting the strain
//    // Control the Strain_rate[1/s] which is the time required to have the final state
//    // Control the Strain between (min_value = 0, max_val= 1)
//    void setStrainAndStrainRate(double strain,double strainRate){
//        // necessity condition (0 < strain_ < 1)
//        if (strain_< 0 || strain_ >1){
//            std::cout <<"Strain not in the Limits!" <<std::endl;
//            exit(EXIT_FAILURE);
//        }
//        strain_ = strain;
//        // strain rate has a unit of [1/s]
//        strainRate_ = strainRate;
//    }
    }
    
    void printTime() const override
    {
        logger(INFO, "t=%, tMax=%, N=%, kE=%, bDx% =%, volume fraction =%", getTime(), getTimeMax(), particleHandler
                       .getSize(), particleHandler.getKineticEnergy() / getElasticEnergy(), scaling_, particleHandler.getMass()
                                                                                                      / cylVol_ * scaling_,
               particleHandler.getVolume() / cylVol_);
    }
    
    void actionsBeforeTimeStep() override
    {
        if (particleHandler.getKineticEnergy() / getElasticEnergy() < threshold_)
        {
            Mdouble bulkDensity = particleHandler.getMass() / cylVol_ * scaling_;
            Mdouble trueDensity = particleHandler.getMass() / particleHandler.getVolume();
            logger(INFO, "Bulk density = %, True Density = %, volume fraction = %", bulkDensity, trueDensity,
                   particleHandler.getVolume() / cylVol_);
            deleteParticles();
            setTimeMax(getTime());
        }
    }
    
    void deleteParticles()
    {
        CubeDeletionBoundary cubedeletionBoundary;
        cubedeletionBoundary.set({-cylRad_, -cylRad_, cylHeight_}, {cylRad_, cylRad_, 2 * cylHeight_});
        boundaryHandler.copyAndAddObject(cubedeletionBoundary);
        checkInteractionWithBoundaries();
        logger(INFO, "Boundary %: CubedeletionBoundary", boundaryHandler.getSize() - 1);
    }
    
    
    void writeEneTimeStep(std::ostream& os) const override
    {
        DPMBase::writeEneTimeStep(os);
        Mdouble bulkDensity = particleHandler.getMass() / cylVol_ * scaling_;
        os << " " << bulkDensity << std::endl;
    }

private:
    // Cylinder height
    Mdouble cylHeight_;
    // Cylinder radius
    Mdouble cylRad_;
    // Cylinder volume
    Mdouble cylVol_;
    // threshold to end the simulation and here also delete particles above a certain limit
    Mdouble threshold_ = 1e-6;
    // scaling = 4 for quarter partitioning of simulation domain and scaling = 2 if only the half of the domain gets
    // computed
    int scaling_ = 1;
};


int main(int argc, char *argv[])
{
    //The following parameters need to be inserted after generating the executable.
    //- Generate the executable: make CylinderParticleInsertion
    //- Run the executable: ./CylinderParticleInsertion 1200
    // runs executable with a particle density equal to 1200
    
    // Particle Density
    Mdouble pDensity = atof(argv[1]);
    Mdouble cScaling = atof(argv[2]);
    // Set up a problem of type CylinderParticleInsertion
    CylinderParticleInsertion CPI(pDensity, cScaling);
    // Set simulation time
    CPI.setTimeMax(10);
    CPI.cg = CPI.cgHandler.copyAndAddObject(CG<CGCoordinates::O>());
    // start the solver
    CPI.solve();
    return 0;
}
