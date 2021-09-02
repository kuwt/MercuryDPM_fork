#include <Mercury3D.h>
#include <Walls/InfiniteWall.h>
#include <Species/LinearPlasticViscoelasticSpecies.h>
#include <Boundaries/PeriodicBoundary.h>

class UniaxialCompression : public Mercury3D{
    
    InfiniteWall* wTop;
    double strain_ =0.0;
    double strainRate_=1.0;
    double loadingTime =0;

private:
    void setupInitialConditions() override{
        setName("UniaxialCompression");
        setGravity(Vec3D(0.0,0.0,-9.81));
        LinearPlasticViscoelasticSpecies* s1;
        s1 = dynamic_cast<LinearPlasticViscoelasticSpecies*>(speciesHandler.getObject(0));
        s1->setPlasticParameters(258.5,8*258.5
                ,258.5,0.1);
        s1->setDissipation(0.00);
        s1->setDensity(300);
//        LinearPlasticViscoelasticSpecies* s2;
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
        Mdouble Vp = 0;
        double N = particleHandler.getNumberOfObjects();
        for (int i=0; i < N; i++) {
            //particleHandler.getObject(i)->setVelocity(Vec3D(0.0, 0.0, 0.0));
            particleHandler.getObject(i)->setSpecies(s1);
            Vp = Vp+particleHandler.getObject(i)->getVolume();
        }
        /// periodic Boundary ///
        PeriodicBoundary* PB;
        // y-axis
        PB = dynamic_cast<PeriodicBoundary*>(boundaryHandler.getObject(1));
        PB->set(Vec3D(0,1,0),getYMin(),getYMax());
        // x-axis
        PB = dynamic_cast<PeriodicBoundary*>(boundaryHandler.getObject(2));
        PB->set(Vec3D(1,0,0),getXMin(),getXMax());
        /// Infinite wall ///
        //z-axis bottom
        InfiniteWall* w;
        w = dynamic_cast<InfiniteWall*>(wallHandler.getObject(0));
        w->set(Vec3D(0,0,-1),Vec3D(0,0,getZMin()));
        //z-axis top
        wTop->set(Vec3D(0,0,1),Vec3D(0,0,getZMax()));
        //const Mdouble Lz0 = (getZMax()-getZMin());
        // final (maximum) top-wall position to satisfy the given strain [loading stage].
        //const Mdouble finalPos = Lz0*(1-strain_);
        
        const Mdouble l0 = (getZMax()-getZMin());
        const Mdouble finalPos = l0*(1-strain_);
        // time required for the wall to cross finalPos distance
        //double loadingTime = strain_ /strainRate_;
        // time required for the wall to cross finalPos distance
        loadingTime = -log(finalPos/l0)/strainRate_;
        // reset Time Max
        // check if the time max is sufficient for both loading and unloading stages
        logger(INFO, "Maximum time for Loading and unloading stages is twice the loading time(approximately)!");
        setTimeMax(2 * loadingTime);
        std::cout << "Time max: "<< loadingTime <<std::endl;
        std::cout << "Sample Volume: "<< Vp <<std::endl;
        // this for setting numbering to start form zero
        getVtkWriter()->setFileCounter(0);
    }
    int countDown = 1; // counter for loading
    int countUp = 1;   // counter for un-loading
    void actionsAfterTimeStep() override {
        // initial length
        Mdouble l0 = (getZMax() - getZMin());
        Mdouble finalPos = l0*(1-strain_);
        // ALWAYS PUT THIS COMMANDS FIRST
        // write a file for stress strain on the top wall
        /// check! saveCount should be the same here
        if (getNumberOfTimeSteps()%2000==0) {
            static std::ofstream myFile("strainStressTopWall.dat");
            myFile << getTime()
                   << '\t' << wallHandler.getObject(1)->getPosition().Z
                   << '\t' << (l0-wallHandler.getObject(1)->getPosition().Z)/l0
                   << '\t' << wallHandler.getObject(1)->getForce()/(getXMax()*getYMax())
                   << '\n';
        }
        //double loadingTime = -log(finalPos/l0)/strainRate_;
        //double loadingTime = strain_ /strainRate_;
        // if the position of the wall now is bigger than finalPos do nothing
        if (getTime() >= loadingTime) {// UNLOADING STAGE
            // here l0 = finalPos
            countUp++;
            // Calculate the next position of the wall
            Mdouble le = finalPos * exp(+strainRate_ * countUp * getTimeStep());
            /*linear calculation*/
            //Mdouble l = l0*(1-strainRate_*countUp*getTimeStep());
            // set the new location of the wall
            wallHandler.getObject(1)->setPosition(Vec3D(0, 0, le));
            //termination criterion >> when the force in z-direction = 1e-6
            if (std::abs(wallHandler.getObject(1)->getForce().Z)<=0){
                setTimeMax(getTime());
            }
            return;
        }
            /// LOADING STAGE ///
        else if (finalPos < wallHandler.getObject(1)->getPosition().Z){
            // why we do not use countDown here ? because time here already start from 1
            // Calculate the next position of the wall
            /*logarithmic calculation*/
            Mdouble l = l0*exp(-strainRate_*getTime());
            /*linear calculation*/
            //Mdouble l = l0*(1-strainRate_*countDown*getTimeStep());
            // set the new location of the wall
            wallHandler.getObject(1)->setPosition(Vec3D(0,0,l));
            return;
        }
        
    }

//public:
//    // constructor
//    explicit UniaxialComp(std::string& fileName){
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
};

int main()
{
    // Set up a problem of type MarbleRun
    UniaxialCompression UC;
    // Set simulation time
    UC.setTimeMax(1);
    // start the solver
    UC.solve();
}