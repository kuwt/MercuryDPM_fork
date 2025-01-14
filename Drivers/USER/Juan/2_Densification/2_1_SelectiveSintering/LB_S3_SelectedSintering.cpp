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

#include "Mercury3D.h"
#include "Species/ThermalSinterLinFrictionReversibleAdhesiveSpecies.h"
#include "Walls/InfiniteWall.h"
//#include "CG/CG.h"

/* This code reads the information generated by S2_Compression file.
 * The temperature increases and heats up all particles. Thus, the Sintering process
 * occurs.
*/
// Main class:
class LB_SelectedSintering : public Mercury3D{

private:
    Mdouble setTConductivity_ = 0.0;
    Mdouble setHCapacity_ = 0.0;

    Mdouble T0_ = 0.0;
    Mdouble maxTemp_ = 0.0;
    Mdouble gradientTemp_ = 0.0; //[C/s] Heating rate
    Mdouble holdingTime_ = 0.0;
    Mdouble meltingTemperature_ = 0.0;

    Mdouble meanCoordinationNumber = 0.0;
    Mdouble deltaR_ = 0.0; //expansion coefficient

    Mdouble scalarNormalForce = 0.0;

    Mdouble newHeight = 0.0;

    Mdouble setDeltaC_ = 0.0;
    Mdouble setC1_ = 0.0;

    Vec3D setLaserPosition_;

    ThermalSinterLinFrictionReversibleAdhesiveSpecies* particleSpecies;
    ThermalSinterLinFrictionReversibleAdhesiveSpecies* particleSpeciesAtMeltingPoint;

public:
    //Constructor
    LB_SelectedSintering(Mdouble deltaR, Mdouble startingTemperature, Mdouble  maxTemp, Mdouble gradientTemp,
                            Mdouble meltingTemperature, Mdouble holdingTime, Mdouble setHCapacity, Mdouble setTConductivity,
                            Mdouble setDeltaC, Mdouble setC1)
    {
        deltaR_ = deltaR;

        maxTemp_ = maxTemp;
        T0_ = startingTemperature;
        gradientTemp_ = gradientTemp;
        meltingTemperature_ = meltingTemperature;
        holdingTime_ = holdingTime;

        setTConductivity_ = setTConductivity;
        setHCapacity_ = setHCapacity;

        setDeltaC_ = setDeltaC;
        setC1_ = setC1;

        //Read file from S2_Compressione
        setName("LB_Compression150L");
        readRestartFile();
        //Set output file
        setRestarted(false);
        setName("LB_SelectedSintering150L");

        //Remove Lid from S2_Compression
        wallHandler.removeLastObject();

        //-------------------------------------------------
        //Get the particle species to set extra parameters to the particles
        particleSpecies = dynamic_cast<ThermalSinterLinFrictionReversibleAdhesiveSpecies*>(speciesHandler.getObject(0));
        //-------------------------------------------------

        const Mdouble YoungM = 1.65e9; //[Pa] Young's Modulus for polyamide12
        Mdouble K1 =  particleSpecies->getLoadingStiffness();

        //Set thermal parameters
        //Thermal parameters:
        particleSpecies->setThermalConductivity(setTConductivity_);
        particleSpecies->setHeatCapacity(setHCapacity_);

        particleSpecies->setSinterType(SINTER_APPROACH::VISCOELASTIC_CONTACT);  //FRENKEL OR VISCOELASTIC_CONTACT

//        particleSpecies->setCohesionStiffness(K1);
        particleSpecies->setPenetrationDepthMax(1.4);
        particleSpecies->setComplianceZero(1.0/(2.0*YoungM)); //Book: Adhesive Particles
        particleSpecies->setSurfTension(0.03); //

        //To control:
//        particleSpecies->setSinterAdhesion(1.0e-8/particleHandler.getMeanRadius());
//        particleSpecies->setSinterAdhesion(0.001*K1);
        particleSpecies->setSinterAdhesion(0.0);

        particleSpecies->setFluidity(setC1_);
        particleSpecies->setSeparationDis(setDeltaC_);

        //Set a "virtual" species to allow change particle properties close to the melting point
        particleSpeciesAtMeltingPoint = new ThermalSinterLinFrictionReversibleAdhesiveSpecies(*particleSpecies);
    }
    //--------------------------------------------------
    //Destructure
    ~LB_SelectedSintering() override
    {
        delete particleSpeciesAtMeltingPoint;
    }
    //--------------------------------------------------
    //Initial conditions
    void setupInitialConditions() override
    {
        particleSpeciesAtMeltingPoint->setLoadingStiffness(particleSpecies->getLoadingStiffness());
        particleSpeciesAtMeltingPoint->setSinterAdhesion(0.0016*particleSpecies->getLoadingStiffness());
        // Set uniform temperature to all particles.
        for (BaseParticle* p : particleHandler)
        {
            auto* tp = dynamic_cast<ThermalParticle*>(p);
            tp->setTemperature(T0_);
            tp->setVelocity(Vec3D(0.0,0.0,0.0));
        }
        logger(INFO, "Set uniform temperature to all particles [C]= %", T0_);

        if(holdingTime_>getTimeMax()){

            logger(INFO, "*******************************");
            logger(INFO, "Set correctly the holding time");
            logger(INFO, "*******************************");
            exit(-1);
        }
    }
    //--------------------------------------------------
    Mdouble getSurfaceHeight() const
    {
        double height = 0;

        for (const BaseParticle* p : particleHandler)
        {
            double newHeight = p->getPosition().Z;
            if (height<newHeight) height = newHeight;
        }
        return height+particleHandler.getMeanRadius();
    }

    //--------------------------------------------------
    Mdouble getTemperature()const
    {
        Mdouble heatingTemperature = T0_ + gradientTemp_*getTime();
        Mdouble coolingTemperature = maxTemp_;

        if(getTime()>=getTimeMax()-holdingTime_){
            //Cooling gradient is less than heating gradient.
            coolingTemperature = T0_ + (gradientTemp_*0.7)*(getTimeMax() - getTime());
        }
        return std::min(std::min(heatingTemperature,coolingTemperature),maxTemp_);
    }

    void setLaserPosition()
    {
        //This vector will change eventually.
        //Now, the position is at the center of the domain
        Vec3D getCenter;
        getCenter.X = (getXMax()-getXMin())/2.0;
        getCenter.Y = (getYMax()-getYMin())/2.0;
        getCenter.Z = (getZMax()-getZMin())/2.0;

        setLaserPosition_.X =  getCenter.X;
        setLaserPosition_.Y = getCenter.Y;
        setLaserPosition_.Z = getSurfaceHeight();
    }

    //--------------------------------------------------
    void setTemperature(Mdouble temperature)
    {
        //Variables to update
        static Mdouble oldTemperature = T0_;
        Mdouble deltaTemperature = temperature-oldTemperature;
        Mdouble deltaCooling = temperature - meltingTemperature_;

        //add temperature above a certain height
        Mdouble lb_z = getSurfaceHeight() - 2.0 * particleHandler.getMeanRadius();
        Mdouble lb_x = setLaserPosition_.X;
        Mdouble lb_y = setLaserPosition_.Y;

        for (BaseParticle* p : particleHandler)
        {
            auto *tp = dynamic_cast<ThermalParticle *>(p);

            if (p->getSpecies()==particleSpecies)
            {
                if (p->getPosition().Z > lb_z and (p->getPosition().X > lb_x*0.9 and p->getPosition().X < lb_x*1.1)
                    and (p->getPosition().Y > lb_y*0.9 and p->getPosition().Y < lb_y*1.1))
                {
                    tp->setTemperature(temperature); //Change the temperature of the particle
                }
            }
        }
        //Set changes of particles heated close to the melting point
        setHeatEffect(temperature,deltaTemperature,deltaCooling);

        //Update temperature
        oldTemperature = temperature;
    }

    //Function to change particle properties based on the temperature
    void setHeatEffect(Mdouble temperature,Mdouble deltaTemperature, Mdouble deltaCooling)
    {
        //Thermal expansion
        Mdouble factorRadius = 1.0 + (deltaR_ * deltaTemperature);

        //change species properties
        Mdouble stiffnessFactor = 0.5*(1.0+tanh((meltingTemperature_-temperature)/gradientTemp_));
        Mdouble oldLoadingStiffness = particleSpecies->getLoadingStiffness();

        for (BaseParticle* p : particleHandler)
        {
            if (p->getSpecies()==particleSpecies)
            {
                auto* tp = dynamic_cast<ThermalParticle*>(p);
                Mdouble pT = tp->getTemperature();
                //Check which particles are 95% closer to the melting point
                if(pT >= 0.95*meltingTemperature_)
                {
                    p->setRadius((factorRadius * p->getRadius()));
                    particleSpecies->setLoadingStiffness(stiffnessFactor*particleSpeciesAtMeltingPoint->getLoadingStiffness());
                    particleSpecies->setSinterAdhesion(particleSpeciesAtMeltingPoint->getSinterAdhesion());
                    particleHandler.getObject(p->getIndex())->setSpecies(particleSpecies);

                }else
                {
                    particleSpecies->setLoadingStiffness(particleSpeciesAtMeltingPoint->getLoadingStiffness());
                    particleSpecies->setSinterAdhesion(0.01*particleSpeciesAtMeltingPoint->getSinterAdhesion());
                    particleHandler.getObject(p->getIndex())->setSpecies(particleSpecies);
                }
            }
        }


        //for decreasing temperature, change the maxOverlap
        if (deltaCooling<0.0)
        {
            for (BaseInteraction* cBase : interactionHandler)
            {
                auto c =  dynamic_cast<SinterLinInteraction*>(cBase);
                Mdouble unloadingStiffness = c->getUnloadingStiffness();
                c->setMaxOverlap(c->getMaxOverlap()
                                 *(unloadingStiffness-oldLoadingStiffness)
                                 /(unloadingStiffness-particleSpecies->getLoadingStiffness())
                );
            }
        }
    }
    //--------------------------------------------------
    void actionsAfterTimeStep() override
    {
        setLaserPosition();
        setTemperature(getTemperature());
        newHeight = getSurfaceHeight();


        //To compute the normal force
        for (auto i : interactionHandler)
        {
            scalarNormalForce += Vec3D::dot(i->getForce(),i->getNormal());
        }

        //To measure the mean coordination number.
        for (int i = particleHandler.getNumberOfObjects()-1; i>=0; i--)
        {
            meanCoordinationNumber += (particleHandler.getObject(i) -> getInteractions()).size();
        }
        meanCoordinationNumber /= particleHandler.getNumberOfObjects();
    }

    //System
    Mdouble getThermalEnergy() const
    {
        Mdouble eneEnergy;

        for (BaseParticle* p : particleHandler)
        {
            if (p->isFixed()) continue;

            ThermalParticle* tp = dynamic_cast<ThermalParticle*>(p);
            eneEnergy += (tp->getTemperature()-T0_) * tp->getMass() * particleSpecies->getHeatCapacity();
//            eneEnergy += deltaTemperature * tp->getMass() * particleSpecies->getHeatCapacity();
        }
        return eneEnergy;
    }

    //--------------------------------------------------
    //Prints the temperature into the data file, such that you can plot it in paraview
    Mdouble getInfo(const BaseParticle& p) const override
    {
        return getTemperature();
//        return (dynamic_cast<const ThermalParticle &>(p).getTemperature()-T0_)/(temperatureSurface_-T0_);
    }

    Mdouble getMeanRelativeContactRadius() const
    {
        Mdouble meanOverlap = interactionHandler.getMeanOverlap();
        Mdouble meanRadius = particleHandler.getMeanRadius();
        return sqrt(meanOverlap/meanRadius); //Todo: Check if it is mean Radius or initial radius
    }
    //--------------------------------------------------
    //Function to override the output file with specific parameters
    void writeFstatHeader(std::ostream &os) const override
    {
        Mdouble volSystem = constants::pi*(std::pow(getXMax()/2,2))*getZMax();

        Mdouble volParticle = (4.0/3.0)*constants::pi*(std::pow(particleHandler.getMeanRadius(),3.0));
        Mdouble densityParticle = speciesHandler.getObject(0)->getDensity();

        Mdouble volTotalParticles = particleHandler.getVolume();
        Mdouble volumeFraction = volTotalParticles/volSystem;

        Mdouble massParticle = densityParticle*volParticle;
        Mdouble massTotalParticles = particleHandler.getMass();

        Mdouble meanRadius = particleHandler.getMeanRadius();
        Mdouble meanOverlap = interactionHandler.getMeanOverlap();
        Mdouble meanContactRadius = sqrt(meanOverlap/meanRadius)*meanRadius;

        os << getTime() //[1] Current time
           << " " << getTimeMax() //[2] Max time
           << " " << volSystem //[3] Volume of the system
           << " " << getZMax() //[4] Height of the system
           << " " << particleHandler.getNumberOfObjects() //[5] Number of particles inserted
           << " " << particleHandler.getMeanRadius() // [6] Mean radius of a particle
           << " " << volParticle //[7] Mean volume of a particle
           << " " << densityParticle //[8] Mean density of a particle
           << " " << volTotalParticles //[9] Volume of particles inserted
           << " " << volumeFraction //[10] Volume fraction
           << " " << massTotalParticles //[11] total mass of particle inserted
           << " " << getSurfaceHeight() // [12] Max height reached by particles.
           << " " << getKineticEnergy() // [13] kinetic energy
           << " " << getElasticEnergy() //[14] Elastic energy
           << " " << scalarNormalForce //[15] Normal force
           << " " << meanCoordinationNumber // [16] meanCoordinationNumber
           << " " << getTemperature()//[17]Current temperature
           << " " << particleSpecies->getLoadingStiffness()// [18] Loading stiffness
           << " " << particleSpeciesAtMeltingPoint->getLoadingStiffness()// [19] Loading stiffness
           << " " << getMeanRelativeContactRadius()//[20] Neckgrowth
           << " " << meanOverlap//[21]
           << " " << meanContactRadius //[22]
           << " " << getThermalEnergy() //[23]
           << std::endl;
    }
    //--------------------------------------------------
    //Function to display data at the console
    void printTime() const override
    {
        Mdouble volSystem = constants::pi*(std::pow(getXMax()/2,2))*getZMax();
        Mdouble volTotalParticles = particleHandler.getVolume();

        std::cout << "t=" << std::setprecision(3)<< std::left << std::setw(3) << std::left<< getTime()
        << ", tmax= " <<  std::setprecision(3) << std::left << std::setw(3) << getTimeMax()
        << ", Temp=" << std::setw(3)<< std::left << std::setw(3) << getTemperature()
        << ", KineticEnergy=" << std::setprecision(3) << std::left << std::setw(6) << getKineticEnergy()
        << ", ThermalEnergy="<< std::setw(3)<< std::left << std::setw(3) << getThermalEnergy()
        << " loading stiffness= " << std::setw(3)<< std::left << std::setw(3) << particleSpecies->getLoadingStiffness()
        << " MeanOverlap= "<< std::setw(3)<< std::left << std::setw(3) << getMeanRelativeContactRadius()
        << std::endl;
    }
};

// Main function:
int main(int argc UNUSED, char *argv[] UNUSED)
{
    Mdouble deltaR = 0.0002; //Expansion coefficient close to the melting point
    Mdouble meltingTemperature = 183.0; // melting point of the material

    Mdouble startingTemperature = 170.0; //Initial temperature
    Mdouble gradientTemp = 840.0; //[C/s] 10.0//Todo:Check!

    Mdouble maxTemp = 183.0; //max temperature

    Mdouble holdingTime = 0.18; //[s] After this time temperature starts decreasing

    Mdouble setHCapacity = 2500.0; //[J/(kgK)]
    Mdouble setTConductivity = 120.25;  //0.25[W/(mK)] //Todo:Check

    Mdouble setDeltaC = 9.55e-07;
    Mdouble setC1 = 1.26;

    LB_SelectedSintering oTest(deltaR, startingTemperature, maxTemp, gradientTemp, meltingTemperature, holdingTime, setHCapacity, setTConductivity,
                                  setDeltaC, setC1);

    oTest.setParticlesWriteVTK(false);

    oTest.setTimeMax(0.3); //[s]
    oTest.setXBallsAdditionalArguments("-solidf -v0 -cmode 8 -cmaxset 100 ");

    //--------------------------------------------------
    //This is for Coarse-Graining live.
    //Live CG to get results over time. For instance, plot stress vs time.
    //    //define coarse-graining resolved in z
    //    CG<CGCoordinates::Z> cgZ; //declare a cg object
    //    cgZ.setN(20);  //set number of mesh points in each spatial direction
    //    cgZ.setWidth(0.5); //set cg width
    //    s.cgHandler.copyAndAddObject(cgZ); // add the CG object to the cgHandler

    oTest.removeOldFiles();
    oTest.solve();

    //Coarse-Graining:
    //--------------------------------------------------
    // It creates the coarse-graining output file at the last iteraction.
    logger(INFO,"Execute 'source S3_Sintering.sh' to get coarse-grained statistics at specific time step");
    helpers::writeToFile("S3_Sintering.sh","../../../../../MercuryCG/fstatistics S3_Sintering -stattype XZ -w 1.0e-5 -h 0.1e-5 -tmin 1.0 -tmax 1.1 -o S3_SinteringScaledMass.XZ.stat");

    //--------------------------------------------------
    // It creates the matlab visualization.
    logger(INFO,"Run 'S3_Sintering.m' in MATLAB/octave to visualise the statistical output");
    helpers::writeToFile("S3_Sintering.m","clc;clear all;close all\n"
                                                 "addpath('../../../../../MercuryCG/')\n"
                                                 "data = loadStatistics('S3_Sintering.XZ.stat');\n"
                                                 "colormap(1-gray)\n"
                                                 "contourf(data.x,data.z,data.Density,20,'EdgeColor','none')\n"
                                                 "c = colorbar\n"
                                                 "c.Label.String = '\\rho';\n"
                                                 "title('Density')\n"
                                                 "xlabel('x')\n"
                                                 "ylabel('z');\n"
                                                 "axis equal\n"
                                                 "%%\n"
                                                 "addpath('/home/juan/Softwares/MercuryDPM_Branch/Trunk/Matlab');\n"
                                                 "particles=read_data('S3_Sintering.data');\n"
                                                 "Cell = particles{637}(1,1); %Specific particle position at 1.119, which is in the cell 24\n"
                                                 "NumPart= Cell.N;\n"
                                                 "Pradius = Cell.Radius;\n"
                                                 "Position = Cell.Position;\n"
                                                 "a=linspace(0,2*pi,40);\n"
                                                 "xCircle = sin(a);\n"
                                                 "zCircle = cos(a);\n"
                                                 "hold on;\n"
                                                 "for i=1:length(Pradius)\n"
                                                 "  plot(Position(i,1) + Pradius(i)*xCircle,Position(i,3)+Pradius(i)*zCircle,'Color',.8*[1 1 1])\n"
                                                 "end\n"
                                                 "hold off");
    return 0;
}