//Copyright (c) 2013-2020, The MercuryDPM Developers Team. All rights reserved.
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
#include <random>
#include "Mercury3D.h"
#include "Walls/InfiniteWall.h"
#include <Species/LinearViscoelasticFrictionSpecies.h>
#include "Particles/SphericalParticle.h"
#include "Boundaries/CubeInsertionBoundary.h"
#include "DPMBase.h"

//This selftTest tests if the setDistributionPhiNormal() function gives the correct D50 and PSD as demanded.
//It outputs the final D50 and the inputted D50, as well as the PSD with paired radius and probability.
//With the outputted PSD, users can check the standard deviation by other tools, e.g., matlab.
class SetDistributionPhiNormalSelfTest : public Mercury3D {

private:
    double D50_, StdDev_;
    int numofPSD_;

public:
    void setupInitialConditions () override {


        auto species = speciesHandler.copyAndAddObject(LinearViscoelasticFrictionSpecies());
        species->setDensity(2650);
        species->setStiffness(1500);
        species->setDissipation(0.002);
        species->setSlidingStiffness(0.29*species->getStiffness());
        species->setSlidingDissipation(0.29*species->getDissipation());
        species->setSlidingFrictionCoefficient(0.4);
        species->setSlidingFrictionCoefficientStatic(0.4);
        species->setRollingStiffness(0.4*species->getStiffness());
        species->setRollingDissipation(0.4*species->getDissipation());
        species->setRollingFrictionCoefficient(0.5);

        setSaveCount(1000);
        setTimeMax(0.001);
        setTimeStep(1e-6);

        SphericalParticle p0;
        p0.setSpecies(speciesHandler.getObject(0));

        CubeInsertionBoundary insertionBoundary;
        unsigned maxFail = 1; //insert as quick as possible: try every time step, until you maxFail=1 particle fails to be insertable (overlaps with another particle or wall)
        Vec3D posMin = {getXMin(),getYMin(),getZMin()};
        Vec3D posMax = {getXMax(),getYMax(),getZMax()};
        Vec3D velMin = {-0.005, 0, -0.005};
        Vec3D velMax = {0.005, 0, 0};
        insertionBoundary.set(&p0, maxFail, posMin, posMax, velMin, velMax);
        //insertionBoundary.setCheckParticleForInteraction(false);

        PSD psd;
        psd.setDistributionPhiNormal(D50_, StdDev_, numofPSD_);
        logger(INFO, "Final distribution");
        for (auto p: psd.getParticleSizeDistribution()){
            logger(INFO, "%, %", p.internalVariable, p.probability);
        }

        insertionBoundary.setPSD(psd);
        insertionBoundary.setInitialVolume((getXMax()-getXMin())*(getYMax()-getYMin())*(getZMax()-getZMin())*0.25);
        boundaryHandler.copyAndAddObject(insertionBoundary);

        InfiniteWall w0;
        w0.setSpecies(speciesHandler.getObject(0));
        w0.set(Vec3D(0.0, 0.0, -1.0), Vec3D(0.0, 0.0, getZMin()));
        wallHandler.copyAndAddObject(w0);
        w0.set(Vec3D(0.0, -1.0, 0.0), Vec3D(0.0, getYMin(), 0.0));
        wallHandler.copyAndAddObject(w0);
        w0.set(Vec3D(0.0, 1.0, 0.0), Vec3D(0.0, getYMax(), 0.0));
        wallHandler.copyAndAddObject(w0);
        w0.set(Vec3D(-1.0, 0.0, 0.0), Vec3D(getXMin(), 0.0, 0.0));
        wallHandler.copyAndAddObject(w0);
        w0.set(Vec3D(1.0, 0.0, 0.0), Vec3D(getXMax(), 0.0, 0.0));
        wallHandler.copyAndAddObject(w0);
    };

    /** Note: std is the standard deviation of diameter in Phi Units **/
    void setD50StandardDeviationAndNumofPSD(double D50, double std, int num){
        D50_ = D50;
        StdDev_ = std;
        numofPSD_ = num;
    }

    void printTime() const override {
        logger(INFO, "t=%3.6", getTime());
    }
};

int main(int argc UNUSED, char *argv[] UNUSED){

    SetDistributionPhiNormalSelfTest st;
    st.setName("SetDistributionPhiNormalSelfTest");
    st.setMax({0.005,0.005,0.005});
    st.setGravity({0.0,0.0,-9.8});
    st.setD50StandardDeviationAndNumofPSD(0.00025, 0.5, 100);

    st.setFileType(FileType::ONE_FILE);
    st.solve();

}
