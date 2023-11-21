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


#include <CMakeDefinitions.h>
#include "Mercury3D.h"
#include "Species/LinearViscoelasticSpecies.h"
#include "Walls/InfiniteWall.h"
#include "Boundaries/CylinderInsertionBoundary.h"


/**
 * An example code on how to use the CylinderInsertionBoundary as well as shift and rotate it in cartesian space
 */

class PSDSelfTest : public Mercury3D
{
public:

    /**
     * Use setupInitialConditions to define a.o. an insertion boundary, which will insert particles in a given region, at a given flow rate and of a given particle size distribution
     */
    void setupInitialConditions() override
    {
        setName("CylinderInsertionBoundaryUnitTest");
        setSystemDimensions(3);
        setGravity(Vec3D(0, 0, 0));
        setTimeStep(1e-4);
        dataFile.setSaveCount(10);
        setTimeMax(2e-2);
        setHGridMaxLevels(2);

        setMin(Vec3D(0, 0, 0));
        setMax(Vec3D(.01, .01, .01));

        setParticlesWriteVTK(true);
        LinearViscoelasticSpecies species;
        species.setDensity(2000);
        species.setStiffness(10000);
        speciesHandler.copyAndAddObject(species);

        // define a cube insertion boundary, which inserts particles in a cylindrical region
        CylinderInsertionBoundary insertionBoundary;

        // define the type of particle you want to insert (radius, position and velocity will be set by the insertion boundary)
        SphericalParticle templateParticle(speciesHandler.getObject(0));

        // insert particles in the whole domain (between getMin and getMax) with initial velocity 0 and a radius between 1 and 2 (uniform number distribution)
        unsigned maxFail = 1; //insert as quick as possible: try every time step, until you maxFail=1 particle fails to be insertable (overlaps with another particle or wall)
        // shift each coordinate by 1
        Vec3D shift = {1,1,1};
        // rotate around X by 90 degrees, then around Y by 90 degrees and then around Z by 90 degrees
        Vec3D rotation = {constants::pi/2.0, constants::pi/2.0, constants::pi/2.0};
        // set the normal to be the x-axis (will be the height coordinate of the cylinder)
        insertionBoundary.set(&templateParticle, maxFail, getXMin(), getXMax(), getZMin(), getZMax(), {1,0,0});
        helpers::check(insertionBoundary.getOrientationMatrix(), Matrix3D(0,0,1,0,1,0,-1,0,0), 1e-3, "aligning the cylinder to its normal direction failed");
        insertionBoundary.shiftBoundary(shift);
        // check shift of the boundary
        logger.assert_always(insertionBoundary.getOrigin() == Vec3D(1,1,1), "The shift of the CylinderInsertionBoundary is not correct");
        // with this rotation the Z-axis will end up being the height coordinate of the cylinder in a right hand system
        insertionBoundary.rotateBoundary(rotation);
        // check rotation of the boundary
        helpers::check(insertionBoundary.getOrientationMatrix(), Matrix3D(-1,0,0,0,1,0,0,0,-1), 1e-3, "The rotation matrix of the CylinderInsertionBoundary is not correctly set");
        PSD psd;
        // laser diffraction measurements usually give volumetric responses. This PSD was already converted to a
        // CUMULATIVE_NUMBER_DISTRIBUTION.
        psd.setPSDFromCSV(getMercuryDPMSourceDir() + "/Drivers/SelfTests/Boundaries/InputData/APAPM.csv",
                          PSD::TYPE::CUMULATIVE_NUMBER_DISTRIBUTION);
        insertionBoundary.setPSD(psd);

        // add the insertion boundary to the handler
        boundaryHandler.copyAndAddObject(insertionBoundary);
    }

    void printTime() const override
    {
        logger(INFO, "t=%, tMax=%, N=%", getTime(), getTimeMax(), particleHandler.getSize());
    }

};

int main(int argc UNUSED, char* argv[] UNUSED)
{
    PSDSelfTest problem;
    problem.solve();

    // check if the particles where inserted correctly into the CylinderInsertionBoundary
    logger.assert_always(problem.particleHandler.getSize() == 209, "The number of particles in the CylinderInsertionBoundary is not correct");
}
