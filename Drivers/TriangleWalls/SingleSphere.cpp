//Copyright (c) 2013-2018, The MercuryDPM Developers Team. All rights reserved.
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

#include "MercuryProb.h"

void MercuryProblem::setupInitialConditions() {
/// define particle positions and velocities
    SphericalParticle p0;
    p0.setSpecies(speciesHandler.getLastObject());
    p0.setRadius(radius);
    static Mdouble kn = dynamic_cast<const LinearViscoelasticSpecies *>(speciesHandler.getObject(
            0))->getStiffness();
    posZ0 = RVESize + p0.getRadius() + p0.getMass() * getGravity().getZ() / kn;
    p0.setPosition(Vec3D(0.5 * (getXMax() + getXMin()), 0.5 * (getYMax() + getYMin()), posZ0));
    p0.setVelocity(Vec3D(0.0, -0.1, 0.0));
    particleHandler.copyAndAddObject(p0);
    p0.setPosition(
            Vec3D(0.5 * (getXMax() + getXMin()) - 0.5 * RVESize, 0.5 * (getYMax() + getYMin()) - 0.25 * RVESize,
                  posZ0));
    p0.setVelocity(Vec3D(0.0, -0.1, 0.0));
    particleHandler.copyAndAddObject(p0);
    createWalls();
}

void MercuryProblem::createWalls() {
    int nFacesX = 2;
    int nFacesY = 5;
    for (int i = 0; i < nFacesX; i++) {
        for (int j = -5; j < nFacesY; j++) {
            std::vector<Vec3D> vertices;
            vertices.push_back(Vec3D((-1 + i) * RVESize, j * RVESize, RVESize));
            vertices.push_back(Vec3D(i * RVESize, j * RVESize, RVESize));
            vertices.push_back(Vec3D(i * RVESize, (j + 1) * RVESize, RVESize));
            vertices.push_back(Vec3D((-1 + i) * RVESize, (j + 1) * RVESize, RVESize));
            Vec3D center = Vec3D((-0.5 + i) * RVESize, (j + 0.5) * RVESize, RVESize);

            unsigned n = 0;
            // create TriangleWalls from oomph face element
            while (n < 4) {
                // get vertices of TriangleWall (multiply vertex position with the length scale of the OomphProblem)
                std::array<Vec3D, 3> vertex;
                // one vertex at the center
                vertex[0] = center;
                // two vertices from the OomphProblem<ELEMENT,TIMESTEPPER>
                vertex[1] = vertices[0];
                vertex[2] = vertices[1];

                // create triangle facet
                TriangleWall *w = createTriangleWall(vertex);
                w->setGroupId(10);
                listOfTriangleWalls.push_back(w);

                // rotate forward by one element
                std::rotate(vertices.begin(), vertices.begin() + 1, vertices.end());
                n++;
            }
        }
    }
}

int main(int argc, char *argv[]) {

    /// RVE size
    double RVESize = 0.4 / units::lUnit;

    /// Setup the coupled problem
    MercuryProblem TriangleWallsTest;
    TriangleWallsTest.setupMercuryProblem("SingleSphere", 3, RVESize, 1, -9.81*100 / units::accUnit(), 20 / units::tUnit, 1000);

    TriangleWallsTest.solve();

} //end of main
