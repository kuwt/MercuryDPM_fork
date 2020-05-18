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

//based on /storage2/usr/people/sluding/MDCC/C3DshearXL30/MU0_LONG2
#include <Species/LinearPlasticViscoelasticSlidingFrictionSpecies.h>
#include <Particles/BaseParticle.h>
#include "Mercury3D.h"
#include "Boundaries/AngledPeriodicBoundary.h"
#include "Walls/AxisymmetricIntersectionOfWalls.h"

class ShearCell3DWalls : public Mercury3D {
public:
    ShearCell3DWalls (Mdouble particleRadius, Mdouble outerRadius, Mdouble innerRadius, Mdouble height, Mdouble openingAngle, Mdouble wallThickness)
        : particleRadius_(particleRadius), outerRadius_(outerRadius), innerRadius_(innerRadius), height_(height), openingAngle_(openingAngle), wallThickness_(wallThickness)
    {
        //check input variables
        if (innerRadius_ <0) {
            std::cerr << "error in constructor of ShearCell3D; inner radius has to be positive" << std::endl;
            exit(-1);
        } else if (innerRadius_ > outerRadius_) {
            std::cerr << "error in constructor of ShearCell3D; inner radius has to be smaller than outer radius" << std::endl;
            exit(-1);
        } else if (height_ <0) {
            std::cerr << "error in constructor of ShearCell3D; domain height has to be positive" << std::endl;
            exit(-1);
        } else if (openingAngle_ <0) {
            std::cerr << "error in constructor of ShearCell3D; opening angle has to be positive" << std::endl;
            exit(-1);
        }

        //set default name
        setName("ShearCell3DWalls");
        setXBallsAdditionalArguments("-v0 -solidf -3dturn 0 -w 1180 -s 20");

        //set gravity (extra strong since we have few particle layers)
        setGravity(Vec3D(0.0,0.0,-9.81));

        //use particles of unit mass and diameter
        species =speciesHandler.copyAndAddObject(LinearPlasticViscoelasticSlidingFrictionSpecies());
        species->setDensity(2000);

        //set inter-particle contact properties
        species->setLoadingStiffness(1e2);
        species->setUnloadingStiffnessMax(1.1*1e2);
        species->setCohesionStiffness(0.0);
        species->setPenetrationDepthMax(0.05);
        species->setSlidingStiffness(2.0/7.0*1e2);
        species->setSlidingFrictionCoefficient(0.50);
        species->setDissipation(2e-3*10.0);
        species->setSlidingDissipation(2e-3*10.0);
        setTimeStep(species->computeTimeStep(species->getMassFromRadius(particleRadius)));
        std::cout << "recommended time step: " << getTimeStep() << std::endl;
        setTimeStep(5.0 * getTimeStep());
        std::cout << "used time step for creating initial conditions: " << getTimeStep() << std::endl;
        std::cout << "used dissipation for creating initial conditions: " << species->getDissipation() << std::endl;

        //set time step
        setTimeMax(1e20);
        setSaveCount(10000);

        //set domain
        setXMin(0.0);
        setYMin(0.0);
        setZMin(0.0);
        setXMax(outerRadius_);
        setYMax(outerRadius_*std::sin(openingAngle_));
        setZMax(height_);

        //set boundaries
        //quarter cell
        AngledPeriodicBoundary b;
        Vec3D normal_left(std::sin(openingAngle),-std::cos(openingAngle),0.0);
        Vec3D normal_right(0.0,-1.0,0.0);
        Vec3D origin(0.0,0.0,0.0);
        b.set(normal_left,normal_right,origin);
        boundaryHandler.copyAndAddObject(b);

        //set walls
        AxisymmetricIntersectionOfWalls w;
        Vec3D AxisDirection(0.,0.,1.);
        Vec3D PointOnAxis(0.,0.,0.);
        Vec3D PrismAxis(0.,-1.,0.);
        w.setPosition(PointOnAxis);
        w.setOrientation(AxisDirection);
        w.setSpecies(species);

        //add points in anti-clockwise direction around the prism axis
        std::vector<Vec3D> Points(2);
        Points[0] = Vec3D(innerRadius_,0.0,getZMin());
        Points[1] = Vec3D(innerRadius_,0.0,getZMax());
        w.createOpenPrism(Points,PrismAxis);
        wallHandler.copyAndAddObject(w);

        Points.resize(2);
        Points[0] = Vec3D(outerRadius_,0.0,getZMax());
        Points[1] = Vec3D(outerRadius_,0.0,getZMin());
        w.createOpenPrism(Points,PrismAxis);
        wallHandler.copyAndAddObject(w);

        InfiniteWall iw;
        iw.set(Vec3D(0.0,0.0,-1.0), Vec3D(0.0,0.0,0.0));
        wallHandler.copyAndAddObject(iw);

    }

    bool continueSolve() const
    {
        static unsigned int counter = 0;
        if (++counter>100)
        {
            counter=0;
            if (getKineticEnergy()<9e-5*getElasticEnergy())
                return false;
        }
        return true;
    }

    void printTime() const
    {
        std::cout << "t=" << getTime() << " Ene " << getKineticEnergy()/getElasticEnergy() << std::endl;
    }

    void createWalls(Mdouble particleRadius, Mdouble outerRadius, Mdouble innerRadius, Mdouble height, Mdouble openingAngle)
    {
        particleRadius_ = particleRadius;
        outerRadius_ = outerRadius;
        innerRadius_ = innerRadius;
        height_ = height;
        openingAngle_ = openingAngle;

        Mdouble dHeight = getZMax() - height - 6.0*particleRadius;
        //only keep wall particles
        for (BaseParticle* p: particleHandler)
        {
            Mdouble r = sqrt(mathsFunc::square(p->getPosition().X) + mathsFunc::square(p->getPosition().Y));
            Mdouble z = p->getPosition().Z;
            //~ if ((r < innerRadius+wallThickness_ && r > innerRadius && z > dHeight)
                //~ || (r > outerRadius-wallThickness_ && r < outerRadius && z > dHeight)
                //~ || (r > innerRadius && r < outerRadius && z > dHeight && z < dHeight+wallThickness_))
            if ((r > outerRadius-wallThickness_ && r < outerRadius && z > dHeight)
                || (r > innerRadius && r < outerRadius && z > dHeight && z < dHeight+wallThickness_))
            {
                p->setPosition(p->getPosition() - Vec3D(0, 0, dHeight));
                p->fixParticle();
                //std::cout << ".";
            }
        }
        for (int i = particleHandler.getNumberOfObjects()-1; i>=0; i--)
        {
            if (!particleHandler.getObject(i)->isFixed())
                particleHandler.removeObject(i);
        }
        std::cout << "kept " << particleHandler.getNumberOfObjects() << " fixed particles" << std::endl;

        //remove walls, reset height
        ///\todo interactions should be cleared when particle gets erased
        wallHandler.clear();
        setXMax(outerRadius);
        setZMax(height);

        //set walls
        AxisymmetricIntersectionOfWalls w;
        Vec3D AxisDirection(0.,0.,1.);
        Vec3D PointOnAxis(0.,0.,0.);
        Vec3D PrismAxis(0.,-1.,0.);
        w.setPosition(PointOnAxis);
        w.setOrientation(AxisDirection);
        w.setSpecies(species);

        //add points in anti-clockwise direction around the prism axis
        std::vector<Vec3D> Points(2);
        Points[0] = Vec3D(innerRadius+wallThickness_,0.0,getZMin());
        Points[1] = Vec3D(innerRadius+wallThickness_,0.0,getZMax());
        w.createOpenPrism(Points,PrismAxis);
        w.setSpecies(speciesInnerWall);
        wallHandler.copyAndAddObject(w);

        Points.resize(2);
        Points[0] = Vec3D(outerRadius,0.0,getZMax());
        Points[1] = Vec3D(outerRadius,0.0,getZMin());
        w.createOpenPrism(Points,PrismAxis);
        w.setSpecies(species);
        wallHandler.copyAndAddObject(w);

        InfiniteWall iw;
        iw.set(Vec3D(0.0,0.0,-1.0), Vec3D(0.0,0.0,0.0));
        wallHandler.copyAndAddObject(iw);

        interactionHandler.clear();

        setName("ShearCell3DInitialConditionsWall");
        writeDataFile();
        writeRestartFile();
    }

private:
    /////////////////////////////////////////////////////////////////////////////////////////////////////
    ///This is were the walls are implemented
    /////////////////////////////////////////////////////////////////////////////////////////////////////
    void setupInitialConditions()
    {
        //hGridRebuild();
        for (unsigned int i = particleHandler.getNumberOfObjects(); i >= 1; i--)
            for (BaseBoundary *it : boundaryHandler)
                it->createPeriodicParticle(particleHandler.getObject(i-1), particleHandler);

        std::cout << "setupInitialConditions()" << std::endl;
        SphericalParticle p;
        p.setSpecies(speciesHandler.getObject(0));
        p.setVelocity(Vec3D(0.0, 0.0, 0.0));
        p.setOrientation(Vec3D(0.0, 0.0, 0.0));
        p.setAngularVelocity(Vec3D(0.0, 0.0, 0.0));
        //this sets the size of the first particle (set it to the largest particle size)
        p.setRadius(particleRadius_*sqrt(2.0));
        // Monodisperse case
        //p.setRadius(particleRadius_);
        Mdouble volumeOfParticlesToInsert = 0.5 * openingAngle_ * height_ * (mathsFunc::square(outerRadius_-1.6*particleRadius_) - mathsFunc::square(innerRadius_+1.6*particleRadius_));
        Mdouble h = getZMin();
        Mdouble Ninit = particleHandler.getNumberOfObjects();
        while (volumeOfParticlesToInsert > 0.0)
        {
            Mdouble r = random.getRandomNumber(innerRadius_+1.6*particleRadius_ , outerRadius_-1.6*particleRadius_);
            Mdouble a = random.getRandomNumber(0, openingAngle_);
            Mdouble z = random.getRandomNumber(getZMin(), h);
            p.setPosition(Vec3D(r * std::cos(a), r * std::sin(a), z));

            if (checkParticleForInteraction(p))
            {
                BaseParticle* q = particleHandler.copyAndAddObject(p);
                //duplicate particle
                for (BaseBoundary *it : boundaryHandler)
                    it->createPeriodicParticle(q, particleHandler);
                //this sets the radius of all other particles
                 p.setRadius(random.getRandomNumber(1.0, 2.0) * particleRadius_ / std::sqrt(2.0));
                
		// for the bidisperse mixture uncomment below
                
                //if (random.getRandomNumber(0.0, 1.0)>0.5)
		//			p.setRadius(particleRadius_*std::sqrt(2.0));
                //else
		//			p.setRadius(particleRadius_/std::sqrt(2.0));
                
                // Monodisperse case
                //p.setRadius(particleRadius_);
                
                volumeOfParticlesToInsert -= 8.0 * mathsFunc::cubic(p.getRadius())*(constants::pi/6.0)/0.66; //here we assume a volume fraction of 0.66 but for friction I changed to 0.6
            }
            else
                h += 0.0002 * particleRadius_;
        }
        for (unsigned int i = particleHandler.getNumberOfObjects(); i >= 1; i--)
            if (particleHandler.getObject(i - 1)->getPeriodicFromParticle() != nullptr)
            {
                while (particleHandler.getObject(i - 1)->getInteractions().size() > 0)
                {
                    interactionHandler.removeObjectKeepingPeriodics(particleHandler.getObject(i - 1)->getInteractions().front()->getIndex());
                }
                particleHandler.removeObject(i - 1);
            }
		
        std::cout << "Inserted " << particleHandler.getNumberOfObjects()-Ninit << " particles" << std::endl;

        setHGridMaxLevels(2);
        std::cout << "simulating " << particleHandler.getNumberOfObjects() << " particles; h=" << h << std::endl;
    }

protected:
    Mdouble wallThickness_;
    Mdouble particleRadius_;
    Mdouble outerRadius_;
    Mdouble innerRadius_;
    Mdouble height_;
    Mdouble openingAngle_;
public:
    LinearPlasticViscoelasticSlidingFrictionSpecies* species;
    LinearPlasticViscoelasticSlidingFrictionSpecies* speciesInnerWall;
};

class ShearCell3DInitialConditions : public ShearCell3DWalls {
public:
    ShearCell3DInitialConditions (Mdouble particleRadius, Mdouble outerRadius, Mdouble innerRadius, Mdouble height, Mdouble openingAngle, Mdouble wallThickness)
        : ShearCell3DWalls (particleRadius, outerRadius, innerRadius, height, openingAngle, wallThickness)
    {
        setName("ShearCell3DInitialConditionsWall");
        readRestartFile();
        setRestarted(false);
        setName("ShearCell3DAbhiInitialConditionsSmoothInnerWall");
        species = dynamic_cast<LinearPlasticViscoelasticSlidingFrictionSpecies*>(speciesHandler.getObject(0));
    }

private:

    void actionsAfterTimeStep()
    {
        // fix particles that fall through the wall
        if (wallHandler.getNumberOfObjects()==0)
            for (BaseParticle* p: particleHandler)
            {
                if (p->isFixed())
                    continue;
                Mdouble rr = mathsFunc::square(p->getPosition().X) + mathsFunc::square(p->getPosition().Y);
                Mdouble z = p->getPosition().Z;
                if ((rr < mathsFunc::square(innerRadius_-particleRadius_-0.5* wallThickness_))
                    || (rr > mathsFunc::square(outerRadius_+particleRadius_+0.5* wallThickness_))
                    || (z < -particleRadius_-0.5* wallThickness_) )
                {
                    p->fixParticle();
                    std::cout << "fixed particle that fell out" << std::endl;
                }
            }
    }
};

int main(int argc, char *argv[]) {
    Mdouble particleRadius = 1.1e-3*1.0;
    Mdouble outerRadius = 110e-3;
    Mdouble innerRadius = std::max(14.7e-3,13*particleRadius);
    Mdouble splitRadius = 85e-3;
    Mdouble height = 1.0*splitRadius;
    //Mdouble height = 38e-3; //it assumes a volume fraction of 0.66, while changing please make sure volumefraction*height is 0.66*0.036.
    //Mdouble height = 36e-3; //it assumes a volume fraction of 0.66, while changing please make sure volumefraction*height is 0.66*0.036.
    Mdouble openingAngle = 30.0*constants::pi/180.0;
    Mdouble wallThickness = 2.2*particleRadius;

    //std::cout << "creating the mold for the bottom and side walls..." << std::endl;
    //ShearCell3DWalls SC(particleRadius, outerRadius+particleRadius, innerRadius-particleRadius, height+6.0*particleRadius+particleRadius, openingAngle, wallThickness);
    //SC.setXBallsAdditionalArguments("-v0 -solidf -3dturn 1 -w 1180 -s 24");
    //SC.setSaveCount(10000);
    //SC.solve();
    //std::cout << "creating the bottom and side walls from the mold and filling it ..." << std::endl;
    //SC.createWalls(particleRadius, outerRadius, innerRadius, height, openingAngle);

    ShearCell3DInitialConditions IC(particleRadius, outerRadius, innerRadius, height, openingAngle, wallThickness);
    IC.solve();

    std::cout << "reset dissipation and time step ..." << std::endl;
    IC.species->setDissipation(2e-3*10);
    IC.species->setSlidingDissipation(2e-3*10);
    IC.setTimeStep(IC.species->computeTimeStep(IC.species->getMassFromRadius(particleRadius)));
    IC.interactionHandler.clear();
    IC.writeRestartFile();
    return 0;
}
