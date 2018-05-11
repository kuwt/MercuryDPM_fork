#include "Mercury3D.h"
#include "TriolietScrew.h"
#include "TriolietBaseScrew.h"
#include "Species/LinearViscoelasticSlidingFrictionSpecies.h"
#include "Walls/RestrictedWall.h"
#include "Walls/InfiniteWall.h"
#include "Walls/IntersectionOfWalls.h"
#include "Walls/AxisymmetricIntersectionOfWalls.h"

#ifndef TRIOLIETDIETFEEDER_H
#define TRIOLIETDIETFEEDER_H

class TriolietDietFeeder : public Mercury3D
{
private:
    /*!
     * \brief The mean radius of the particles in the feeder
     */
    Mdouble particleRadius;
    /*!
     * \brief Pointer to the left screw
     */
    TriolietScrew* leftScrew = nullptr;
    /*!
     * \brief Pointer to the right screw
     */
    TriolietScrew* rightScrew = nullptr;
    /*!
     * \brief The rotation speed of the screw
     */
    Mdouble rotationSpeed;
    /*!
     * \brief The time where the screw begins to spin
     */
    Mdouble timeMin;
    /*!
     * \brief The index number of all mounted blades (the blades mounts are numbered 0-11, with the i-th blade mounts at relative height q=(1+2i)/25)
     */
    std::vector<unsigned> bladeMounts_;

public:

    /*!
     * \brief Constructor, turns off fstat output by default
     */
    TriolietDietFeeder(Mdouble particleRadius, Mdouble rotationSpeed, Mdouble timeMin) : particleRadius(particleRadius), rotationSpeed(rotationSpeed), timeMin(timeMin)
    {
        fStatFile.setFileType(FileType::NO_FILE);
    }

    /*!
     * \brief Set function for the bladeMounts index.
     * \detailed The blades are actually inserted into the screw by the function setBlades.
     */
    void setBlades(const std::vector<unsigned> bladeMounts) {
        logger(INFO,"Adding % blades to each screw", bladeMounts.size());
        bladeMounts_ = bladeMounts;
    }

    /*!
     * \brief Inserts blades into the screw.
     * \detailed This function is called in setupInitialConditions, right after the screws are defined.
     */
    void setBlades()
    {
        // radial width of the screw
        Mdouble bladeWidth = 150e-3;
        // maximum value is 1/25, which is the relative height between two blade mounts.
        Mdouble bladeLength = 1.0/25.0;
        // now transform the index of each mount into the actual relative height on the screw.
        std::vector<Mdouble> bladeMounts;
        for (auto n : bladeMounts_) {
            bladeMounts.push_back(static_cast<Mdouble>(2 * n + 1) / 25.0);
        }
        leftScrew->setBlades(bladeWidth, bladeLength, bladeMounts);
        rightScrew->setBlades(bladeWidth, bladeLength, bladeMounts);
    }

    ///sets four walls, leftScrew, rightScrew, leftBaseScrew, rightBaseScrew
    void setScrewWalls(const ParticleSpecies* s, Mdouble screwCenter,
        Mdouble screwBaseHeight, Mdouble screwBaseRadius, Mdouble screwTopHeight, 
        Mdouble windingLength, Mdouble minR, Mdouble lowerR, Mdouble diffR, Mdouble thickness)
    {
        //Screws are dx=2200 apart, outerTopHeight height 1236, bottom height 110 (plus cap)
        //winding width 411
        Vec3D leftStart = Vec3D(-screwCenter, 0, screwBaseHeight);
        Vec3D rightStart = Vec3D(screwCenter, 0, screwBaseHeight);
        Mdouble length = screwTopHeight-screwBaseHeight;
        Mdouble numTurns = length/windingLength;
        //create a screw that increases its length clockwise, just like (sin(t),cos(t), t).
        // screw starts at height 0.11 until 1.197, distance 1.1 from the center in negative x-direction, making one turn each 0.411.
        // rotation speed 1.0, thickness as small as possible, radius function is piecewise linear, decreasing from 0.6 to 0.4.
        // screw starts at height 0.11 until 1.197, making one turn each 0.411.
        leftScrew = wallHandler.copyAndAddObject
           (TriolietScrew(leftStart, length, minR, lowerR, diffR, numTurns, -rotationSpeed/2.0/constants::pi, thickness, s));       
        // same as left screw, but increases its length anticlockwise, just like (-sin(t),cos(t), t), 
        // and rotating the other way.
        // distance 0.6 from the center in positive x-direction
        rightScrew = wallHandler.copyAndAddObject
            (TriolietScrew(rightStart, length, minR, lowerR, diffR, numTurns, -rotationSpeed/2.0/constants::pi, thickness, s));
        //todo: change orientation of one screw

        Mdouble sinA2Max = 0.25;
        auto leftBaseScrew = wallHandler.copyAndAddObject
            (TriolietBaseScrew(Vec3D(-screwCenter,0,0),Vec3D(0,0,1),
                {{Vec3D(0,0,-1),Vec3D(0,0,screwBaseHeight)},
                {Vec3D(-1,0,0),Vec3D(screwBaseRadius,0,0)}},s,sinA2Max,timeMin));

        auto rightBaseScrew = wallHandler.copyAndAddObject
            (TriolietBaseScrew(Vec3D(screwCenter,0,0),Vec3D(0,0,1),
                {{Vec3D(0,0,-1),Vec3D(0,0,screwBaseHeight)},
                {Vec3D(-1,0,0),Vec3D(screwBaseRadius,0,0)}},s,sinA2Max,timeMin));
    }

    ///sets four walls, leftScrewCore, rightScrewCore, leftScrewBottom, rightScrewBottom
    void setScrewCore(const ParticleSpecies* s, Mdouble screwCenter, Mdouble screwBaseHeight,
        Mdouble coreTopHeight, Mdouble coreTopRadius, 
        Mdouble coreBottomHeight, Mdouble coreBottomRadius)
    {
        //the inner, thinner core cylinder in the screw
        auto leftScrewCore = wallHandler.copyAndAddObject
            (AxisymmetricIntersectionOfWalls(Vec3D(-screwCenter,0,0),Vec3D(0,0,1),
                {{Vec3D(-.3,0,-1),Vec3D(coreTopRadius,0,coreTopHeight)},//slightly slant top, so particles roll off
                {Vec3D(-1,0,0),Vec3D(coreTopRadius,0,0)}},s));

        auto rightScrewCore = wallHandler.copyAndAddObject
            (AxisymmetricIntersectionOfWalls(Vec3D(screwCenter,0,0),Vec3D(0,0,1),
                {{Vec3D(-.3,0,-1),Vec3D(0,0,coreTopHeight)},
                {Vec3D(-1,0,0),Vec3D(coreTopRadius,0,0)}},s));

        //the bottom, thicker core cylinder in the screw
        //set angle on top
        double Angle=constants::pi/3;
        double Xtan=sin(Angle)/cos(Angle);

        auto leftScrewBottom = wallHandler.copyAndAddObject
            (AxisymmetricIntersectionOfWalls(Vec3D(-screwCenter,0,0),Vec3D(0,0,1),
                {{Vec3D(-Xtan,0,-1),Vec3D(0,0,coreBottomHeight+coreBottomRadius*Xtan)},
                {Vec3D(-1,0,0),Vec3D(coreBottomRadius,0,0)}},s));

        auto rightScrewBottom = wallHandler.copyAndAddObject
            (AxisymmetricIntersectionOfWalls(Vec3D(screwCenter,0,0),Vec3D(0,0,1),
                {{Vec3D(-Xtan,0,-1),Vec3D(0,0,coreBottomHeight+coreBottomRadius*Xtan)},
                {Vec3D(-1,0,0),Vec3D(coreBottomRadius,0,0)}},s));        
    }

    ///sets other walls that define the outer boundary
    void setOuterWalls(const ParticleSpecies* s, Mdouble outerBaseRadius, Mdouble screwCenter, Mdouble outerTopCenter,
        Mdouble outerTopRadius, Mdouble outerTopHeight)
    {
        Vec3D orientation = (Vec3D(outerTopCenter,0,outerTopHeight)-Vec3D(screwCenter,0,0));
        orientation.normalize();
        Vec3D position = Vec3D(screwCenter,0,0);
        Vec3D normalWall = Vec3D(outerTopHeight,0,outerBaseRadius-outerTopRadius);
        normalWall.normalize();
        Vec3D positionWall = Vec3D(outerBaseRadius,0,0);
        auto fullRightSide = AxisymmetricIntersectionOfWalls(position, orientation, {{normalWall,positionWall}},s);
        auto rightRestriction = InfiniteWall(Vec3D(orientation.Z,0,-orientation.X),position,s);      
        auto rightSide = wallHandler.copyAndAddObject
            (RestrictedWall(&fullRightSide, &rightRestriction));

        orientation.X = -orientation.X;
        position = -position;
        auto fullLeftSide = AxisymmetricIntersectionOfWalls(position, orientation, {{normalWall,positionWall}},s);
        auto leftRestriction = InfiniteWall(Vec3D(-orientation.Z,0,orientation.X),position,s);      
        auto leftSide = wallHandler.copyAndAddObject
            (RestrictedWall(&fullLeftSide, &leftRestriction));
        
        position = Vec3D(0,outerBaseRadius,0);
        normalWall = Vec3D(0,outerTopHeight,outerBaseRadius-outerTopRadius);
        auto backPlate = wallHandler.copyAndAddObject
            (InfiniteWall(normalWall,position,s));
        position = -position;
        normalWall.Y = -normalWall.Y;
        auto frontPlate = wallHandler.copyAndAddObject
            (InfiniteWall(normalWall,position,s));
        
        auto bottomPlate = wallHandler.copyAndAddObject
            (InfiniteWall(Vec3D(0,0,-1),Vec3D(0,0,getZMin()),s));  

        //narrowing in the center
        std::vector<Vec3D> lPoints;
        std::vector<Vec3D> rPoints;
        //from 0 to 80 degrees, 7 points on circle, expect 4 points
        Mdouble n = 7;
        for (Mdouble i=0; i<n; i++) {
            Mdouble a = 80.0*i/n;
            Mdouble s=sin(constants::pi/180.0*a);
            Mdouble c=cos(constants::pi/180.0*a);
            lPoints.push_back(Vec3D(screwCenter-s*outerBaseRadius,c*outerBaseRadius,0));
            if (rPoints.size()<=n-4) {
                rPoints.push_back(Vec3D(-screwCenter+s*outerBaseRadius,c*outerBaseRadius,0));
            } 
        }
        //4 points on one side on line
        Vec3D right = rPoints.back();
        rPoints.push_back((lPoints.back()+3.0*right)/4.0);
        rPoints.push_back((lPoints.back()+right)/2.0);
        rPoints.push_back((3.0*lPoints.back()+right)/4.0);
        rPoints.push_back(lPoints.back());
        //define points on outerTopHeight
        Vec3D lTop = Vec3D(-0.3,outerTopRadius,0.8*outerTopHeight);
        Vec3D rTop = Vec3D( 0.3,outerTopRadius,0.8*outerTopHeight);

        //add narrowing in the center for y>0
        for (unsigned i=1; i<lPoints.size(); i++)
        {
            auto side = wallHandler.copyAndAddObject(IntersectionOfWalls({},s));  
            side->addObject(Vec3D::cross(lTop-lPoints[i-1],lTop-lPoints[i]),lTop);
            side->addObject(Vec3D::cross(rTop-lPoints[i],rTop-rPoints[i]),rTop);
            side->addObject(Vec3D::cross(rTop-rPoints[i],rTop-rPoints[i-1]),rTop);
        }

        //add narrowing in the center for y<0
        for (auto& p : lPoints) p = -p;
        for (auto& p : rPoints) p = -p;
        lTop.Y = -lTop.Y;
        rTop.Y = -rTop.Y;
        lTop.X = -lTop.X;
        rTop.X = -rTop.X;
        for (unsigned i=1; i<lPoints.size(); i++)
        {
            auto side = wallHandler.copyAndAddObject(IntersectionOfWalls({},s));  
            side->addObject(Vec3D::cross(lTop-lPoints[i-1],lTop-lPoints[i]),lTop);
            side->addObject(Vec3D::cross(rTop-lPoints[i],rTop-rPoints[i]),rTop);
            side->addObject(Vec3D::cross(rTop-rPoints[i],rTop-rPoints[i-1]),rTop);
        }
    }

    void setupInitialConditions() override
    {        
        ///Here wall properties are set
        const ParticleSpecies* s = speciesHandler.getObject(0);
        
        // First, we define the screw
        // bottom radius 0.974, outerTopHeight 1210, center-bottom x=1.21, center-outerTopHeight x=1.3385
        Mdouble screwCenter = 1.1; //changed from 1.1
        Mdouble screwBaseHeight = 0.11; //same
        Mdouble screwTopHeight = 1.155; //changed from 1.155 //max. Height of the screw/somehow, screw is still felt on top of the core
        Mdouble windingLength = 0.411; //changed from 0.411
        Mdouble minR = 0.4;//change from 0.4 //radius in upper part
        Mdouble lowerR = 0.6; //change from 0.6 ??//radius in upper part
        Mdouble diffR = -0.3; //changed from -0.3 //radius in upper part
        Mdouble thickness = 0.5*particleRadius; //needs change?

        //bottom screw core radius 270, outerTopHeight height 336
        //outerTopHeight screw core radius ~170, outerTopHeight height 1236
        Mdouble coreTopHeight = 1.236; //changed from 1236
        Mdouble coreTopRadius = 0.17; //same
        Mdouble coreBottomHeight = .356; //changed from 0.356
        Mdouble coreBottomRadius = 0.27; //same
        Mdouble screwBaseRadius = 0.9; //changed from 0.9

        //the cylindrical container:
        Mdouble outerBaseRadius  =0.974; //changed from 0.974
        Mdouble outerTopCenter = screwCenter;// should be 1.3385, but for that I need to correct the implementation of orientation
        Mdouble outerTopRadius  =1.21; //changed from 1.21
        Mdouble outerTopHeight = 2.252; // changed from 2.252

        setXMax(outerTopCenter+outerTopRadius);
        setYMax(outerTopRadius);
        setZMax(outerTopHeight);
        setXMin(-getXMax());
        setYMin(-getYMax());
        setZMin(0.0);

        setScrewWalls(s, screwCenter, screwBaseHeight, screwBaseRadius, screwTopHeight, windingLength, minR, lowerR, diffR, thickness);
        setBlades();
        setScrewCore(s, screwCenter, screwBaseHeight, coreTopHeight, coreTopRadius, coreBottomHeight, coreBottomRadius);
        setOuterWalls(s, outerBaseRadius, screwCenter, outerTopCenter, outerTopRadius, outerTopHeight);

        //introduceSingleParticle(Vec3D(screwCenter,0.5*coreTopRadius,coreTopHeight));
        //introduceParticlesAtWall();
        introduceParticlesInDomain(1.3);
    }

    void introduceSingleParticle(Vec3D p) {
                /* Simple run settings
         * Nx*Ny*Nz particles are created evenly spaced between [xmin,xmax]*[ymin,ymax]*[zmin,zmax] and checked for contact with the screw
         */
        const ParticleSpecies* s = speciesHandler.getObject(0);
        particleHandler.clear();
        BaseParticle p0;
        p0.setSpecies(s);
        p0.setVelocity(Vec3D(0.0, 0.0, 0.0));
        p0.setRadius(particleRadius);
        p0.setPosition(p+Vec3D(0,0,particleRadius));
        particleHandler.copyAndAddObject(p0);
        std::cout << "Inserted single particle " << std::endl;
    }

    void introduceParticlesAtWall() {
        /* Simple run settings
         * Nx*Ny*Nz particles are created evenly spaced between [xmin,xmax]*[ymin,ymax]*[zmin,zmax] and checked for contact with the screw
         */
        const ParticleSpecies* s = speciesHandler.getObject(0);
        particleHandler.clear();
        BaseParticle p0;
        p0.setSpecies(s);
        p0.setVelocity(Vec3D(0.0, 0.0, 0.0));
        p0.setRadius(2.0*particleRadius);
        BaseParticle p1;
        p1.setSpecies(s);
        p1.setVelocity(Vec3D(0.0, 0.0, 0.0));
        p1.setRadius(particleRadius);

        //domain size d
        Vec3D d = Vec3D(getXMax() - getXMin(), getYMax() - getYMin(), getZMax() - getZMin());
        //number of particles that fit in domain
        Mdouble Nx = floor(d.X / (2.0 * particleRadius));
        Mdouble Ny = floor(d.Y / (2.0 * particleRadius));
        Mdouble Nz = floor(d.Z / (2.0 * particleRadius));
        
        Mdouble distance;
        Vec3D normal;
        Vec3D p;
        Mdouble minDistance;
        unsigned counter = 0;
        for (p.X = getXMin() + particleRadius; p.X < getXMax(); p.X += 2.0 * particleRadius)
            for (p.Y = getYMin() + particleRadius; p.Y < getYMax(); p.Y += 2.0 * particleRadius)
                for (p.Z = getZMin() + particleRadius; p.Z < getZMax(); p.Z += 2.0 * particleRadius)
                {
                    minDistance = p0.getRadius();
                    p0.setPosition(p);
                    for (auto w : wallHandler) {
                        //if touching the wall
                        if (w->getDistanceAndNormal(p0, distance, normal) && distance<minDistance)
                        {
                            minDistance=distance;
                            if (distance<0) break; 
                            p1.setPosition(p0.getPosition()+(distance-particleRadius)*normal);
                        }
                    }
                    if (minDistance<p0.getRadius() && minDistance>0)
                    {
                        particleHandler.copyAndAddObject(p1);
                        counter++;
                    }
                }
        std::cout << "Inserted particles: " << counter << std::endl;
    }

    void introduceParticlesInDomain(Mdouble polydispersity=1) {
        /* Simple run settings
         * Nx*Ny*Nz particles are created evenly spaced between [xmin,xmax]*[ymin,ymax]*[zmin,zmax] and checked for contact with the screw
         */
        const ParticleSpecies* s = speciesHandler.getObject(0);
        particleHandler.clear();
        BaseParticle p0;
        p0.setSpecies(s);
        p0.setVelocity(Vec3D(0.0, 0.0, 0.0));
        p0.setRadius(particleRadius);
        
        //domain size d
        Vec3D d = Vec3D(getXMax() - getXMin(), getYMax() - getYMin(), getZMax() - getZMin());
        //number of particles that fit in domain
        Mdouble Nx = floor(d.X / (2.0 * particleRadius));
        Mdouble Ny = floor(d.Y / (2.0 * particleRadius));
        Mdouble Nz = floor(d.Z / (2.0 * particleRadius));

        //CHANGED BY BERT need ~ 2.3 times cubic packing fraction to get HCP packing, assume particles are more HCP packed than cubic
        Mdouble distance;
        Vec3D normal;
        Vec3D p;
        Mdouble minDistance;
        unsigned counter = 0;
        for (p.X = getXMin() + particleRadius; p.X < getXMax(); p.X += 2.0 * particleRadius)
            for (p.Y = getYMin() + particleRadius; p.Y < getYMax(); p.Y += 2.0 * particleRadius)
                for (p.Z = getZMin() + particleRadius; p.Z < 2.3*getZMax(); p.Z += 2.0 * particleRadius) //Changed Cubic to HCP here
                {
                    bool touch = false;
                    p0.setPosition(p);
                    for (auto w : wallHandler) {
                        //if touching the wall
                        if (w->getDistanceAndNormal(p0, distance, normal))
                        {
                            touch = true;
                            break;
                        }
                    }
                    if (!touch) {
                        particleHandler.copyAndAddObject(p0);
                        p0.setRadius(random.getRandomNumber(particleRadius/polydispersity,particleRadius));
                        counter++;
                    }
                }
        std::cout << "Inserted particles: " << counter << std::endl;
    }

    void outputXBallsData (std::ostream& os) const override {
        DPMBase::outputXBallsData(os);
        wallHandler.writeVTK();
    }

    void printTime () const override {
        std::cout << "t=" << std::setprecision(3) << std::left << std::setw(6) << getTime()
            << ", tmax=" << std::setprecision(3) << std::left << std::setw(6) << getTimeMax()
            << ", EneRatio=" << std::setprecision(3) << std::left << std::setw(6) << getKineticEnergy()/getElasticEnergy()
            << std::endl;
        std::cout.flush();
    }    

    std::vector<Quaternion> orientations;

    void actionsBeforeTimeLoop () override {
        for (const auto w: wallHandler) {
            if (dynamic_cast<RestrictedWall*>(w) != nullptr) break;
            orientations.push_back(w->getOrientation());
        }
    }

    void actionsBeforeTimeStep () override {
        auto w = wallHandler.begin();
        for (auto o = orientations.begin(); o!=orientations.end(); w++, o++) {
            (*w)->setOrientation(*o);
        }
    }

    void actionsAfterTimeStep () override {
        if (getTime()>timeMin)
        {
            if (leftScrew) leftScrew->move_time(getTimeStep());
            if (rightScrew) rightScrew->move_time(getTimeStep());
            if (getTime()-getTimeStep()<=timeMin) {
                std::cout << "start rotation" << std::endl;
                auto w = wallHandler.begin();
                for (auto o = orientations.begin(); o!=orientations.end(); w++, o++) {
                    (*w)->setAngularVelocity(Vec3D(0,0,rotationSpeed));
                }
            }
        }
    }

};

#endif
