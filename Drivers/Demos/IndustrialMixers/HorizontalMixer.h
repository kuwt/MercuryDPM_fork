#include "Mercury3D.h"
#include "Species/LinearViscoelasticSlidingFrictionSpecies.h"
#include "Walls/RestrictedWall.h"
#include "Walls/InfiniteWall.h"
#include "Walls/IntersectionOfWalls.h"
#include "Walls/AxisymmetricIntersectionOfWalls.h"
#include "Walls/HorizontalScrew.h"
#include "Walls/HorizontalBaseScrew.h"

#ifndef HorizontalMixer_h
#define HorizontalMixer_h

class HorizontalMixer : public Mercury3D
{
private:
    /*!
     * \brief The mean radius of the particles in the feeder
     */
    Mdouble particleRadius;
    /*!
     * \brief Pointer to the right screw
     */
    HorizontalScrew* screw = nullptr;
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
    /**
     * Fill Height
     */
    Mdouble fillHeight_;

public:
    /**
     * determines if outerWalls are included (useful for printing)
     */
    bool haveOuterWalls = true;


public:

    /*!
     * \brief Constructor, turns off fstat output by default
     */
    HorizontalMixer(Mdouble particleRadius, Mdouble rotationSpeed, Mdouble timeMin, Mdouble fillHeight)
            : particleRadius(particleRadius), rotationSpeed(rotationSpeed), timeMin(timeMin), fillHeight_(fillHeight)
    {
        fStatFile.setFileType(FileType::NO_FILE);
        setParticlesWriteVTK(true);
        setWallsWriteVTK(FileType::MULTIPLE_FILES);
    }

    ///sets four walls, leftScrew, rightScrew, leftBaseScrew, rightBaseScrew
    void setScrewWalls(const ParticleSpecies* s, Mdouble screwCenter,
        Mdouble screwBaseHeight, Mdouble screwBaseRadius, Mdouble screwTopHeight, 
        Mdouble windingLength, Mdouble minR, Mdouble lowerR, Mdouble diffR, Mdouble thickness)
    {
        //Screws are dx=2200 apart, outerTopHeight height 1236, bottom height 110 (plus cap)
        //winding width 411
        Vec3D start = Vec3D(screwCenter, 0, screwBaseHeight);
        Mdouble length = screwTopHeight-screwBaseHeight;
        Mdouble numTurns = length/windingLength;
        //create a screw that increases its length clockwise, just like (sin(t),cos(t), t).
        // screw starts at height 0.11 until 1.197, distance 1.1 from the center in negative x-direction, making one turn each 0.411.
        // rotation speed 1.0, thickness as small as possible, radius function is piecewise linear, decreasing from 0.6 to 0.4.
        // screw starts at height 0.11 until 1.197, making one turn each 0.411.
        screw = wallHandler.copyAndAddObject
           (HorizontalScrew(start, length, minR, lowerR, diffR, numTurns, rotationSpeed/2.0/constants::pi, thickness, s));
        //todo: change orientation of one screw
        if (screw) screw->setAngularVelocity(Vec3D(0,0,rotationSpeed));
        if (screw) screw->move_time(40);

        Mdouble sinA2Max = 0.25;
        auto baseScrew = wallHandler.copyAndAddObject
            (HorizontalBaseScrew(Vec3D(screwCenter,0,0),Vec3D(0,0,1),
                {{Vec3D(0,0,-1),Vec3D(0,0,screwBaseHeight)},
                {Vec3D(-1,0,0),Vec3D(screwBaseRadius,0,0)}},s,sinA2Max,timeMin));
        baseScrew->setAngularVelocity(Vec3D(0,0,rotationSpeed));
        //baseScrew->rotate(Vec3D(0,0,100));

    }

    ///sets four walls, leftScrewCore, rightScrewCore, leftScrewBottom, rightScrewBottom
    virtual void setScrewCore(const ParticleSpecies* s, Mdouble screwCenter, Mdouble screwBaseHeight,
        Mdouble coreTopHeight, Mdouble coreTopRadius, 
        Mdouble coreBottomHeight, Mdouble coreBottomRadius)
    {
        //the inner, thinner core cylinder in the screw
        auto screwCore = wallHandler.copyAndAddObject
            (AxisymmetricIntersectionOfWalls(Vec3D(screwCenter,0,0),Vec3D(0,0,1),
                {{Vec3D(-.3,0,-1),Vec3D(coreTopRadius,0,coreTopHeight)},//slightly slant top, so particles roll off
                {Vec3D(-1,0,0),Vec3D(coreTopRadius,0,0)}},s));
        screwCore->setAngularVelocity(Vec3D(0,0,rotationSpeed));

        //the bottom, thicker core cylinder in the screw
        auto screwBottom = wallHandler.copyAndAddObject
            (AxisymmetricIntersectionOfWalls(Vec3D(screwCenter,0,0),Vec3D(0,0,1),
                {{Vec3D(0,0,-1),Vec3D(0,0,coreBottomHeight)},
                {Vec3D(-1,0,0),Vec3D(coreBottomRadius,0,0)}},s));
        screwBottom->setAngularVelocity(Vec3D(0,0,rotationSpeed));
    }

    ///sets other walls that define the outer boundary
    virtual void setOuterWalls(const ParticleSpecies* s, Mdouble outerBaseRadius, Mdouble screwCenter, Mdouble outerTopCenter,
        Mdouble outerTopRadius, Mdouble outerTopHeight)
    {
        //cylindrical wall
        Vec3D normal = (Vec3D(outerTopCenter,0,outerTopHeight)-Vec3D(screwCenter,0,0));
        normal.normalise();
        Vec3D position = Vec3D(screwCenter,0,0);
        Vec3D normalWall = Vec3D(outerTopHeight,0,outerBaseRadius-outerTopRadius);
        normalWall.normalise();
        Vec3D positionWall = Vec3D(outerBaseRadius,0,0);
        AxisymmetricIntersectionOfWalls outerWall(position, normal, {{normalWall,positionWall}},s);
        auto rightSide = wallHandler.copyAndAddObject(outerWall);

        //bottom plate
        auto bottomPlate = wallHandler.copyAndAddObject
            (InfiniteWall(Vec3D(0,0,-1),Vec3D(0,0,getZMin()),s));
    }

    void setupInitialConditions() override
    {        
        ///Here wall properties are set
        const ParticleSpecies* s = speciesHandler.getObject(0);
        
        // First, we define the screw
        // bottom radius 0.974, outerTopHeight 1210, center-bottom x=1.21, center-outerTopHeight x=1.3385
        Mdouble screwCenter = 1.21; //changed from 1.1
        Mdouble screwBaseHeight = 0.11; //same
        Mdouble screwTopHeight = 1.5; //changed from 1.155 //max. Height of the screw/somehow, screw is still felt on top of the core
        Mdouble windingLength = 0.5; //changed from 0.411
        Mdouble minR = 0.54;//change from 0.4 //radius in upper part
        Mdouble lowerR = 1.195; //change from 0.6 ??//radius in upper part
        Mdouble diffR = -0.6; //changed from -0.3 //radius in upper part
        Mdouble thickness = 0.5*particleRadius; //needs change?

        //bottom screw core radius 270, outerTopHeight height 336
        //outerTopHeight screw core radius ~170, outerTopHeight height 1236
        Mdouble coreTopHeight = 1.5; //changed from 1236
        Mdouble coreTopRadius = 0.17; //same
        Mdouble coreBottomHeight = .491; //changed from 0.356
        Mdouble coreBottomRadius = 0.27; //same
        Mdouble screwBaseRadius = 1.0; //changed from 0.9

        //the cylindrical container:
        Mdouble outerBaseRadius  =1.208; //changed from 0.974
        Mdouble outerTopCenter = screwCenter;// should be 1.3385, but for that I need to correct the implementation of orientation
        Mdouble outerTopRadius  =1.329; //changed from 1.21
        Mdouble outerTopHeight = 1.995; // changed from 2.252

        setXMax(screwCenter+outerTopRadius);
        setYMax(outerTopRadius);
        setZMax(outerTopHeight);
        setXMin(screwCenter-outerTopRadius);
        setYMin(-getYMax());
        setZMin(0.0);

        setScrewWalls(s, screwCenter, screwBaseHeight, screwBaseRadius, screwTopHeight, windingLength, minR, lowerR, diffR, thickness);
        setScrewCore(s, screwCenter, screwBaseHeight, coreTopHeight, coreTopRadius, coreBottomHeight, coreBottomRadius);
        if (haveOuterWalls) setOuterWalls(s, outerBaseRadius, screwCenter, outerTopCenter, outerTopRadius, outerTopHeight);

        //introduceSingleParticle(Vec3D(screwCenter,0.5*coreTopRadius,coreTopHeight));
        //introduceParticlesAtWall();
        introduceParticlesInDomain(1.3);

        setGravity(Vec3D(0,0,-9.8));
    }

    void introduceSingleParticle(Vec3D p) {
                /* Simple run settings
         * Nx*Ny*Nz particles are created evenly spaced between [xmin,xmax]*[ymin,ymax]*[zmin,zmax] and checked for contact with the screw
         */
        const ParticleSpecies* s = speciesHandler.getObject(0);
        particleHandler.clear();
        SphericalParticle p0;
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
        SphericalParticle p0;
        p0.setSpecies(s);
        p0.setVelocity(Vec3D(0.0, 0.0, 0.0));
        p0.setRadius(2.0*particleRadius);
        SphericalParticle p1;
        p1.setSpecies(s);
        p1.setVelocity(Vec3D(0.0, 0.0, 0.0));
        p1.setRadius(particleRadius);

        //number of particles that fit in domain
        Mdouble distance;
        Vec3D normal;
        Vec3D p;
        Mdouble minDistance;
        unsigned counter = 0;
        for (p.X = getXMin() + particleRadius; p.X < getXMax(); p.X += 2.0 * particleRadius)
            for (p.Y = getYMin() + particleRadius; p.Y < getYMax(); p.Y += 2.0 * particleRadius)
                for (p.Z = getZMin() + particleRadius; p.Z < fillHeight_; p.Z += 2.0 * particleRadius)
                {
                    minDistance = p0.getRadius();
                    p0.setPosition(p);
                    for (auto w : wallHandler) {
                        //if touching the wall
                        if (w->getDistanceAndNormal(p0, distance, normal) && distance<minDistance)
                        {
                            minDistance=distance;
                            if (distance<0) break; 
                            p1.setPosition(p0.getPosition()+(distance-1.0001*particleRadius)*normal);
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
        SphericalParticle p0;
        p0.setSpecies(s);
        p0.setVelocity(Vec3D(0.0, 0.0, 0.0));
        p0.setRadius(particleRadius);

        //CHANGED BY BERT need ~ 1.5 times cubic packing fraction to get HCP packing, assume particles are more HCP packed than cubic
        Mdouble distance;
        Vec3D normal;
        Vec3D p;
        Mdouble minDistance;
        unsigned counter = 0;
        for (p.X = getXMin() + particleRadius; p.X < getXMax(); p.X += 2.0 * particleRadius)
            for (p.Y = getYMin() + particleRadius; p.Y < getYMax(); p.Y += 2.0 * particleRadius)
                for (p.Z = getZMin() + particleRadius; p.Z < fillHeight_; p.Z += 2.0 * particleRadius) //Changed Cubic to HCP here
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

    void printTime () const override {
        std::cout << "t=" << std::setprecision(3) << std::left << std::setw(6) << getTime()
            << ", tmax=" << std::setprecision(3) << std::left << std::setw(6) << getTimeMax()
            << ", EneRatio=" << std::setprecision(3) << std::left << std::setw(6) << getKineticEnergy()/getElasticEnergy()
            << std::endl;
    }

    void actionsAfterTimeStep () override {
        if (screw) screw->move_time(getTimeStep());
    }

    void writeScript() {
        logger(INFO,"Writing matlab script % to color the particles",getName() + ".m");
        helpers::writeToFile(getName() + ".m", "cd " + helpers::getPath() + "\n"
                "%% read in first file, to get the initial positions\n"
                "f = fopen('" + getName() + "_7.vtu');\n"
                "% header\n"
                "line = textscan(f,'%s %s %s %s %s %s %s',1,'Delimiter','\\n');\n"
                "% number of particles\n"
                "N = textscan(line{5}{1}(24:end),'%d',1); N=N{1};\n"
                "% positions\n"
                "P = textscan(f,'%f %f %f',N);\n"
                "%scatter(P{1},P{2})\n"
                "fclose(f);\n"
                "%% define a new speciesIndex, based on position, to color particles\n"
                "index = 1000*P{1};\n"
                "%% read in second file, a write out again with modified index\n"
                "f = fopen('" + getName() + "_250.vtu');\n"
                "g = fopen('Particle.vtu','w');\n"
                "% header\n"
                "line = textscan(f,'%s',3*N+15,'Delimiter','\\n');\n"
                "for i=1:length(line{1}), fprintf(g,'%s\\n',line{1}{i}); end\n"
                "% i/o indSpecies\n"
                "textscan(f,'%f',N,'Delimiter','\\n');\n"
                "fprintf(g,'%f\\n',index);\n"
                "% footer\n"
                "line = textscan(f,'%s','Delimiter','\\n');\n"
                "for i=1:length(line{1}), fprintf(g,'%s\\n',line{1}{i}); end\n"
                "fclose(f);\n"
                "fclose(g);");
    }
};

#endif
