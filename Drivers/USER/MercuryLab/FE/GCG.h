#ifndef MERCURY_GCG_H
#define MERCURY_GCG_H

#include <random>
#include <Walls/Screw.h>
#include "Mercury3D.h"
#include "Walls/AxisymmetricIntersectionOfWalls.h"
#include "Walls/BasicIntersectionOfWalls.h"
using constants::pi;
using mathsFunc::cubic;

class GCG : public Mercury3D
{
public:

    /**
     * All fixed parameters (defining the geometry) are set in the constructor; the variable parameters
     * (particleSpecies, particleWallSpecies, inflow rate, rotation rate, particle radius, polydispersity)
     * are passed into the constructor.
     *
     * The dimensions of the system are:
     *  - |y|,|z|<cylinderRadius
     *  - 0 < z < mixerLength,
     * The inflow is at z=caplength, outflow at z=mixerLength-caplength.
     *
     * The geometry is split into several walls:
     *  - mixer: the cylinder forming the outer mixer wall
     *    parameters: mixerLength, capLength, mixerRadius
     *  - inflow/outflow: the cylinders forming the in- and outflow
     *    implemented as an intersection with the mixer
     *    parameters: inflowRadius, inflowHeight, outflowHeight
     *  - axis: cylinder forming the inner beam on which the blades are attached
     *    parameters: axisRadius
     *  - inflowScrew, outflowScrew: the screwd below the inflow/ above the outflow
     *    parameters: bladeRadius, screwPitch, inflowScrewLength, outflowScrewLength
     *  - blade: the blades between the screws, positioned at a uniform distance from each other, and half that distance from the screws
     *    parameters: bladeRadius, bladeNumber, bladeInclination, bladeOpening
     *
     *  @param particleSpecies species describing the contact between two particles
     *  @param particleWallSpecies     species describing the contact between a particle and a wall
     *  @param fillRate        rate (in m^3/s) of particles to be inserted
     *  @param rpm             rotation rate of the mixer (in revolutions per minute)
     *  @param particleRadius  mean particle radius (in m)
     *  @param polydispersity  relative standard deviation of the log of the particle radius, assuming log-normal distribution
     *  @param batchFill true/false turns on/off batch filling - if set to false, particles all emptied in one step
	     */
    GCG(Mdouble maxFillRate, Mdouble rpm, Mdouble particleRadius, Mdouble polydispersity,
            ParticleSpecies* particleSpecies, BaseSpecies* particleWallSpecies)
    : maxFillRate_(maxFillRate), particleRadius_(particleRadius), polydispersity_(polydispersity)
    {
        logger(INFO,"Mean particle radius: % +/- % m", particleRadius, polydispersity*particleRadius);
        logger(INFO,"Max fill rate: % l/s", maxFillRate_*1e3);

        //random.randomise();

        // angular velocity of axis, screw and blades
        const Vec3D angularVelocity = {2.*pi*rpm/60.,0,0};

        // BaseInteractable::rotate(quarterturn) rotates an object by 90 degrees (¡empirical value, not exact!)
        const Vec3D quarterTurn = {2,0,0};

        //set name, gravity, dimensions
        setName("GCG");
        setGravity({0,0,-9.8});
        setDomain({0,-mixerRadius,-outflowHeight},{mixerLength,mixerRadius,inflowHeight});
        setXBallsAdditionalArguments("-s 5.5 -oh 0 -o 100 -h 530 -w 1500");

        //add species
        //use this pointer to set the species of the particles
        auto pSpecies = speciesHandler.copyAndAddObject(particleSpecies);
        //use this pointer to set the species of the walls
        auto wSpecies = speciesHandler.copyAndAddObject(particleSpecies);
        speciesHandler.getMixedObject(0,1)->mixAll(particleWallSpecies,particleWallSpecies);

        // add mixer
        AxisymmetricIntersectionOfWalls mixer;
        mixer.setSpecies(wSpecies);
        mixer.setPosition({0,0,0});
        mixer.setAxis({1,0,0});
        mixer.addObject({1,0,0},{mixerRadius,0,0});
        //mixer.addObject({0,0,1},{0.5*mixerRadius,0,0});
        //wallHandler.copyAndAddObject(mixer);

        // add caps
        InfiniteWall cap;
        cap.setSpecies(wSpecies);
        cap.set({1,0,0},{mixerLength,0,0});
        AxisymmetricIntersectionOfWalls capRender;
        capRender.setSpecies(wSpecies);
        capRender.setAxis({1,0,0});
        capRender.addObject({0,0,1},{0,0,mixerLength});
        capRender.addObject({-1,0,0},{mixerRadius,0,0});
        cap.addRenderedWall(capRender.copy());
        wallHandler.copyAndAddObject(cap);

        // add caps
        InfiniteWall cap2;
        cap2.setSpecies(wSpecies);
        cap2.set({-1,0,0},{0,0,0});
        AxisymmetricIntersectionOfWalls capRender2;
        capRender2.setSpecies(wSpecies);
        capRender2.setAxis({1,0,0});
        capRender2.addObject({0,0,-1},{0,0,0});
        capRender2.addObject({-1,0,0},{mixerRadius,0,0});
        cap2.addRenderedWall(capRender2.copy());
        wallHandler.copyAndAddObject(cap2);

        // add hopper above inflow
        AxisymmetricIntersectionOfWalls hopper;
        hopper.setSpecies(wSpecies); //do we need species here?
        hopper.setPosition({capLength,0,0});
        hopper.setAxis({0,0,1});
        // conical wall
        hopper.addObject({1,0,-1},{inflowRadius,0,mixerRadius});
        // other side of conical wall
        hopper.addObject({0,0,1},{0,0,mixerRadius});
        AxisymmetricIntersectionOfWalls hopperRender = hopper;
        hopperRender.addObject({-1,0,1.1},{1.1*inflowRadius,0,mixerRadius});
        hopper.addRenderedWall(hopperRender.copy());
        wallHandler.copyAndAddObject(hopper);

        // add inflow
        AxisymmetricIntersectionOfWalls inflow;
        inflow.setSpecies(wSpecies); //do we need species here?
        inflow.setPosition({capLength,0,0});
        inflow.setAxis({0,0,1});
        // cylindrical wall
        inflow.addObject({1,0,0},{inflowRadius,0,0});
        // end of cylinder at z=0
        inflow.addObject({0,0,1},{0,0,0});
        // end of cylinder at z=0
        inflow.addObject({0,0,-1},{0,0,mixerRadius});
        AxisymmetricIntersectionOfWalls inflowRender = inflow;
        //render only above z=0.8*mixerHeight
        inflowRender.addObject({0,0,1},{0,0,0.8*mixerRadius});
        //render a wall thickness of 0.1*inflowRadius
        inflowRender.addObject({-1,0,0},{1.1*inflowRadius,0,0});

        // add outflow
        AxisymmetricIntersectionOfWalls outflow;
        outflow.setSpecies(wSpecies); //do we need species here?
        outflow.setPosition({mixerLength-capLength,0,0});
        outflow.setAxis({0,0,1});
        outflow.addObject({1,0,0},{inflowRadius,0,0});
        outflow.addObject({0,0,-1},{0,0,0});
        AxisymmetricIntersectionOfWalls outflowRender = outflow;
        outflowRender.addObject({0,0,-1},{0,0,-0.8*mixerRadius});
        outflowRender.addObject({-1,0,0},{1.01*inflowRadius,0,0});

        //intersect mixer and inflow
        BasicIntersectionOfWalls mixerInflow;
        mixerInflow.setSpecies(wSpecies);
        mixerInflow.add(mixer);
        mixerInflow.add(inflow);
        //mixerInflow.addRenderedWall(&mixer);
        mixerInflow.addRenderedWall(inflowRender.copy());
        wallHandler.copyAndAddObject(mixerInflow);

        //intersect mixer and outflow
        BasicIntersectionOfWalls mixerOutflow;
        mixerOutflow.setSpecies(wSpecies);
        mixerOutflow.add(mixer);
        mixerOutflow.add(outflow);
        mixerOutflow.addRenderedWall(outflowRender.copy());
        wallHandler.copyAndAddObject(mixerOutflow);

        // add axis
        AxisymmetricIntersectionOfWalls axis;
        axis.setSpecies(wSpecies);
        axis.setPosition({0,0,0});
        axis.setAxis({1,0,0});
        axis.addObject({-1,0,0},{axisRadius,0,0});
        axis.setAngularVelocity(angularVelocity);
        wallHandler.copyAndAddObject(axis);

        //Screw(Vec3D start, Mdouble l, Mdouble r, Mdouble n, Mdouble omega, Mdouble thickness);
        Screw inflowScrew({0,0,0}, inflowScrewLength, bladeRadius, -inflowScrewLength/screwPitch, rpm/60.0, 0);
        inflowScrew.setSpecies(wSpecies);
        inflowScrew.setPosition({0,0,0});
        inflowScrew.setOrientationViaNormal({1,0,0});
        wallHandler.copyAndAddObject(inflowScrew);

        Screw outflowScrew({0,0,0}, outflowScrewLength, bladeRadius, outflowScrewLength/screwPitch, -rpm/60.0, 0);
        outflowScrew.setSpecies(wSpecies);
        outflowScrew.setPosition({mixerLength,0,0});
        outflowScrew.setOrientationViaNormal({-1,0,0});
        wallHandler.copyAndAddObject(outflowScrew);

        Mdouble bladePitch = 8.0*bladeInclination*bladeRadius;
        Mdouble bladeLength = bladeOpening*bladePitch/(2.0*pi);
        Mdouble bladeDistance = (mixerLength-inflowScrewLength-outflowScrewLength)/bladeNumber;
        logger(INFO,"p% l% n%",bladePitch,bladeLength,bladeLength/bladePitch);
        Screw bladeUp({0,0,0}, bladeLength, bladeRadius, -bladeLength/bladePitch, -rpm/60.0, 0);
        bladeUp.setSpecies(wSpecies);
        bladeUp.setPosition({mixerLength-outflowScrewLength-0.5*bladeDistance,0,0});
        bladeUp.setOrientationViaNormal({-1,0,0});
        bladeUp.setAngularVelocity({0,0,0});

        Screw bladeDown = bladeUp;
        bladeDown.rotate(quarterTurn);
        bladeDown.rotate(quarterTurn);

        Screw bladeLeft = bladeUp;
        bladeLeft.move({-bladeDistance,0,0});
        bladeLeft.rotate(quarterTurn);

        Screw bladeRight = bladeUp;
        bladeRight.move({-bladeDistance,0,0});
        bladeRight.rotate(-quarterTurn);

        for (unsigned i=0; i<bladeNumber; i+=2) {
            wallHandler.copyAndAddObject(bladeUp);
            wallHandler.copyAndAddObject(bladeDown);
            bladeUp.move({-2.0*bladeDistance,0,0});
            bladeDown.move({-2.0*bladeDistance,0,0});
        }
        for (unsigned i=1; i<bladeNumber; i+=2) {
            wallHandler.copyAndAddObject(bladeLeft);
            wallHandler.copyAndAddObject(bladeRight);
            bladeLeft.move({-2.0*bladeDistance,0,0});
            bladeRight.move({-2.0*bladeDistance,0,0});
        }

        logger(INFO,"\nAdded % walls", wallHandler.getNumberOfObjects());

        //if (maxFillRate_==0)
            insertParticles();
        logger(INFO,"\nInitially added % particles", particleHandler.getNumberOfObjects());
    }

    void printTime() const override
    {
        std::cout << "t " << std::setprecision(3) << std::left << std::setw(6) << getTime()
                  << " tMax " << std::setprecision(3) << std::left << std::setw(6) << getTimeMax()
                  << " N " << std::setprecision(3) << std::left << std::setw(6) << particleHandler.getNumberOfObjects()
                  << " FillRate [ml/s] " << std::setprecision(3) << std::left << std::setw(6) << insertedVolume/getTime()*1e6
                  << std::endl;
    }

    Mdouble generateRadius() {
        //set up normal distribution
        static std::random_device rd;
        static std::mt19937 gen(rd());
        static std::normal_distribution<> d(particleRadius_, polydispersity_*particleRadius_);
        static const Mdouble radiusMin = (1.-1.5*polydispersity_)*particleRadius_;
        static const Mdouble radiusMax = (1.+1.5*polydispersity_)*particleRadius_;
        return std::max(std::min(d(gen),radiusMax),radiusMin);
    }

    //remove particles according to (max.) outflow rate
    void insertParticles()
    {
        //add particles
        BaseParticle p;
        p.setSpecies(speciesHandler.getObject(0));
        p.setRadius(generateRadius());
        unsigned maxFailed = 200;
        unsigned failed = 0;
        while (failed<=maxFailed) {
            const Mdouble x = random.getRandomNumber(getXMin(),getXMax());
            const Mdouble y = 1.2*random.getRandomNumber(getYMin(),getYMax());
            const Mdouble z = 1.2*random.getRandomNumber(getZMin(),getZMax());
            p.setPosition({x, y, z});
            // check if particle can be inserted
            if (checkParticleForInteraction(p)) {
                particleHandler.copyAndAddObject(p);
                p.setRadius(generateRadius());
                failed = 0;
            } else failed++;
        }
        logger(INFO,"Inserted % particles",particleHandler.getNumberOfObjects());
    }


    //remove particles according to (max.) outflow rate
    void actionsAfterTimeStep() override
    {
        // Functions after this statement only get executed every 100th timestep (if counter==100)
        static unsigned counter = 0;
        if (++counter!=100) return;
        else counter = 0;

        //add particles
        BaseParticle p;
        p.setSpecies(speciesHandler.getObject(0));
        p.setRadius(generateRadius());
        unsigned maxFailed = 100;
        unsigned failed = 0;
        while (insertedVolume<=maxFillRate_*getTime() && failed<=maxFailed) {
            const Mdouble r = inflowRadius+(inflowHeight-mixerRadius)-p.getRadius();
            const Mdouble x = capLength + random.getRandomNumber(-r,r);
            const Mdouble y = random.getRandomNumber(-r,r);
            const Mdouble z = inflowHeight * random.getRandomNumber(1,1.2);
            p.setPosition({x, y, z});
            // check if particle can be inserted
            if (checkParticleForInteraction(p)) {
                particleHandler.copyAndAddObject(p);
                insertedVolume += p.getVolume();
                p.setRadius(generateRadius());
            } else failed++;
        }

        //delete particles
        particleHandler.removeIf([this](BaseParticle* p){return p->getPosition().Z<-outflowHeight;});
        //logger(INFO,"Time %: Inserted % particles",getTime(), particleHandler.getNumberOfObjects());
    }

private:
    const Mdouble maxFillRate_;
    const Mdouble particleRadius_;
    const Mdouble polydispersity_;
    Mdouble insertedVolume = 0;

    //Fixed parameters
    const Mdouble mixerLength = 636.5e-3; //distance between in- and outflow
    const Mdouble capLength = 0.5*(mixerLength-495e-3); //length of caps, i.e. mixer length beyond in and outflow
    const Mdouble mixerRadius = 50e-3; //radius of the mixer
    const Mdouble inflowRadius = 30e-3; //radius of the in- and outflow
    const Mdouble inflowHeight = mixerRadius-inflowRadius+sqrt(3)*inflowRadius; //determines where the particles flow in
    // (the value in the sqrt is the factor by which the density is increased in the hopper)
    const Mdouble outflowHeight = 1.2*mixerRadius; //radial distance of outflow from mixer axis
    const Mdouble axisRadius = 0.4*mixerRadius; //radius of the in- and outflow
    const Mdouble bladeRadius = 0.95*mixerRadius; //radius of the in- and outflow
    const Mdouble screwPitch = (capLength+inflowRadius)/2; //radius of the screws at the in- and outflow
    const Mdouble inflowScrewLength = 2.5*screwPitch; //length of the inflow screw
    const Mdouble outflowScrewLength = 1.5*screwPitch; //length of the outflow screw
    const unsigned bladeNumber = 7; //number of blades
    const Mdouble bladeInclination = 20.0*pi/180.0; //inclination angle between the blade and the yz-plane
    const Mdouble bladeOpening = 60.0*pi/180.0; //opening angle of the arc-shaped blade

};

#endif //MERCURY_GCG_H
