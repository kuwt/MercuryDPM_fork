#ifndef MERCURY_CSK_H
#define MERCURY_CSK_H

#include "Mercury3D.h"
#include "Species/LinearViscoelasticSlidingFrictionSpecies.h"
#include "Walls/AxisymmetricIntersectionOfWalls.h"
#include "Walls/BasicIntersectionOfWalls.h"
using constants::pi;
using mathsFunc::cubic;

class CSK : public Mercury3D
{
public:

    /**
     * All fixed parameters (defining the geometry) are set in the constructor; the variable parameters
     * (particleSpecies, particleWallSpecies, fill volume, outflow rate, rotation rate, particle radius, polydispersity)
     * are passed into the constructor.
     *
     * The dimensions of the system are:
     *  - |x|,|y|<cylinderRadius
     *  - (-outflowHeight) < z < mixerHeight,
     * i.e. the origin is at the tip of the cone forming the base of the mixer
     *
     * The geometry is split into several walls:
     *  - cylinder: the cylinder forming part of the outer mixer wall
     *    parameters: cylinderRadius: 900 mm
     *  - cone: the cone forming the base of the outer mixer wall;
     *    parameters: coneInclination: 30 deg
     *  - outflow: the cylinder forming the outflow at the base
     *    implemented only if outflowRate>0
     *    implemented as an intersection with the cone
     *    \todo there is a second part of the outflow that I do not understand
     *    parameters: outflowRadius: 102 mm, outflowCenter (x-position of the cylindrical axis): 130 mm
     *  - beam: cylinder forming the  inner beam onto which the blades are attached
     *    parameters: beamRadius: 70 mm
     *  - anchor1: the blade at the base of the mixer, initially oriented in +x direction
     *    parameters:
     *    . anchorLength, anchorWidth, anchorThick (dimensions of rectangular box): 900 mm x 26 mm x 170 mm,
     *    . anchorOffset (shift in angular direction): 19 deg (where center of blade intersects the central beam)
     *    . anchorAngle (radial rotation, in yz-plane): 45 deg
     *    . blade is parallel to the cone surface (no parameter, hard-coded)
     *    . anchorConeDistance (space between anchor and cone): 30.67 mm
     *    . Note: The upmost edge of the blade is chipped off by about 30.14 mm (not implemented)
     *  - anchor2: same as anchor1, but rotated by 180 degree
     *  - cAnchor1: the rectangular blade at the center of the mixer, initially oriented in +y direction
     *    parameters:
     *    . cAnchorLength, cAnchorWidth, cAnchorThick (dimensions of rectangular box): 692 mm x 20 mm x 255 mm,
     *    . radial rotation in xz-plane: -45 deg, i.e. in opposite direction to the anchor blade (no parameter, hard-coded)
     *    . anchorConeDistance (space between anchor and cone): 30.67 mm
     *    . the upmost edge of the blade is chipped off by about 30.14 mm (not implemented)
     *    . cAnchorHeight: (vertical distance from anchorOffset.Z): (75 + 850 + 90) = ~1015 mm
     *  - cAnchor2: same as cAnchor1, but rotated by 180 degree
     *  - guillo1: the "guillotine" blade attached to cAnchor1
     *    . Axisymmetric wall forming a torus piece between
     *      z = anchorOffset.Z+cAnchorHeight pm guilloThick/2, 0 < r-guilloRadius < guilloWidth, and theta = pm guilloAngle/2
     *    . Rotated radially by -45 degree (i.e. opposite to the anchor blade)
     *    . guilloThickness: 20 mm
     *    . guilloRadius: 692-20 = 672 mm
     *    . guilloWidth: 170 mm
     *    . guilloAngle = 69 * pi/180 rad
     *  - guillo2: same as guillo1, but rotated by 180 degree
     *  - disk1: a small disk between the guillo1 and cAnchor1
     *    parameters: diskRadius: 85 mm, diskThick: 20 mm
     *  - disk2: same as disk1, but rotated by 180 degree
     *  - case: a small disk around the beam where anchor1 is placed
     *    parameters: caseLength: 150 mm, caseRadius: 22 mm + beamRadius
     *  - cCase: a small disk around the beam where cAnchor1 is placed
     *    parameters: cCaseLength: 180 mm, same radius as case
     *
     *  @param particleSpecies species describing the contact between two particles
     *  @param particleWallSpecies     species describing the contact between a particle and a wall
     *  @param fillVolume      absolute volume (not bulk volume, in m^3) of particles to be inserted
     *  @param outflowRate     rate (in m^3/s) at which particles are deleted at the end of the outflow pipe
     *                         (in absolute volume)
     *  @param rpm             rotation rate of the mixer (in revolutions per minute)
     *  @param particleRadius  mean particle radius (in m)
     *  @param polydispersity  false/true true turns on the true size distribution. Note it also turns on the correct hand filling process
     *  @param batchFill true/false turns on/off batch filling - if set to false, particles all emptied in one step
	     */
    CSK(ParticleSpecies* particleSpecies, BaseSpecies* particleWallSpecies, Mdouble fillVolume, Mdouble outflowRate,
        Mdouble rpm, Mdouble particleRadius, bool delayedStartOfRotation=true, bool polydispersity = false, bool batchFill = true, Mdouble bagsPerRevolution = 1.0)
     : bagsEmptied_(0), outflowRate_(outflowRate), outflowVolume_(0), delayedStartOfRotation_(delayedStartOfRotation), rpm_(rpm), particleRadius_(particleRadius),batchFill_(batchFill),bagsPerRevolution_(bagsPerRevolution)
    {
        logger(INFO,"Mean particle radius: % m", particleRadius);
        logger(INFO,"Polydispersity: %", polydispersity);
        logger(INFO,"Fill volume: % m^3", fillVolume);
        logger(INFO,"DelayedStartOfRotation: %", delayedStartOfRotation);
        logger(INFO,"Polydispersity: %", polydispersity);
        logger(INFO,"BatchFill: %", batchFill);
        logger(INFO,"BagsPerRevolution: %", bagsPerRevolution);

        //random.randomise();

        //Fixed parameters
        const Mdouble mixerHeight = 2040e-3; //upper end of the domain
        const Mdouble outflowHeight = 200e-3; //lower end of the domain
        const Mdouble cylinderRadius = 900e-3; //radius of the domain
        const Mdouble coneAngle = 30.*pi/180.;
        const Mdouble outflowRadius = 102e-3;
        const Mdouble outflowCenter = 130e-3;
        const Mdouble beamRadius = 50e-3;
        const Mdouble anchorLength = 900e-3;
        const Mdouble anchorWidth = 170e-3;
        const Mdouble anchorThick = 26e-3;
        const Mdouble anchorAngle = 45.*pi/180.;
        const Mdouble anchorOffsetAngle = -0*pi/180; //should have been 19 deg, why not??
        const Mdouble anchorConeDistance = 30.67e-3;
        const Mdouble cAnchorHeight = 1015e-3+170e-3;
        const Mdouble cAnchorLength = 692e-3;
        const Mdouble cAnchorWidth = 255e-3;
        const Mdouble cAnchorThick = 20e-3;
        const Mdouble guilloRadius = 918e-3;
        const Mdouble guilloRadialDistance = 672e-3;
        const Mdouble guilloAngle = 69.*pi/180.;
        const Mdouble guilloWidth = 170e-3;
        const Mdouble guilloThick = 20e-3;
        const Mdouble diskRadius = 85e-3;
        const Mdouble diskThick = 20e-3;
        const Mdouble caseLength = 150e-3;
        const Mdouble cCaseLength = 180e-3;
        const Mdouble caseRadius = 22e-3 + beamRadius;

        // BaseInteractable::rotate(quarterturn) rotates an object by 90 degrees (¡empirical value, not exact!)
        const Vec3D quarterTurn = {0,0,.63662*pi};

        // don't rotate if delayed rotation is set to true
        // (a non-zero value is used as an indicator of which walls should be rotated )
        const Vec3D angularVelocity = {0,0,delayedStartOfRotation_?1e-20:(rpm*pi/30.)};

        //set name, gravity, dimensions
        setName("CSK");
        setGravity({0,0,-9.8});
        setMax({cylinderRadius,cylinderRadius,mixerHeight});
        setMin({-cylinderRadius,-cylinderRadius,-outflowHeight});
        setXBallsAdditionalArguments("-solidf -v0 -w 1200 -h 800 -s 1 -moh 500");

	
        //add species
        //use this pointer to set the species of the particles
        auto pSpecies = speciesHandler.copyAndAddObject(particleSpecies);
        //use this pointer to set the species of the walls
        auto wSpecies = speciesHandler.copyAndAddObject(particleSpecies);
        speciesHandler.getMixedObject(0,1)->mixAll(particleWallSpecies,particleWallSpecies);

        // add cylinder
        AxisymmetricIntersectionOfWalls cylinder;
        cylinder.setSpecies(wSpecies);
        cylinder.setPosition({0,0,0});
        cylinder.setAxis({0,0,1});
        cylinder.addObject({1,0,0},{cylinderRadius,0,0});
        wallHandler.copyAndAddObject(cylinder);

        // add cone
        AxisymmetricIntersectionOfWalls cone;
        cone.setSpecies(wSpecies);
        cone.setPosition({0,0,0});
        cone.setAxis({0,0,1});
        cone.addObject({sin(coneAngle),0,-cos(coneAngle)},{0,0,0});
	if (outflowRate==0.0) {
            wallHandler.copyAndAddObject(cone);
        } else {
            // add outflow
            AxisymmetricIntersectionOfWalls outflow;
            outflow.setSpecies(wSpecies); //do we need species here?
            outflow.setPosition({outflowCenter,0,0});
            outflow.setAxis({0,0,1});
            outflow.addObject({1,0,0},{outflowRadius,0,0});
            //intersect cone and outflow
            BasicIntersectionOfWalls coneOutflow;
            coneOutflow.setSpecies(wSpecies);
            coneOutflow.add(cone);
            coneOutflow.add(outflow);
            wallHandler.copyAndAddObject(coneOutflow);

            // put wall at bottom of outflow
            InfiniteWall base;
            base.setSpecies(wSpecies);
            base.set({0,0,-1},{0,0,-outflowHeight});
            wallHandler.copyAndAddObject(base);
        }

        // add anchor1
        IntersectionOfWalls anchor1;
        anchor1.setSpecies(wSpecies);
        // rotated by anchorAngle in yz plane, then by coneAngle in the xz-plane
        // normalLength=(1,0,0) for un-inclined blade
        Vec3D normalLength = {cos(coneAngle),0,sin(coneAngle)};
        // normalWidth=(0,1,0) for un-inclined blade
        Vec3D normalWidth  = {-sin(anchorAngle)*sin(coneAngle),cos(anchorAngle),sin(anchorAngle)*cos(coneAngle)};
        // normalThick=(0,0,1) for un-inclined blade
        Vec3D normalThick  = {-cos(anchorAngle)*sin(coneAngle),-sin(anchorAngle),cos(anchorAngle)*cos(coneAngle)};
        // position of the blade: centered in width and thickness, leftmost in length
        Vec3D anchorOffset = beamRadius*Vec3D(cos(anchorOffsetAngle),sin(anchorOffsetAngle),cos(anchorOffsetAngle)*tan(coneAngle));
        Vec3D pos = anchorOffset + Vec3D(0,0,anchorConeDistance/cos(coneAngle) + anchorWidth*sin(anchorAngle));
        anchor1.addObject(-normalLength,pos+anchorLength*normalLength);
        anchor1.addObject(normalLength,pos);
        anchor1.addObject(-normalWidth,pos+0.5*anchorWidth*normalWidth);
        anchor1.addObject(normalWidth,pos-0.5*anchorWidth*normalWidth);
        anchor1.addObject(-normalThick,pos+0.5*anchorThick*normalThick);
        anchor1.addObject(normalThick,pos-0.5*anchorThick*normalThick);
        anchor1.setAngularVelocity(angularVelocity);
        wallHandler.copyAndAddObject(anchor1);

        auto anchor2 = wallHandler.copyAndAddObject(anchor1);
        anchor2->rotate(quarterTurn);
        anchor2->rotate(quarterTurn);

        // add beam
        AxisymmetricIntersectionOfWalls beam;
        beam.setSpecies(wSpecies);
        beam.setPosition({0,0,pos.Z});
        beam.setAxis({0,0,1});
        beam.addObject({-1,0,0},{beamRadius,0,0});
        beam.addObject({0,0,1},{0,0,0});
        beam.setAngularVelocity(angularVelocity);
        wallHandler.copyAndAddObject(beam);

        // create case_ (case is a C++ keyword, thus cannot be used)
        AxisymmetricIntersectionOfWalls case_;
        case_.setSpecies(wSpecies);
        case_.setAxis({0,0,1});
        case_.setPosition({0,0,pos.Z});
        case_.addObject({0,0,-1},{0,0,0.5*caseLength});
        case_.addObject({0,0,1},{0,0,-0.5*caseLength});
        case_.addObject({-1,0,0},{caseRadius,0,0});
        case_.setAngularVelocity(angularVelocity);
        wallHandler.copyAndAddObject(case_);

        // create cCase
        AxisymmetricIntersectionOfWalls cCase;
        cCase.setSpecies(wSpecies);
        cCase.setAxis({0,0,1});
        cCase.setPosition({0,0,anchorOffset.Z+cAnchorHeight});
        cCase.addObject({0,0,-1},{0,0,0.5*cCaseLength});
        cCase.addObject({0,0,1},{0,0,-0.5*cCaseLength});
        cCase.addObject({-1,0,0},{caseRadius,0,0});
        cCase.setAngularVelocity(angularVelocity);
        wallHandler.copyAndAddObject(cCase);

        // add cAnchor1
        IntersectionOfWalls cAnchor1;
        cAnchor1.setSpecies(wSpecies);
        // rotated by -anchorAngle in xz plane
        // rotated by -45 deg in xz plane
        normalLength = Vec3D::getUnitVector({0,1,0});
        normalWidth  = Vec3D::getUnitVector({1,0,1});
        normalThick  = Vec3D::getUnitVector({1,0,-1});
        // position of the blade: centered in width and thickness, leftmost in length
        pos = {0,0,anchorOffset.Z+cAnchorHeight};
        cAnchor1.addObject(-normalLength,pos+cAnchorLength*normalLength);
        cAnchor1.addObject(normalLength,pos+0.5*beamRadius*normalLength);
        cAnchor1.addObject(-normalWidth,pos+0.5*cAnchorWidth*normalWidth);
        cAnchor1.addObject(normalWidth,pos-0.5*cAnchorWidth*normalWidth);
        cAnchor1.addObject(-normalThick,pos+0.5*cAnchorThick*normalThick);
        cAnchor1.addObject(normalThick,pos-0.5*cAnchorThick*normalThick);
        cAnchor1.setAngularVelocity(angularVelocity);
        wallHandler.copyAndAddObject(cAnchor1);

        auto cAnchor2 = wallHandler.copyAndAddObject(cAnchor1);
        cAnchor2->rotate(quarterTurn);
        cAnchor2->rotate(quarterTurn);

        // create torus of guillo1
        AxisymmetricIntersectionOfWalls torus1;
        torus1.setSpecies(wSpecies);
        // rotated by -45 deg in xz plane
        torus1.setAxis({1,0,1});
        // position of the guillotine: centered in width and thickness, leftmost in length
        torus1.setPosition({0,guilloRadialDistance-guilloRadius,anchorOffset.Z+cAnchorHeight});
        torus1.addObject({1,0,0},{guilloRadius,0,0});
        torus1.addObject({-1,0,0},{guilloRadius+guilloWidth,0,0});
        torus1.addObject({0,0,1},{0,0,-0.5*guilloThick});
        torus1.addObject({0,0,-1},{0,0,0.5*guilloThick});
        // create angleCut of guillo1
        IntersectionOfWalls angleCut1;
        angleCut1.setSpecies(wSpecies);
        // rotated by pm 45 in xz plane, then by guilloAngle in the xy-plane
        Mdouble sin45 = sqrt(0.5);
        // normalFront=(1,0,0) for un-inclined guillo
        Vec3D normalFront = {cos(0.5*guilloAngle)*sin45,sin(0.5*guilloAngle),-cos(5.5*guilloAngle)*sin45};
        // normalBack=(-1,0,0) for un-inclined guillo
        Vec3D normalBack = {-cos(0.5*guilloAngle)*sin45,sin(0.5*guilloAngle),cos(5.5*guilloAngle)*sin45};
        pos.Y += guilloRadialDistance-guilloRadius;
        angleCut1.addObject(normalFront,pos);
        angleCut1.addObject(normalBack,pos);
        //intersect angleCut and torus to create guillo1
        BasicIntersectionOfWalls guillo1;
        guillo1.setSpecies(wSpecies);
        guillo1.add(torus1);
        guillo1.add(angleCut1);
        //wallHandler.copyAndAddObject(angleCut1);
        //wallHandler.copyAndAddObject(torus1);
        guillo1.setAngularVelocity(angularVelocity);
        wallHandler.copyAndAddObject(guillo1);

        auto guillo2 = wallHandler.copyAndAddObject(guillo1);
        guillo2->rotate(quarterTurn);
        guillo2->rotate(quarterTurn);

        // create disk1
        AxisymmetricIntersectionOfWalls disk1;
        disk1.setSpecies(wSpecies);
        disk1.setAxis({0,1,0});
        disk1.setPosition({0,0,pos.Z});
        disk1.addObject({0,0,-1},{0,0,guilloRadialDistance+diskThick});
        disk1.addObject({0,0,1},{0,0,guilloRadialDistance});
        disk1.addObject({-1,0,0},{diskRadius,0,0});
        disk1.setAngularVelocity(angularVelocity);
        wallHandler.copyAndAddObject(disk1);

        auto disk2 = wallHandler.copyAndAddObject(disk1);
        disk2->rotate(quarterTurn);
        disk2->rotate(quarterTurn);



        BaseParticle p;
        p.setSpecies(pSpecies);
        if (polydispersity==false)
        {
            logger(INFO,"Creating particles for mono-dispersed distribution");

            //\todo TW: Anthony, the 1.1 polydispersity has been taken out; this will create chrystallisation artifacts, why?
            Mdouble rMin = particleRadius;
            Mdouble rMax = particleRadius;
            p.setRadius(rMax);
            Mdouble fillHeight = getZMin();
            while (fillVolume>0)
            {
                Mdouble x = random.getRandomNumber(getXMin(), getXMax());
                Mdouble y = random.getRandomNumber(getYMin(), getYMax());
                Mdouble z = random.getRandomNumber(getZMin(), fillHeight);
                p.setPosition({x, y, z});
                // check if particle can be inserted
                if (checkParticleForInteraction(p))
                {
                    particleHandler.copyAndAddObject(p);
                    fillVolume -= p.getVolume();
                    p.setRadius(random.getRandomNumber(rMin,rMax));
                    //double pol = random.getRandomNumber(-polydispersity,polydispersity);
                    //p.setRadius(radius*(1+pol));
                    if (particleHandler.getNumberOfObjects()%100==0) std::cout << '.' << std::flush;
                    if (particleHandler.getNumberOfObjects()%1000==0) std::cout << ' ';
                    if (particleHandler.getNumberOfObjects()%10000==0) std::cout << fillVolume << '\n';
                }
                else
                {
                    fillHeight += 0.01*particleRadius;// the lower this number, the lower the initial filling height,
                    // but the insertion might take longer. Thus, increase this number, if the insertion takes too long.
                }
            }
        }
        else //You are here if it is polydispersited (this is Anthony's superhack)
        {
            //This is crap code (as it is full of magic numbers) for FC and could be generised and hence used for a lot of problems

            logger(INFO,"You are creating a polydispersed mix");

            //add a second species but convert non. We will need a very strange combination of species for the true CSK filling.

            auto p2Species = speciesHandler.copyAndAddObject(particleSpecies);
            speciesHandler.getMixedObject(0,1)->mixAll(particleWallSpecies,particleWallSpecies);
            speciesHandler.getMixedObject(1,2)->mixAll(particleSpecies, particleSpecies);


        }
        logger(INFO,"\nInserted % particles", particleHandler.getNumberOfObjects());

    }

    // Constructor used for restarting at mixing stage (stage 2)
    // By setting the second argument to true, can also just restart a current mixing file
    CSK(std::string name, bool justRestart = false) : outflowRate_(0), particleRadius_(0) {
        setName(name);
        readRestartFile();
        logger(INFO,"Restarted %",getName());
        if (justRestart == false) {
	    setName("Mix");
            setTime(0.0);
	    }
        //note, this code does not insert particles because batchFill_ is uninitialized.
    }

    // Constructor used for restarting for stage 3 (outflow) from stage 2 (mixing)
    //second argument determines outflow rate
    CSK(std::string name, Mdouble rate) : outflowRate_(rate), particleRadius_(0) {
        setName(name);
        readRestartFile();
        logger(INFO,"Restarted %",getName());
        setName("Empty");
        setTime(0.0);

        const Mdouble outflowRadius = 102e-3;
        const Mdouble outflowCenter = 130e-3;
        const Mdouble coneAngle = 30.*pi/180.;
	    const Mdouble outflowHeight = 200e-3; //lower end of the domain
        
    	// add cone
	    AxisymmetricIntersectionOfWalls cone;
        cone.setSpecies(speciesHandler.getObject(1));
        cone.setPosition({0,0,0});
        cone.setAxis({0,0,1});
        cone.addObject({sin(coneAngle),0,-cos(coneAngle)},{0,0,0});
        // add outflow
        AxisymmetricIntersectionOfWalls outflow;
        outflow.setSpecies(speciesHandler.getObject(1)); //do we need species here?
        outflow.setPosition({outflowCenter,0,0});
        outflow.setAxis({0,0,1});
        outflow.addObject({1,0,0},{outflowRadius,0,0});
    	//intersect cone and outflow
	    BasicIntersectionOfWalls coneOutflow;
        coneOutflow.setSpecies(speciesHandler.getObject(1));
        coneOutflow.add(cone);
        coneOutflow.add(outflow);
        wallHandler.removeObject(1);
        wallHandler.copyAndAddObject(coneOutflow);

        // put wall at bottom of outflow
        InfiniteWall base;
        base.setSpecies(speciesHandler.getObject(1));
        base.set({0,0,-1},{0,0,-outflowHeight});
        wallHandler.copyAndAddObject(base);
    }



    // This is a function for testing that removes all particles and instead
    // creates a slice of particles in the xz plane that can be viewed with xballs.
    // Only one time step is executed, allowing us to view the walls, as they change the velocity of the particles.
    void wallTest(Mdouble r) {
        particleHandler.clear();
        BaseParticle p;
        p.setSpecies(speciesHandler.getObject(0));
        p.setRadius(r);
        for (Mdouble x=getXMin()+r; x<getXMax(); x+=r+r) {
            for (Mdouble z=getZMin()+r; z<getZMax(); z+=r+r) {
                p.setPosition({x,0,z});
                particleHandler.copyAndAddObject(p);
            }
        }
        setTimeStep(1e-10);
        setTimeMax(getTimeStep());
        setXBallsAdditionalArguments("-solid -v0 -w 1200 -h 800 -s 1 -moh 500 -cmaxset .01");
        logger(INFO,"\nInserted % particles in a plane", particleHandler.getNumberOfObjects());
    }



    // This is a function for testing that removes all particles and instead
    // creates a slice of particles in the xz plane that can be viewed with xballs.
    // Only one time step is executed, allowing us to view the walls, as they change the velocity of the particles.
    void guilloTest(Mdouble r) {
        particleHandler.clear();
        BaseParticle p;
        p.setSpecies(speciesHandler.getObject(0));
        p.setRadius(r);
        for (Mdouble x=-0.5+r; x<0.5; x+=r+r)
        {
            for (Mdouble y=getYMin()+r; y<-0.5; y+=r+r)
            {
                for (Mdouble z=0.75+r; z<1.7; z+=r+r)
                {
                    p.setPosition({x, y, z});
                    Mdouble dist;
                    Vec3D normal;
                    if (wallHandler.getObject(8)->getDistanceAndNormal(p,dist,normal))
                    {
                        p.setPosition(Vec3D(x, y, z)-(dist-0.001)*normal);
                        particleHandler.copyAndAddObject(p);
                    }
                    if (wallHandler.getObject(10)->getDistanceAndNormal(p,dist,normal))
                    {
                        p.setPosition(Vec3D(x, y, z)-(dist-0.001)*normal);
                        particleHandler.copyAndAddObject(p);
                    }
                    //else
                    //    logger(INFO,"\nskip particle");
                }
            }
        }
        setGravity({0,0,0});
        setTimeStep(1e-10);
        setTimeMax(getTimeStep());
        setXBallsAdditionalArguments("-solid -v0 -w 1200 -h 800 -s 1 -moh 500 -cmaxset .01");
        logger(INFO,"\nInserted % particles at guillo", particleHandler.getNumberOfObjects());
    }

    void printTime() const override
    {
        std::cout << "t " << std::setprecision(3) << std::left << std::setw(6) << getTime()
                  << " tMax " << std::setprecision(3) << std::left << std::setw(6) << getTimeMax()
         << " EneRatio " << std::setprecision(3) << std::left << std::setw(6) << getKineticEnergy()/getElasticEnergy()
         << " TotalMass " << std::setprecision(3) << std::left << std::setw(6) << particleHandler.getMass()
                  << std::endl;
    }

    //to get rid of batch process
    //comment line 439 (checkIfBag...) from
    //change the fill volume from 10kg 600kg
    //OR in actionsBeforeTimeLoop, take out "checkIfBagSHouldBeAdded" and replace with sequence:
    //e.g. addWhiteBag();
    //addWhiteBag();
    //addOrangeBag();
    //addOrangeBag();
    //addWhiteBag();
    //addWhiteBag();

    //alternatively in line ~590, change the upper and lower limits of the size distribution to narrow the total size distribution (see Thomas' notes on  skype)
    //i.e. in first line, change 0.4 to 0.25, so no particle can be smaller than 80% of the smallest ('2mm') size
    void actionsBeforeTimeLoop() override
    {
        if (batchFill_) {
            checkIfBagShouldBeAdded();
        } else {
            for (int i = 0; i <= 5; i++) {
                addWhiteBag();
            }
            for (int i = 0; i <= 5; i++) {
                addOrangeBag();
            }
            for (int i = 0; i <= 47; i++) {
                addWhiteBag();
            }
        }
    }

    void setRPM(Mdouble rpm) {
        rpm_ = rpm;
        for (auto& w : wallHandler)
        {
            if (w->getAngularVelocity().Z!=0)
                w->setAngularVelocity({0,0,rpm_*pi/30.});
        }
    }

    //remove particles according to (max.) outflow rate
    void actionsAfterTimeStep() override
    {
        if (batchFill_) {
            checkIfBagShouldBeAdded();
        }

        // Functions after this statement only get executed every 100th timestep (if counter==100)
        static unsigned counter = 0;
        if (++counter!=100) return;
        else counter = 0;

        if (delayedStartOfRotation_ == true && getKineticEnergy() < 1e-2 * getElasticEnergy())
        {
            logger(INFO,"\nStart Rotation at time %", getTime());
            delayedStartOfRotation_ = false;
            for (auto& w : wallHandler)
            {
                if (w->getAngularVelocity().Z!=0)
                    w->setAngularVelocity({0,0,rpm_*pi/30.});
            }
        }

        // outflowVolume_: volume of particles already deleted
        // outflowRate_*getTime() : volume of outflow requested
        // check if the volume of particles deleted is smaller than the requested outflow
        if (outflowVolume_<outflowRate_*getTime())
        {
            for (auto p : particleHandler) {
                if (p->getPosition().Z<0.5*getZMin()) {
                    outflowVolume_ += p->getVolume();
                    particleHandler.removeObject(p->getIndex());
                    logger(INFO,"Deleted particle because % < %",outflowVolume_,outflowRate_*getTime());
                    return;
                }
            }
        }
    }

    /**
     * Adds a second particle species; assumes one existing particle and wall species.
     * Comment by Ant. Actually this functions does a lot more than it claims as it puts some of the particles in this species based on initial spactial distrinution.
     * @param particle2Species          the new particle species
     * @param particle2particle1Species species describing the contact between particles of different species
     * @param particle2WallSpecies      species describing the contact between particles of the new species and a wall
     * @param conversionCondition       condition for which particles should be converted to the second species
     *                                  (e.g. [] (const BaseParticle* p const) {return p->getPosition().X<0;} )
     */
    void addSecondParticleSpecies(ParticleSpecies* particle2Species,
                                  BaseSpecies* particle2Particle1Species,
                                  BaseSpecies* particle2WallSpecies,
                                  std::function<bool(const BaseParticle* const)> conversionCondition)
    {
        if (speciesHandler.getNumberOfObjects()!=2)
            logger(ERROR,"Cannot add second particle species, as this assumes one existing particle and wall species");
        //add species
        auto p2Species = speciesHandler.copyAndAddObject(particle2Species); //use this pointer to set particles to the new species
        //set the mixed species
        speciesHandler.getMixedObject(0,2)->mixAll(particle2Particle1Species,particle2Particle1Species);
        speciesHandler.getMixedObject(1,2)->mixAll(particle2WallSpecies,particle2WallSpecies);

        //convert certain particles to the new species
        unsigned counter = 0;
        for (auto p : particleHandler) {
            if (conversionCondition(p)) {
                p->setSpecies(p2Species);
                counter++;
            }
        }
        logger(INFO,"Converted % particles to the new species",counter);
    }

    /**
     * This is the function which checks if it is time to add another bag of material to the mixing; if it is in the filling stage.
     */
    void checkIfBagShouldBeAdded()
    {
        if (bagsPerRevolution_*getTime()*(rpm_/60.0) >= (bagsEmptied_))
            //true if is time to empty a new bag into the mixer (one bag per revolution)
        {
            bagsEmptied_++;
            if (bagsEmptied_<=5)
            {
                addWhiteBag();
            }
            else if (bagsEmptied_<=11)
            {
                addOrangeBag();
            }
            else if (bagsEmptied_<=60)
            {
                addWhiteBag();
            }
            else
            {
                logger(INFO,"Adding no bag at time %",getTime());
            }
        }

    }

    /**
     * This function adds a 10kg bag of orange particles to the mixer
     */
    void addOrangeBag()
    {
        //(volume,species,{fraction of each particles size FOR CURRENT SPECIES})
        //volume is in litres -- need to convert from density
        // Sizes: {< 2mm, 2-4 mm, 4-6 mm, 6-8 mm, >8 mm}
        // ENSURE THAT PARTICLE FRACTION ADDS UP TO UNITY!
        addParticlesToMixer(10,0,{0.5,2,6,1.5,0});
    }

    /**
     * This function adds a 10kg bag of white particles to the mixer
     */
    void addWhiteBag()
    {
        addParticlesToMixer(10,2,{7.2,54,27,1.8,0});
        //addParticlesToMixer(10.1,2,{0,0,1,0,0});
    }

    /**
     * This functions insert particles in the mixer
     * Inserts particles of different size into the domain, at given percentages:
     * \param[in] fillMass               total mass of  particles that are inserted
     * \param[in] SpeciesInd         species of the particles that are inserted
     * \param[in] percentageInBin[0] volume-percentage of particles with size 0.25<radius/particleRadius_<0.5
     * \param[in] percentageInBin[1] volume-percentage of particles with size 0.5<radius/particleRadius_<1.0
     * \param[in] percentageInBin[2] volume-percentage of particles with size 1.0<radius/particleRadius_<1.5
     * \param[in] percentageInBin[3] volume-percentage of particles with size 1.5<radius/particleRadius_<2.0
     * \param[in] percentageInBin[4] volume-percentage of particles with size 2.0<radius/particleRadius_<4.0
     *
     * So for full-scale simulations, use particleRadius_=4mm
     */
    void addParticlesToMixer(Mdouble fillMass, int SpeciesInd, std::vector<Mdouble> percentageInBin)
    {
        BaseParticle p;
        p.setSpecies(speciesHandler.getObject(SpeciesInd));

        //check if percentages add up to pne, otherwise normalise
        Mdouble totalPercentage = 0;
        for (const Mdouble p : percentageInBin) {
            totalPercentage += p;
        }
        if (totalPercentage!=1.0) {
            //logger(WARN, "Percentages add up to %, normalizing.", totalPercentage);
            for (Mdouble& p : percentageInBin) {
                p /= totalPercentage;
            }
        }

        //compute fillVolume
        double fillVolume = fillMass/p.getSpecies()->getDensity();

        //Setup size ranges
        std::vector<double> sizeRanges(6);
        //originally implemented distribution
/*	    sizeRanges[0]=0.25*particleRadius_;
            sizeRanges[1]=0.5*particleRadius_;
            sizeRanges[2]=1.0*particleRadius_;
            sizeRanges[3]=1.5*particleRadius_;
            sizeRanges[4]=2.0*particleRadius_;
            sizeRanges[5]=4.0*particleRadius_;
*/
        //narrowed distribution for faster simulations (hopefully)
        sizeRanges[0]=0.4*particleRadius_;
        sizeRanges[1]=0.5*particleRadius_;
        sizeRanges[2]=1.0*particleRadius_;
        sizeRanges[3]=1.5*particleRadius_;
        sizeRanges[4]=2.0*particleRadius_;
        sizeRanges[5]=2.05*particleRadius_;

        // setup mean volume of particles in each bin
        // barV = int_r0^r1 4/3*pi*r^3 dr /(r1-r0) = pi/3 * (r1^4-r0^4)/(r1-r0)
        std::vector<double> meanVolume(5);
        for (int i =0; i<5; i++) {
            meanVolume[i] = pi/3.0 * (pow(sizeRanges[i+1],4)-pow(sizeRanges[i],4))/(sizeRanges[i+1]-sizeRanges[i]);
        }

        // compute mean radius
        Mdouble totalMeanVolume = 0;
        for (int i =0; i<5; i++) {
            totalMeanVolume += percentageInBin[i]*meanVolume[i];
        }
        //V=4/3*pi*r^3 -> r=cbrt(3/4*V/pi)
        Mdouble meanRadius = cbrt(0.75*totalMeanVolume/pi);

        // compute number of particles in each bin
        std::vector<unsigned> particlesPerBin(5);
        for (int i =0; i<5; i++) {
            particlesPerBin[i] = lround(percentageInBin[i]*fillVolume/meanVolume[i]);
        }

        unsigned numberOfParticles = 0;
        for (int i =0; i<5; i++) {
            numberOfParticles += particlesPerBin[i];
        }
        logger(INFO,"Creating % particles of mean radius %mm. Number per bin: %, %, %, %, %",
               numberOfParticles,meanRadius*1000, particlesPerBin[0],particlesPerBin[1],
               particlesPerBin[2],particlesPerBin[3],particlesPerBin[4]);

        unsigned numberOfParticlesToInsert = numberOfParticles;
        while (numberOfParticlesToInsert>0)
        {
            int rand=random.getRandomNumber(0,numberOfParticlesToInsert);

            int pickedBin=0;
            while (rand>=particlesPerBin[pickedBin])
            {
                rand -= particlesPerBin[pickedBin];
                pickedBin++;
            }

            if (particlesPerBin[pickedBin]==0)
                logger(ERROR,"Error in Particle insertion");

            p.setRadius(random.getRandomNumber(sizeRanges[pickedBin],sizeRanges[pickedBin+1]));

            Mdouble fillHeight = getZMax();
            unsigned count = 0;
            do
            {
                count++;

                Mdouble x = random.getRandomNumber(getXMin(), getXMax());
                Mdouble y = random.getRandomNumber(getYMin(), getYMax());
                Mdouble z = random.getRandomNumber(getZMax(), fillHeight);
                p.setPosition({x, y, z});

                fillHeight+=p.getRadius();

                if (count==100)
                {
                    logger(INFO, "Trying to insert particle with radius % position % (bin % tries % remaining %)",
                           p.getRadius(), p.getPosition(), pickedBin, count, numberOfParticlesToInsert);
                }
            } while (checkParticleForInteraction(p)==false);

            particleHandler.copyAndAddObject(p);
            particlesPerBin[pickedBin]--;
            numberOfParticlesToInsert--;

            //logger(INFO,"Inserted particle with radius % (bin % tries % remaining %)",p.getRadius(),pickedBin,count,numberOfParticlesToInsert);
        }
        //logger(INFO,"Insertion complete");
        logger(INFO,"Added bag % at time % #particles %, total mass % kg",bagsEmptied_,getTime(),particleHandler.getNumberOfObjects(),particleHandler.getMass());

    }



private:
    const Mdouble particleRadius_;
    const Mdouble outflowRate_;
    Mdouble outflowVolume_;
    bool delayedStartOfRotation_;
    bool batchFill_;
    Mdouble rpm_;
    //This is the number of bags that have been emptied into the mixer.
    int bagsEmptied_;
    Mdouble bagsPerRevolution_;
};

#endif //MERCURY_CSK_H
