//Copyright (c) 2013-2020, The MercuryDPM Developers Team. All rights reserved.
//For the list of developers, see <http://www.MercuryDPM.org/Team>.
#include "HorizontalMixer.h"

//adds more details to the walls (double thick, end at z=0 and z=1.4, circular base plate)
class HorizontalMixerWalls : public HorizontalMixer {

    public:

    HorizontalMixerWalls(Mdouble particleRadius, Mdouble rotationSpeed, Mdouble timeMin, Mdouble fillHeight)
            : HorizontalMixer(particleRadius,rotationSpeed, timeMin, fillHeight)
    {}

    ///sets four walls, leftScrewCore, rightScrewCore, leftScrewBottom, rightScrewBottom
    void setScrewCore(const ParticleSpecies* s, Mdouble screwCenter, Mdouble screwBaseHeight,
                      Mdouble coreTopHeight, Mdouble coreTopRadius,
                      Mdouble coreBottomHeight, Mdouble coreBottomRadius) override
    {
        //the inner, thinner core cylinder in the screw
        auto screwCore = wallHandler.copyAndAddObject
                (AxisymmetricIntersectionOfWalls(Vec3D(screwCenter,0,0),Vec3D(0,0,1),
                                                 {{Vec3D(-.3,0,-1),Vec3D(coreTopRadius,0,coreTopHeight)},//slightly slant top, so particles roll off
                                                  {Vec3D(-1,0,0),Vec3D(coreTopRadius,0,0)},{Vec3D(0,0,1),Vec3D(0,0,0)}},s));

        //the bottom, thicker core cylinder in the screw
        auto screwBottom = wallHandler.copyAndAddObject
                (AxisymmetricIntersectionOfWalls(Vec3D(screwCenter,0,0),Vec3D(0,0,1),
                                                 {{Vec3D(0,0,-1),Vec3D(0,0,coreBottomHeight)},
                                                  {Vec3D(-1,0,0),Vec3D(coreBottomRadius,0,0)},{Vec3D(0,0,1),Vec3D(0,0,0)}},s));
    }

    ///sets other walls that define the outer boundary
    void setOuterWalls(const ParticleSpecies* s, Mdouble outerBaseRadius, Mdouble screwCenter, Mdouble outerTopCenter,
                       Mdouble outerTopRadius, Mdouble outerTopHeight) override
    {
        //cylindrical wall
        Vec3D normal = (Vec3D(outerTopCenter,0,outerTopHeight)-Vec3D(screwCenter,0,0));
        normal.normalise();
        Vec3D position = Vec3D(screwCenter,0,0);
        Vec3D normalWall = Vec3D(outerTopHeight,0,outerBaseRadius-outerTopRadius);
        normalWall.normalise();
        Vec3D positionWall = Vec3D(outerBaseRadius,0,0);
        AxisymmetricIntersectionOfWalls outerWall(position, normal, {{normalWall,positionWall},{Vec3D(0,0,-1),Vec3D(0,0,1.4)},{Vec3D(0,0,1),Vec3D(0,0,0)},{-normalWall,positionWall+Vec3D(.02,0,0)}},s);
        wallHandler.copyAndAddObject(outerWall);

        //bottom plate
        AxisymmetricIntersectionOfWalls bottomPlate(position, normal, {{Vec3D(0,0,-1),Vec3D(0,0,0)},{Vec3D(0,0,1),Vec3D(0,0,-0.02)},{-normalWall,positionWall}},s);
        wallHandler.copyAndAddObject(bottomPlate);
    }

    void setupInitialConditions() override {
        HorizontalMixer::setupInitialConditions();
        particleHandler.clear();
    }

};

int main()
{
    Mdouble revolutionsPerSecond = 0.25;
    Mdouble rotationSpeed = 0.25*constants::pi*revolutionsPerSecond;
    Mdouble particleRadius = 0.2;
    Mdouble timeMin = 0.0;
    Mdouble fillHeight = 0;
    HorizontalMixerWalls mixer(particleRadius, rotationSpeed, timeMin, fillHeight);

    //set name
    mixer.setName("HorizontalMixerWalls");

    //remove old files (note this is a bit dangerous)
    //mixer.removeOldFiles();

    //set species
    auto s = mixer.speciesHandler.copyAndAddObject(LinearViscoelasticSlidingFrictionSpecies());
    s->setDensity(2000);
    Mdouble mass = s->getMassFromRadius(particleRadius);
    s->setCollisionTimeAndRestitutionCoefficient(0.01, 0.5, mass);
    s->setSlidingFrictionCoefficient(0.5);
    s->setSlidingStiffness(2.0/7.0*s->getStiffness());
    s->setSlidingDissipation(2.0/7.0*s->getDissipation());

    //set timestep
    mixer.setTimeStep(0.2 * s->getCollisionTime(mass));
    //save every 2 collision times for smooth viewable output
    mixer.setSaveCount((unsigned)(20.0*s->getCollisionTime(mass)/mixer.getTimeStep()));
    std::cout << "Savecount: " << mixer.dataFile.getSaveCount() << std::endl;

    mixer.setTimeMax(100);

    mixer.solve();
    return 0;
}
