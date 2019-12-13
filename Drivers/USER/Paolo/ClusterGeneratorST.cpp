//
// Created by paolo on 29-6-19.
//
#include "ClusterGenerator.h"

/*!
 * ClusterGeneratorSelfTest
 * In this file three different clusters are created:
 *      the first one is created imposing a size dispersity greater than zero, setConstantRestitution(true), and imposing collision time
 *  and restitution coefficient. setConstantRestitution is true because, having a size dispersity, this is the only way to define
 *  all parameters correctly (the mass of the  actual particle that will be created are random, and so unknown). For this reasons,
 *  in order to have a good value of computational accuracy collision time has to be computed taking into account the smallest mass
 *  possible (that is the mass of a particle having radius R = radiusParticle * (1 - sizeDispersityParticle) );
 *      the second one is created imposing a size dispersity greater than zero, setConstantRestitution(true) and imposing
 *  loading stiffness and restitution coefficient. In this case loading stiffness has to be set as Kp / m, where Kp is the desired value of loading
 *  stiffness and m is (just like in the first cluster) the smallest mass possible, so the mass of a particle having radius
 *  R = radiusParticle * (1 - sizeDispersityParticle);
 *      the third one is set with monodispersity (in order to do that is enough not setting size dispersity, given that its default value is 0.0),
 *  with setConstantRestitution(false) (also default) and imposing loading stiffness and restitution coefficient. In this case loading stiffness
 *  is the desired value of loading stiffness, but it is necessary to specify the mass (which in the previous cases was 1 because of constant
 *  restitution). This time massParticle is the mass of a particle having radius = radiusParticle. Accuracy in this case is guaranteed because
 *  collision time will be computed over the smallest particle's mass.
 *  In the three clusters the random seed is fixed;
 *  In all three cases it is fundamental to define radiusParticle and species, otherwise the program will give an error.
 *
 * @return
 */

int main(){

    Mdouble radiusParticle = 1e-6;
    Mdouble sizeDispersityParticle = 0.1;
    Mdouble densityParticle = 2500;
    Mdouble restitutionCoefficient = 0.5;
    Mdouble penetrationDepthMax = 0.05;
    Mdouble slidingFrictionCoefficient = 0.5;
    Mdouble slidingRollingCoefficient = 0.3;
    Mdouble slidingTorsionCoefficient = 0.0;
    Mdouble loadingStiffness = 1.0e5;
    Mdouble unLoadingStiffnessMax = 5.0e5;
    Mdouble cohesionStiffness = 0.5e5;
    //! smallest mass possible (that is the mass of a particle having radius R = radiusParticle * (1 - sizeDispersityParticle) )
    Mdouble smallestMass = densityParticle * 4 * constants::pi * pow(radiusParticle*(1-sizeDispersityParticle), 3) / 3;
    //! collision time computed over the smallest mass possible
    Mdouble collisionTimeSmallestMass = sqrt( smallestMass * ( pow(constants::pi, 2) + pow(log(restitutionCoefficient), 2) ) / ( 2 * loadingStiffness ) );



    //CREATION OF THE FIRST CLUSTER

    LinearPlasticViscoelasticFrictionSpecies* species1;

    species1 = new LinearPlasticViscoelasticFrictionSpecies;
    species1 -> setConstantRestitution(true);
    species1 -> setDensity(densityParticle);
    // mass is set to 1 (any value would be good because of setConstantRestitution(true) ).
    species1 -> setCollisionTimeAndRestitutionCoefficient(collisionTimeSmallestMass, restitutionCoefficient, 1);
    species1 -> setUnloadingStiffnessMax(species1 -> getLoadingStiffness() * unLoadingStiffnessMax / loadingStiffness);
    species1 -> setCohesionStiffness(species1 -> getLoadingStiffness() * cohesionStiffness / loadingStiffness);
    species1 -> setPenetrationDepthMax(penetrationDepthMax);

    species1 -> setSlidingFrictionCoefficient(slidingFrictionCoefficient);
    species1 -> setSlidingStiffness(species1 -> getLoadingStiffness()*2.0/7.0);
    species1 -> setSlidingDissipation(species1 -> getDissipation()*2.0/7.0);
    species1 -> setRollingFrictionCoefficient(slidingRollingCoefficient);
    species1 -> setRollingStiffness(species1 -> getLoadingStiffness()*2.0/7.0);
    species1 -> setRollingDissipation(species1 -> getDissipation()*2.0/7.0);
    species1 -> setTorsionFrictionCoefficient(slidingTorsionCoefficient);
    species1 -> setTorsionStiffness(species1 -> getLoadingStiffness()*2.0/7.0);
    species1 -> setTorsionDissipation(species1 -> getDissipation()*2.0/7.0);



    ClusterGenerator a;

    a.clusterProperties.setRadiusParticle(radiusParticle);
    a.clusterProperties.setNumberOfParticles(3);
    a.clusterProperties.setClusterId(1);
    a.clusterProperties.setSizeDispersityParticle(sizeDispersityParticle);
    a.clusterProperties.setParticleSpecies(species1);
    a.create();
    helpers::check(a.clusterProperties.getAverageOverlap(), penetrationDepthMax, penetrationDepthMax / 100, "Average overlap");
    std::cout << std::endl << std::endl;


    //CREATION OF THE SECOND CLUSTER

    LinearPlasticViscoelasticFrictionSpecies* species2;

    species2 = new LinearPlasticViscoelasticFrictionSpecies;
    species2 -> setConstantRestitution(true);
    species2 -> setDensity(densityParticle);
    // mass is set to 1 (any value would be good because of setConstantRestitution(true) ).
    species2 -> setStiffnessAndRestitutionCoefficient(loadingStiffness/smallestMass, restitutionCoefficient,1);
    species2 -> setUnloadingStiffnessMax(species2 -> getLoadingStiffness() * unLoadingStiffnessMax / loadingStiffness);
    species2 -> setCohesionStiffness(species2 -> getLoadingStiffness() * cohesionStiffness / loadingStiffness);
    species2 -> setPenetrationDepthMax(penetrationDepthMax);

    species2 -> setSlidingFrictionCoefficient(slidingFrictionCoefficient);
    species2 -> setSlidingStiffness(species2 -> getLoadingStiffness()*2.0/7.0);
    species2 -> setSlidingDissipation(species2 -> getDissipation()*2.0/7.0);
    species2 -> setRollingFrictionCoefficient(slidingRollingCoefficient);
    species2 -> setRollingStiffness(species2 -> getLoadingStiffness()*2.0/7.0);
    species2 -> setRollingDissipation(species2 -> getDissipation()*2.0/7.0);
    species2 -> setTorsionFrictionCoefficient(slidingTorsionCoefficient);
    species2 -> setTorsionStiffness(species2 -> getLoadingStiffness()*2.0/7.0);
    species2 -> setTorsionDissipation(species2 -> getDissipation()*2.0/7.0);



    ClusterGenerator b;

    b.clusterProperties.setRadiusParticle(radiusParticle);
    b.clusterProperties.setNumberOfParticles(3);
    b.clusterProperties.setClusterId(2);
    b.clusterProperties.setSizeDispersityParticle(sizeDispersityParticle);
    b.clusterProperties.setParticleSpecies(species2);
    b.create();
    helpers::check(b.clusterProperties.getAverageOverlap(), penetrationDepthMax, penetrationDepthMax / 100, "Average overlap");
    std::cout << std::endl << std::endl;



    //CREATION OF THE THIRD CLUSTER

    //! mass of a particle having radius = radiusParticle.
    Mdouble massParticle = densityParticle * 4 * constants::pi * pow(radiusParticle, 3) / 3;

    LinearPlasticViscoelasticFrictionSpecies* species3;

    species3 = new LinearPlasticViscoelasticFrictionSpecies;
    species3 -> setDensity(densityParticle);
    // Here mass particle has to be set correctly.
    species3 -> setStiffnessAndRestitutionCoefficient(loadingStiffness, restitutionCoefficient, massParticle);
    species3 -> setUnloadingStiffnessMax(species3 -> getLoadingStiffness() * unLoadingStiffnessMax / loadingStiffness);
    species3 -> setCohesionStiffness(species3 -> getLoadingStiffness() * cohesionStiffness / loadingStiffness);
    species3 -> setPenetrationDepthMax(penetrationDepthMax);

    species3 -> setSlidingFrictionCoefficient(slidingFrictionCoefficient);
    species3 -> setSlidingStiffness(species3 -> getLoadingStiffness()*2.0/7.0);
    species3 -> setSlidingDissipation(species3 -> getDissipation()*2.0/7.0);
    species3 -> setRollingFrictionCoefficient(slidingRollingCoefficient);
    species3 -> setRollingStiffness(species3 -> getLoadingStiffness()*2.0/7.0);
    species3 -> setRollingDissipation(species3 -> getDissipation()*2.0/7.0);
    species3 -> setTorsionFrictionCoefficient(slidingTorsionCoefficient);
    species3 -> setTorsionStiffness(species3 -> getLoadingStiffness()*2.0/7.0);
    species3 -> setTorsionDissipation(species3 -> getDissipation()*2.0/7.0);



    ClusterGenerator c;

    c.clusterProperties.setRadiusParticle(radiusParticle);
    c.clusterProperties.setNumberOfParticles(3);
    c.clusterProperties.setClusterId(3);
    c.clusterProperties.setParticleSpecies(species3);
    c.create();
    helpers::check(c.clusterProperties.getAverageOverlap(), penetrationDepthMax, penetrationDepthMax / 100, "Average overlap");




    return 0;
}
