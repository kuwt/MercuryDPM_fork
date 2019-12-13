//
// Created by paolo on 29-6-19.
//
//
// Created by paolo on 29-6-19.
//
#include "ClusterGenerator.h"

/*!
 * ClusterGeneratorRestartSelfTest
 *  In this file restart process for a cluster is testd.
 *  In orde to do so the progress is prematurely stopped thanks to doRestartSelfTest(true) and than restarted imposing
 *  setRestarted(true). ClusterID is 0 (default value).
 *  The random seed is fixed;
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



    LinearPlasticViscoelasticFrictionSpecies* species;

    species = new LinearPlasticViscoelasticFrictionSpecies;
    species -> setConstantRestitution(true);
    species -> setDensity(densityParticle);
    // mass is set to 1 (any value would be good because of setConstantRestitution(true) ).
    species -> setCollisionTimeAndRestitutionCoefficient(collisionTimeSmallestMass, restitutionCoefficient, 1);
    species -> setUnloadingStiffnessMax(species -> getLoadingStiffness() * unLoadingStiffnessMax / loadingStiffness);
    species -> setCohesionStiffness(species -> getLoadingStiffness() * cohesionStiffness / loadingStiffness);
    species -> setPenetrationDepthMax(penetrationDepthMax);

    species -> setSlidingFrictionCoefficient(slidingFrictionCoefficient);
    species -> setSlidingStiffness(species -> getLoadingStiffness()*2.0/7.0);
    species -> setSlidingDissipation(species -> getDissipation()*2.0/7.0);
    species -> setRollingFrictionCoefficient(slidingRollingCoefficient);
    species -> setRollingStiffness(species -> getLoadingStiffness()*2.0/7.0);
    species -> setRollingDissipation(species -> getDissipation()*2.0/7.0);
    species -> setTorsionFrictionCoefficient(slidingTorsionCoefficient);
    species -> setTorsionStiffness(species -> getLoadingStiffness()*2.0/7.0);
    species -> setTorsionDissipation(species -> getDissipation()*2.0/7.0);



    ClusterGenerator a;

    a.clusterProperties.setRadiusParticle(radiusParticle);
    a.clusterProperties.setNumberOfParticles(3);
    a.clusterProperties.setSizeDispersityParticle(sizeDispersityParticle);
    a.clusterProperties.setParticleSpecies(species);
    //! quin prima c'era a.clusterProperties.isRestartSelfTest(true) che però è stata eliminata.
    a.create();
    

    ClusterGenerator b;
    b.clusterProperties.setRestarted(true);
    b.clusterProperties.setName("Cluster_ID_0"); // Very important in order to restart the calculation correctly.
    b.create();
    helpers::check(b.clusterProperties.getAverageOverlap(), penetrationDepthMax, penetrationDepthMax * 1e-4, "Average overlap");


    return 0;
}
