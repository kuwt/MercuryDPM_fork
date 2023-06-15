//
// Created by mitchel on 1/23/20.
//

#include "DieFilling.h"

void DieFilling::setupInitialConditions()
{
    //logger(ERROR,"MercuryDieFilling::setBC() not overwritten correctly");
    logger(INFO,"Call to MercuryDieFilling::setupInitialConditions()");
    
    double tc = 50.* getTimeStep();
    
    LinearPlasticViscoelasticFrictionSpecies steel;
    steel.setHandler(&speciesHandler);
    steel.setDensity(rhoSteel);
    Mdouble mSteel = steel.getMassFromRadius(dSteel/2.0);
    steel.setCollisionTimeAndRestitutionCoefficient(tc,rSteel,mSteel);
    steel.setUnloadingStiffnessMax(5.0*steel.getLoadingStiffness());
    steel.setPenetrationDepthMax(0.05*dSteel);
    steel.setSlidingFrictionCoefficient(slidingFrictionSteel);
    steel.setSlidingStiffness(2.0 / 7.0 * steel.getLoadingStiffness());
    steel.setSlidingDissipation(2.0 / 7.0 * steel.getDissipation());
    steel.setRollingFrictionCoefficient(rollingFrictionSteel);
    steel.setRollingStiffness(0.2 * steel.getLoadingStiffness());
    steel.setRollingDissipation(0.2 * steel.getDissipation());
    pSteel = speciesHandler.copyAndAddObject(steel);
    
    LinearPlasticViscoelasticFrictionSpecies glass;
    glass.setHandler(&speciesHandler);
    glass.setDensity(rhoGlass);
    Mdouble mGlass = glass.getMassFromRadius(dGlass/2.0);
    glass.setCollisionTimeAndRestitutionCoefficient(tc,rGlass,mGlass);
    glass.setUnloadingStiffnessMax(5.0*glass.getLoadingStiffness());
    glass.setPenetrationDepthMax(0.05*dGlass);
    glass.setSlidingFrictionCoefficient(slidingFrictionGlass);
    glass.setSlidingStiffness(2.0 / 7.0 * glass.getLoadingStiffness());
    glass.setSlidingDissipation(2.0 / 7.0 * glass.getDissipation());
    glass.setRollingFrictionCoefficient(rollingFrictionGlass);
    glass.setRollingStiffness(0.2 * glass.getLoadingStiffness());
    glass.setRollingDissipation(0.2 * glass.getDissipation());
    pGlass = speciesHandler.copyAndAddObject(glass);
    
    LinearPlasticViscoelasticFrictionSpecies MCC;
    MCC.setHandler(&speciesHandler);
    MCC.setDensity(rhoMCC);
    Mdouble mMCC = MCC.getMassFromRadius(dMCC/2.0);
    MCC.setCollisionTimeAndRestitutionCoefficient(tc,rMCC,mMCC);
    MCC.setUnloadingStiffnessMax(5.0*MCC.getLoadingStiffness());
    MCC.setPenetrationDepthMax(0.05*dMCC);
    MCC.setSlidingFrictionCoefficient(slidingFrictionMCC);
    MCC.setSlidingStiffness(2.0 / 7.0 * MCC.getLoadingStiffness());
    MCC.setSlidingDissipation(2.0 / 7.0 * MCC.getDissipation());
    MCC.setRollingFrictionCoefficient(rollingFrictionMCC);
    MCC.setRollingStiffness(0.2 * MCC.getLoadingStiffness());
    MCC.setRollingDissipation(0.2 * MCC.getDissipation());
    pMCC = speciesHandler.copyAndAddObject(MCC);
    
    
    SphericalParticle p0;
    p0.setSpecies(pSteel);
    p0.setRadius(dSteel/2.);
    
    double percXMin = 0.0;
    double percXMax = 1.0;
    double percYMin = 0.0;
    double percYMax = 1.0;
    double percZMin = 0.65;
    double percZMax = 0.95;

    const unsigned int numToInsert = 1;
    for (unsigned int i=0; i<numToInsert; i++)
    {
        Vec3D getRandomPos;
        Vec3D getRandomVel;
        int failcounter = 0;

        do
        {
            getRandomPos.X = random.getRandomNumber(getXMin() + p0.getRadius() + percXMin*(getXMax()-getXMin()), percXMax*getXMax() - p0.getRadius());
            getRandomPos.Y = random.getRandomNumber(getYMin() + p0.getRadius() + percYMin*(getYMax()-getYMin()), percYMax*getYMax() - p0.getRadius());
            getRandomPos.Z = random.getRandomNumber(getZMin() + p0.getRadius() + percZMin*(getZMax()-getZMin()), percZMax*getZMax() - p0.getRadius());

            getRandomPos = {0.5*getXMax(), 0.5*getYMax(), 0.9*getZMax()};

            if (i==1)
            {
            }
            //getRandomPos.X = 10;

            p0.setPosition(getRandomPos);
            p0.setVelocity(getRandomVel);
            failcounter++;
            if (failcounter == 1000) {logger(INFO,"Failcounter reached"); break;}

        } while (!checkParticleForInteraction(p0));

        if (failcounter != 1000)
        {
            particleHandler.copyAndAddObject(p0);
        }
        else {break;}
    }
    
    //std::cout << "Placed " << particleHandler.getNumberOfObjects() << " out of " << numToInsert <<  " particles" << std::endl;
    
    InfiniteWall wall;
    wall.setSpecies(pSteel);
    wall.set(Vec3D(-1,0,0),{getXMin(),0.0,0.0});
    wallHandler.copyAndAddObject(wall);
    wall.set(Vec3D( 1,0,0),{getXMax(),0.0,0.0});
    wallHandler.copyAndAddObject(wall);
    wall.set(Vec3D(0,-1,0),{0.0,getYMin(),0.0});
    wallHandler.copyAndAddObject(wall);
    wall.set(Vec3D(0, 1,0),{0.0,getYMax(),0.0});
    wallHandler.copyAndAddObject(wall);
    wall.set(Vec3D(0, 0,-1),{0.0,0.0,getZMin()});
    wallHandler.copyAndAddObject(wall);

//    wall.set(Vec3D(0, 0, 1),getMax());
//    wallHandler.copyAndAddObject(wall);
    
    std::cout<< "Completed call to MercuryDieFilling::setupInitialConditions" << std::endl;
}


void DieFilling::pinBC()
{
    logger(DEBUG,"Call to pinBC()");
    
    // Boundaries are numbered:
    // 0 is at the bottom
    // 1 2 3 4 from the front  proceeding anticlockwise
    // 5 is at the top
    
    // ---------------------------------------------------------------------------------------------- //
    // FOR UNIFORM INFLOW in z-ir
    //
    //   ^   ^   ^  ^   ^
    //   |   |   |  |   |
    //
    for (unsigned iBound = 0; iBound < mesh_pt()->nboundary(); iBound++)
    {
        if (iBound == 0)
        {
            unsigned long int nNode = mesh_pt()->nboundary_node(iBound);
            for (unsigned iNode = 0; iNode < nNode; iNode++)
            {
                //Pin u and v and w
                mesh_pt()->boundary_node_pt(iBound, iNode)->pin(0);
                mesh_pt()->boundary_node_pt(iBound, iNode)->pin(1);
                mesh_pt()->boundary_node_pt(iBound, iNode)->pin(2);
            }
        }
        else if (iBound == 5)
        {
            unsigned long int nNode = mesh_pt()->nboundary_node(iBound);
            for (unsigned iNode = 0; iNode < nNode; iNode++)
            {
                //Pin u and v and w
                mesh_pt()->boundary_node_pt(iBound, iNode)->pin(0);
                mesh_pt()->boundary_node_pt(iBound, iNode)->pin(1);
            }
        }
        else
        {
            unsigned long int nNode = mesh_pt()->nboundary_node(iBound);
            for (unsigned iNode = 0; iNode < nNode; iNode++)
            {
                //Pin all velocities
                for (unsigned i = 0; i < 3; i++)
                {
                    mesh_pt()->boundary_node_pt(iBound, iNode)->pin(i);
                }
            }
        }
    }
}


void DieFilling::setBC()
{
    //logger(ERROR,"UnderResolvedCoupling::setBC() not overwritten correctly");
    
    // Pin redudant pressure dofs
    oomph::RefineableAndersonJacksonEquations<3>::pin_redundant_nodal_pressures(mesh_pt()->element_pt());
    //Fix 3-th pressure value in first element to 0.0.
    dynamic_cast<oomph::RefineableAJQCrouzeixRaviartElement<3> *>(mesh_pt()->element_pt(0))->fix_pressure(1,0.0);
    
    // Boundaries are numbered:
    // 0 is at the bottom
    // 1 2 3 4 from the front  proceeding anticlockwise
    // 5 is at the top
    
    // Periodic flow
    //double zVelIn = 0.05+ 0.3*sin(9.*getTime());
    
    // Uniform inflow
    double zVelIn = 0.0;// getInflowVel(getTime());// getInflowVel(getTime());
    
    for (unsigned iBound = 0; iBound < mesh_pt()->nboundary(); iBound++)
    {
        if (iBound == 0)
        {
            unsigned long int nNode = mesh_pt()->nboundary_node(iBound);
            for (unsigned iNode = 0; iNode < nNode; iNode++)
            {
                mesh_pt()->boundary_node_pt(iBound, iNode)->set_value(0,0.0);
                mesh_pt()->boundary_node_pt(iBound, iNode)->set_value(1,0.0);
                mesh_pt()->boundary_node_pt(iBound, iNode)->set_value(2,zVelIn);
            }
        }
        else if (iBound == 5)
        {
            unsigned long int nNode = mesh_pt()->nboundary_node(iBound);
            for (unsigned iNode = 0; iNode < nNode; iNode++)
            {
                mesh_pt()->boundary_node_pt(iBound, iNode)->set_value(0,0.0);
                mesh_pt()->boundary_node_pt(iBound, iNode)->set_value(1,0.0);
            }
        }
        else if (iBound == 1 || iBound == 3)
        {
            unsigned long int nNode = mesh_pt()->nboundary_node(iBound);
            for (unsigned iNode = 0; iNode < nNode; iNode++)
            {
                mesh_pt()->boundary_node_pt(iBound, iNode)->set_value(0,0.0);
                mesh_pt()->boundary_node_pt(iBound, iNode)->set_value(1,0.0);
                mesh_pt()->boundary_node_pt(iBound, iNode)->set_value(2,0.0);
            }
        }
        else if (iBound == 2)
        {
            unsigned long int nNode = mesh_pt()->nboundary_node(iBound);
            for (unsigned iNode = 0; iNode < nNode; iNode++) {
                mesh_pt()->boundary_node_pt(iBound, iNode)->set_value(0,0.0);
                mesh_pt()->boundary_node_pt(iBound, iNode)->set_value(1,0.0);
                mesh_pt()->boundary_node_pt(iBound, iNode)->set_value(2,0.0);
            }
        }
        else if (iBound == 4)
        {
            unsigned long int nNode = mesh_pt()->nboundary_node(iBound);
            for (unsigned iNode = 0; iNode < nNode; iNode++) {
                mesh_pt()->boundary_node_pt(iBound, iNode)->set_value(0,0.0);
                mesh_pt()->boundary_node_pt(iBound, iNode)->set_value(1,0.0);
                mesh_pt()->boundary_node_pt(iBound, iNode)->set_value(2,0.0);
            }
        }
    }
}





