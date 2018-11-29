#ifndef MERCURY_ADDSPECIES_H
#define MERCURY_ADDSPECIES_H

#include<vector>
#include <Species/LinearViscoelasticFrictionReversibleAdhesiveSpecies.h>
#include "Species/LinearPlasticViscoelasticFrictionSpecies.h"
#include "Math/PSD.h"
#include "Math/ExtendedMath.h"
#include "DPMBase.h"

using constants::pi;
using mathsFunc::cubic;

enum class SpeciesType : unsigned {
    Glass,
    Cocoa,
    RubberWall,
    SteelWall,
    MPT, //Metroprolol
    APAPP, //APAP powder
    PH101, //MCC, Avicel
    SD100, //Mannitol, Pearlitol
};

/*!
 * write to file
 */
std::ostream& operator<<(std::ostream& os, SpeciesType type)
{
    switch (type) {
        case SpeciesType::Glass: os << "Glass"; break;
        case SpeciesType::Cocoa: os << "Cocoa"; break;
        case SpeciesType::RubberWall: os << "RubberWall"; break;
        case SpeciesType::SteelWall: os << "SteelWall"; break;
        case SpeciesType::MPT: os << "MPT"; break;
        case SpeciesType::APAPP: os << "APAPP"; break;
        case SpeciesType::PH101: os << "PH101"; break;
        case SpeciesType::SD100: os << "SD100"; break;
        default: logger(ERROR,"SpeciesType % not found", static_cast<unsigned>(type));
    }
    return os;
}

/*!
 * write from file
 */
std::istream& operator>>(std::istream& is, SpeciesType& type)
{
    std::string s;
    is >> s;
    switch (type) {
        case SpeciesType::Glass: type = SpeciesType::Glass; break;
        case SpeciesType::Cocoa: type = SpeciesType::Cocoa; break;
        case SpeciesType::RubberWall: type = SpeciesType::RubberWall; break;
        case SpeciesType::SteelWall: type = SpeciesType::SteelWall; break;
        case SpeciesType::MPT: type = SpeciesType::MPT; break;
        case SpeciesType::APAPP: type = SpeciesType::APAPP; break;
        case SpeciesType::PH101: type = SpeciesType::PH101; break;
        case SpeciesType::SD100: type = SpeciesType::SD100; break;
        default: logger(ERROR,"SpeciesType % not found", s);
    }
    return is;
}

const std::vector<PSD> PH101={{0.000674885,0},{0.000702057,0.00155301}, {0.000730324,0.0153487}, {0.000759728,0.0300549}, {0.000790317,0.0431186}, {0.000822137,0.0527894},  {0.000855238,0.059662}, {0.000889672,0.0650039}, {0.000925492,0.0690714},  {0.000962754,0.073889},  {0.00100152,0.0803085},  {0.00104184,0.0898127},   {0.00108379,0.102899},   {0.00112742,0.119774},   {0.00117282,0.140761},    {0.00122004,0.16621},   {0.00126916,0.195914},   {0.00132026,0.230241},    {0.00137341,0.26924},   {0.00142871,0.313096},   {0.00148623,0.361877},    {0.00154607,0.41568},   {0.00160832,0.473937},   {0.00167308,0.536015},   {0.00174044,0.600539},   {0.00181051,0.665645},    {0.00188341,0.72911},   {0.00195924,0.788703},   {0.00203812,0.842212},   {0.00212018,0.888111},   {0.00220554,0.925377},   {0.00229434,0.953944},   {0.00238672,0.974418},   {0.00248281,0.988042},   {0.00258278,0.996183},   {0.00268677,0.999729},          {0.00279494,1},};
const std::vector<PSD> APAPP={{0.000347498,0},{0.000359119,0.00444603}, {0.000371128,0.0104884}, {0.000383539,0.0177878}, {0.000396365,0.0293614},  {0.00040962,0.0428434}, {0.000423318,0.0591303}, {0.000437474,0.0775758}, {0.000452104,0.0998587},  {0.000467223,0.125095},  {0.000482847,0.152533},  {0.000498994,0.181535},  {0.000515681,0.210816},  {0.000532926,0.239386},   {0.000550748,0.26712},  {0.000569165,0.293364},   {0.000588199,0.31866},  {0.000607869,0.342496},  {0.000628196,0.365338},  {0.000649204,0.387162},  {0.000670914,0.408641},   {0.00069335,0.430263},  {0.000716537,0.452372},  {0.000740498,0.475192},  {0.000765261,0.498855},  {0.000790852,0.523625},  {0.000817299,0.549461},  {0.000844631,0.576458},  {0.000872876,0.604324},  {0.000902066,0.632797},   {0.000932232,0.66139},  {0.000963407,0.689945},  {0.000995624,0.718008},   {0.00102892,0.745323},   {0.00106333,0.771528},   {0.00109889,0.796434},   {0.00113563,0.819843},   {0.00117361,0.841626},    {0.00121286,0.86165},   {0.00125342,0.879897},   {0.00129533,0.896382},   {0.00133865,0.911061},   {0.00138342,0.924049},   {0.00142968,0.935429},   {0.00147749,0.945293},    {0.0015269,0.953796},   {0.00157796,0.961107},   {0.00163073,0.967351},   {0.00168526,0.972664},   {0.00174162,0.977186},   {0.00179986,0.981017},   {0.00186005,0.984265},    {0.00192225,0.98699},   {0.00198653,0.989288},   {0.00205296,0.991204},     {0.00212162,0.9928},   {0.00219257,0.994118},   {0.00226589,0.995198},    {0.00234166,0.99608},   {0.00241997,0.996799},    {0.0025009,0.997398},   {0.00258453,0.997899},   {0.00267096,0.998326},   {0.00276028,0.998693},   {0.00285258,0.999013},    {0.00294798,0.99929},   {0.00304656,0.999523},   {0.00314844,0.999714},   {0.00325373,0.999858},   {0.00336254,0.999961},          {0.00347498,1},};
const std::vector<PSD> SD100={{0.000583505,0},{0.000606998,0.00919727}, {0.000631437,0.0602607},   {0.00065686,0.124673},  {0.000683307,0.195593},  {0.000710818,0.263604},  {0.000739438,0.324019},  {0.000769209,0.374863},  {0.000800179,0.416516},  {0.000832396,0.450843},   {0.00086591,0.480544},   {0.000900774,0.50728},  {0.000937041,0.532593},  {0.000974769,0.557022},   {0.00101402,0.580942},   {0.00105484,0.604163},   {0.00109731,0.626542},   {0.00114149,0.647632},   {0.00118745,0.667134},   {0.00123526,0.684868},     {0.001285,0.700742},   {0.00133673,0.714844},   {0.00139055,0.727753},   {0.00144654,0.740155},   {0.00150478,0.752983},   {0.00156537,0.767261},   {0.00162839,0.783993},   {0.00169395,0.803724},    {0.00176216,0.82656},    {0.0018331,0.851939},   {0.00190691,0.878674},   {0.00198369,0.905289},    {0.00206355,0.93016},   {0.00214664,0.951864},   {0.00223307,0.969482},    {0.00232297,0.98265},    {0.0024165,0.991633},    {0.0025138,0.997072},   {0.00261501,0.999617},          {0.00272029,1},};
const std::vector<PSD> MPT={{0.000373296,0},{0.000385779,0.00253278},  {0.00039868,0.0117119}, {0.000412012,0.0236669},  {0.00042579,0.0392079}, {0.000440029,0.0579821}, {0.000454744,0.0811775},  {0.000469951,0.108148},  {0.000485667,0.138613},  {0.000501908,0.171966},  {0.000518693,0.208177},   {0.000536038,0.24665},  {0.000553964,0.287495},  {0.000572489,0.330315},  {0.000591634,0.374904},  {0.000611419,0.421188},  {0.000631865,0.468742},  {0.000652995,0.517312},  {0.000674832,0.566286},   {0.000697399,0.61516},  {0.000720721,0.663325},  {0.000744823,0.710131},   {0.00076973,0.754771},  {0.000795471,0.796516},  {0.000822072,0.834666},  {0.000849563,0.868697},  {0.000877974,0.898241},   {0.000907334,0.92316},  {0.000937676,0.943532},  {0.000969033,0.959633},    {0.00100144,0.97194},   {0.00103493,0.980992},    {0.00106954,0.98741},    {0.0011053,0.991798},   {0.00114227,0.994701},   {0.00118046,0.996557},   {0.00121994,0.997719},   {0.00126074,0.998426},    {0.0013029,0.998853},   {0.00134647,0.999107},   {0.00139149,0.999269},   {0.00143803,0.999366},   {0.00148612,0.999433},   {0.00153581,0.999493},   {0.00158717,0.999557},   {0.00164025,0.999623},     {0.0016951,0.99969},   {0.00175179,0.999757},   {0.00181037,0.999819},   {0.00187091,0.999869},   {0.00193348,0.999909},   {0.00199813,0.999941},   {0.00206495,0.999957},    {0.00213401,0.99998},          {0.00220537,1},};


Mdouble getDensity(SpeciesType type)
{
    switch (type)
    {
        case SpeciesType::Cocoa:
            return 1427;
        case SpeciesType::Glass:
            return 1500;
        case SpeciesType::MPT:
            return 1226.2;
        case SpeciesType::APAPP:
            return 1199.8;
        case SpeciesType::PH101:
            return 1473.7;
        case SpeciesType::SD100:
            return 1585.5;
        default:
            logger(ERROR,"Species density not defined");
    }
}

const bool hasPSD(SpeciesType type)
{
    switch (type)
    {
        case SpeciesType::MPT:
            return true;
        case SpeciesType::APAPP:
            return true;
        case SpeciesType::PH101:
            return true;
        case SpeciesType::SD100:
            return true;
        default:
            return false;
    }
}

const std::vector<PSD>& getPSD(SpeciesType type)
{
    switch (type)
    {
        case SpeciesType::MPT:
            return MPT;
        case SpeciesType::APAPP:
            return APAPP;
        case SpeciesType::PH101:
            return PH101;
        case SpeciesType::SD100:
            return SD100;
        default:
            logger(ERROR,"Species type not defined");
    }
}

Mdouble getMedianParticleRadius(const std::vector<PSD>&  psd)
{
    //const std::vector<PSD>&  psd = getPSD(type);
    double probability = 0.5;
    auto high = std::lower_bound(psd.begin(),psd.end(),probability);
    auto low = std::max(psd.begin(),high-1);
    Mdouble a = (probability - low->probability)/(high->probability -low->probability);
    return a*low->radius + (1-a)*high->radius;
}

Mdouble getMinParticleRadius(const std::vector<PSD>&  psd)
{
    return psd.front().radius;
}

Mdouble getMaxParticleRadius(const std::vector<PSD>&  psd)
{
    return psd.back().radius;
}

ParticleSpecies* addSingleSpecies(SpeciesType type, Mdouble medianParticleRadius, Mdouble minParticleRadius, DPMBase& dpm, bool adjustTimeStep=false)
{
    auto&  speciesHandler = dpm.speciesHandler;
    
    if (type == SpeciesType::Glass)
    {
        //define the particle properties
        const Mdouble density = getDensity(type); //kg/m^3
        const Mdouble collisionTime = 0.003; // s
        const Mdouble restitution = 0.2;
        const Mdouble friction = 0.5;
        const Mdouble minMass = 4. / 3. * density * pi * cubic(minParticleRadius);
        
        LinearViscoelasticFrictionReversibleAdhesiveSpecies species;
        species.setDensity(density);
        species.setCollisionTimeAndRestitutionCoefficient(collisionTime, restitution, minMass);
        species.setSlidingFrictionCoefficient(friction);
        species.setSlidingStiffness(2.0 / 7.0 * species.getStiffness());
        species.setSlidingDissipation(2.0 / 7.0 * species.getDissipation());
        species.setRollingFrictionCoefficient(friction);
        species.setRollingStiffness(2.0 / 5.0 * species.getStiffness());
        species.setRollingDissipation(2.0 / 5.0 * species.getDissipation());
        species.setAdhesionStiffness(species.getStiffness());
        const Mdouble mass = 4. / 3. * density * pi * cubic(medianParticleRadius);
        species.setAdhesionForceMax(1 * mass);
        
        if (adjustTimeStep) {
            dpm.setTimeStep(0.04 * collisionTime);
            logger(INFO, "Set time step to %", dpm.getTimeStep());
        }
        logger(INFO,"Species % is Glass",speciesHandler.getSize());
        return speciesHandler.copyAndAddObject(species);
    } else if (type == SpeciesType::Cocoa) {
        //scale stiffness by mass
        const Mdouble density = getDensity(type); //kg/m^3
        const Mdouble unloadingStiffnessMax =
                24067 * mathsFunc::cubic(medianParticleRadius / 0.0025); //for 2.5mm particles
        const Mdouble loadingStiffness = 0.2 * unloadingStiffnessMax; // s
        const Mdouble cohesionStiffness = 0.873 * unloadingStiffnessMax; // s
        const Mdouble penetrationDepthMax = 0.05; //check
        const Mdouble friction = 0.5;
        const Mdouble mass = 4. / 3. * density * pi * cubic(medianParticleRadius);
        const Mdouble minMass = 4. / 3. * density * pi * cubic(minParticleRadius);
        
        LinearPlasticViscoelasticFrictionSpecies species;
        species.setDensity(density);
        species.setPlasticParameters(loadingStiffness, unloadingStiffnessMax, cohesionStiffness,
                                     penetrationDepthMax);
        species.setRestitutionCoefficient(0.45, mass);
        species.setSlidingFrictionCoefficient(friction);
        species.setSlidingStiffness(2.0 / 7.0 * species.getLoadingStiffness());
        species.setSlidingDissipation(2.0 / 7.0 * species.getDissipation());
        species.setRollingFrictionCoefficient(friction);
        species.setRollingStiffness(0.2 * species.getLoadingStiffness());
        species.setRollingDissipation(0.2 * species.getDissipation());
        
        if (adjustTimeStep) {
            dpm.setTimeStep(0.04 * species.getCollisionTime(minMass));
            logger(INFO, "Set time step to %", dpm.getTimeStep());
        }
        logger(INFO, "Species % is Cocoa", speciesHandler.getSize());
        return speciesHandler.copyAndAddObject(species);
    } else if (type == SpeciesType::RubberWall) {
        // check the particle species has been defined already
        logger.assert_always(dpm.speciesHandler.getSize() != 0, "RubberWall can only be used for walls, not particles");
        // make a copy of the particle species
        auto species = dpm.speciesHandler.copyAndAddObject(dpm.speciesHandler.getObject(0));
        // we need the particle mass
        const Mdouble mass = 4. / 3. * species->getDensity() * pi * cubic(minParticleRadius);
        // now change its properties: r=0.1, mur=0.5
        auto cSpecies = dynamic_cast<LinearPlasticViscoelasticFrictionSpecies*>(species);
        auto gSpecies = dynamic_cast<LinearViscoelasticFrictionReversibleAdhesiveSpecies*>(species);
        if (cSpecies != nullptr)
        {
            cSpecies->setRestitutionCoefficient(0.1, mass);
        }
        else if (gSpecies != nullptr)
        {
            gSpecies->setRestitutionCoefficient(0.1, mass);
        }
        else
        {
            logger(ERROR, "Particle species must be cocoa or glass");
        }
        logger(INFO, "Species % is Rubber", speciesHandler.getSize() - 1);
        return species;
    }
    else if (type == SpeciesType::SteelWall)
    {
        // check the particle species has been defined already
        logger.assert_always(speciesHandler.getSize() != 0, "SteelWall can only be used for walls, not particles");
        // make a copy of the particle species
        auto species = speciesHandler.copyAndAddObject(speciesHandler.getObject(0));
        // now change its properties: mu=0.1, mur=0
        auto cSpecies = dynamic_cast<LinearPlasticViscoelasticFrictionSpecies*>(species);
        auto gSpecies = dynamic_cast<LinearViscoelasticFrictionReversibleAdhesiveSpecies*>(species);
        if (cSpecies != nullptr)
        {
            cSpecies->setSlidingFrictionCoefficient(0.2);
            cSpecies->setRollingFrictionCoefficient(0);
        }
        else if (gSpecies != nullptr)
        {
            gSpecies->setSlidingFrictionCoefficient(0.2);
        }
        else
        {
            logger(ERROR, "Particle species must be cocoa or glass");
        }
        logger(INFO, "Species % is Steel", speciesHandler.getSize() - 1);
        return species;
    } else {
        LinearPlasticViscoelasticFrictionSpecies species;
        species.setHandler(&dpm.speciesHandler);
        species.setDensity(getDensity(type));
        species.setConstantRestitution(true);
        const Mdouble unloadingStiffnessMax = 24067/species.getMassFromRadius(0.0025);  //scaled to 2.5 mm particles
        const Mdouble loadingStiffness = 0.2 * unloadingStiffnessMax; // s
        const Mdouble penetrationDepthMax = 0.05*0.0025; //scaled to 2.5 mm particles
        //default values are nearly cohesionless cocoa
        Mdouble restitution = 0.45;
        Mdouble cohesionStiffness = 1e-2 * unloadingStiffnessMax; //0.873 * unloadingStiffnessMax;
        Mdouble slidingFriction = 0.5;
        Mdouble rollingFriction = 0.1;
        if (type == SpeciesType::APAPP) {
            restitution = 6.918283e-01;
            cohesionStiffness = 1.865953e+00 * unloadingStiffnessMax;
            slidingFriction = 7.735157e-01;
            rollingFriction = 4.569588e-01;
            cohesionStiffness = 0.2 * unloadingStiffnessMax;
        } else if (type == SpeciesType::MPT) {
            restitution = 7.793266e-01;
            cohesionStiffness = 1.909827e+00 * unloadingStiffnessMax;
            slidingFriction = 6.889219e-01;
            rollingFriction = 6.140204e-01;
            cohesionStiffness = 0.2 * unloadingStiffnessMax;
        } else if (type == SpeciesType::PH101) {
            restitution = 8.963799e-01;
            cohesionStiffness = 2.457027e-01 * unloadingStiffnessMax;
            slidingFriction = 4.787953e-01;
            rollingFriction = 3.424929e-01;
        } else if (type == SpeciesType::SD100) {
            restitution = 9.726155e-01;
            cohesionStiffness = 6.858071e-04 * unloadingStiffnessMax;
            slidingFriction = 5.222745e-01;
            rollingFriction = 9.557100e-01;
        }
        species.setPlasticParameters(loadingStiffness, unloadingStiffnessMax, cohesionStiffness,
                                     penetrationDepthMax);
        species.setRestitutionCoefficient(restitution, 1.0);
        species.setSlidingFrictionCoefficient(slidingFriction);
        species.setSlidingStiffness(2.0 / 7.0 * species.getLoadingStiffness());
        species.setSlidingDissipation(2.0 / 7.0 * species.getDissipation());
        species.setRollingFrictionCoefficient(rollingFriction);
        species.setRollingStiffness(0.2 * species.getLoadingStiffness());
        species.setRollingDissipation(0.2 * species.getDissipation());
        species.setTorsionFrictionCoefficient(0.0);
        species.setTorsionStiffness(0.2 * species.getLoadingStiffness());
        species.setTorsionDissipation(0.2 * species.getDissipation());
        
        if (adjustTimeStep) { //if particle species
            dpm.setTimeStep(0.04 * species.getCollisionTime(1.0));
            logger(INFO, "Set time step to %", dpm.getTimeStep());
        }
        logger(INFO, "Species % is % (r=%, eps=%)", speciesHandler.getSize(),type,medianParticleRadius,species.getRestitutionCoefficient(1.0));
        return speciesHandler.copyAndAddObject(species);
    }
}

void modifySpecies(ParticleSpecies* s, Mdouble restitutionCoefficient, Mdouble relativeCohesionStiffness, Mdouble slidingFriction, Mdouble rollingFriction, SpeciesType speciesType) {
    const auto species = dynamic_cast<LinearPlasticViscoelasticFrictionSpecies*>(s);
    logger.assert_always(species!= nullptr, "Species % could not be cast to LinearPlasticViscoelasticFrictionSpecies", s->getName());
    const Mdouble radius = getMinParticleRadius(getPSD(speciesType));
    const Mdouble mass = s->getMassFromRadius(radius);
    //logger(INFO, "mass %", mass);
    species->setRestitutionCoefficient(restitutionCoefficient, mass);
    species->setCohesionStiffness(relativeCohesionStiffness*species->getUnloadingStiffnessMax());
    species->setSlidingFrictionCoefficient(slidingFriction);
    species->setRollingFrictionCoefficient(rollingFriction);
}

#endif //MERCURY_ADDSPECIES_H
