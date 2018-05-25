#include "HeaterBoundary.h"
#include "ParticleHandler.h"
#include "Particles/BaseParticle.h"

HeaterBoundary::HeaterBoundary() : BaseBoundary()
{
    posMin_ = Vec3D(0, 0, 0);
    posMax_ = Vec3D(0, 0, 0);
    specificHeatStrength_ = Vec3D(0, 0, 0);
}

HeaterBoundary::HeaterBoundary(const HeaterBoundary& other) : BaseBoundary(other)
{
    posMin_ = other.posMin_;
    posMax_ = other.posMax_;
    specificHeatStrength_ = other.specificHeatStrength_;
}

HeaterBoundary::~HeaterBoundary()
{
    logger(VERBOSE, "A HeaterBoundary has been destroyed.");
}

HeaterBoundary* HeaterBoundary::copy() const
{
    return new HeaterBoundary(*this);
}

void HeaterBoundary::set(Vec3D posMin, Vec3D posMax, Vec3D specificHeatStrength)
{
    posMin_ = posMin;
    posMax_ = posMax;
    setStrength(specificHeatStrength);
}


void HeaterBoundary::set2D(Vec3D posMin, Vec3D posMax, Mdouble specificHeatStrength)
{
    posMin_ = posMin;
    posMax_ = posMax;
    setStrength2D(specificHeatStrength);
}

void HeaterBoundary::set3D(Vec3D posMin, Vec3D posMax, Mdouble specificHeatStrength)
{
    posMin_ = posMin;
    posMax_ = posMax;
    setStrength3D(specificHeatStrength);
}

void HeaterBoundary::setStrength(Vec3D specificHeatStrength)
{
    specificHeatStrength_ = specificHeatStrength;
}

void HeaterBoundary::setStrength2D(Mdouble specificHeatStrength)
{
    specificHeatStrength_.X = specificHeatStrength / sqrt(2.);
    specificHeatStrength_.Y = specificHeatStrength / sqrt(2.);
    specificHeatStrength_.Z = 0;
}

void HeaterBoundary::setStrength3D(Mdouble specificHeatStrength)
{
    specificHeatStrength_.X = specificHeatStrength / sqrt(3.);
    specificHeatStrength_.Y = specificHeatStrength / sqrt(3.);
    specificHeatStrength_.Z = specificHeatStrength / sqrt(3.);
}

inline Mdouble HeaterBoundary::getVolume() const
{
    return ((posMax_.X - posMin_.X) * (posMax_.Y - posMin_.Y) * (posMax_.Z - posMin_.Z));
}

/*! 
 * \todo JMFT: Calculate the distance properly, not just 1 or -1.
 */
Mdouble HeaterBoundary::getDistance(const Vec3D& position) const
{
    // std::cerr << "Checking a particle position " << position << " against " << posMin_ << " and " << posMax_ << std::endl;
    if (posMin_.X <= position.X && position.X <= posMax_.X
        && posMin_.Y <= position.Y && position.Y <= posMax_.Y
        && posMin_.Z <= position.Z && position.Z <= posMax_.Z
            )
    {
        // std::cerr << "Yes\n";
        return -1;
    }
    else
    {
        // std::cerr << "No\n";
        return 1;
    }
}

void HeaterBoundary::checkBoundaryAfterParticlesMove(ParticleHandler& pH)
{
    for (auto p = pH.begin(); p != pH.end(); ++p)
        checkBoundaryAfterParticleMoved(*p, pH);
}


/*!
 * \todo JMFT: Make a 2D version of this 
 */
bool HeaterBoundary::checkBoundaryAfterParticleMoved(BaseParticle* p, ParticleHandler& pH)
{
    if (!(p->isFixed()) && (getDistance(p->getPosition()) < 0))
    {
        // fprintf(stderr, "heating a particle\n");
        // Mdouble noise = heatStrength_ / p->getMass() * pow(pH.getDPMBase()->getTimeStep(),0.5) ;
        // We need a factor of sqrt(3) in front because we are generating random
        // numbers using a uniformly distributed random variable, not a normal one.
        Mdouble noiseX = sqrt(2. * pH.getDPMBase()->getTimeStep() * specificHeatStrength_.X);
        Mdouble noiseY = sqrt(2. * pH.getDPMBase()->getTimeStep() * specificHeatStrength_.Y);
        Mdouble noiseZ = sqrt(2. * pH.getDPMBase()->getTimeStep() * specificHeatStrength_.Z);
        
        Vec3D brownianMotion = Vec3D(
                noiseX * pH.getDPMBase()->random.getRandomNumber(-1, 1),
                noiseY * pH.getDPMBase()->random.getRandomNumber(-1, 1),
                noiseZ * pH.getDPMBase()->random.getRandomNumber(-1, 1));
        p->setVelocity(p->getVelocity() + brownianMotion);
    }
    return true;
}

/*!
 * \details Reads a number of boundary properties from the given std::istream.
 * \param[in,out] is   the istream
 */
void HeaterBoundary::read(std::istream& is)
{
    BaseBoundary::read(is);
    std::string dummy;
    is >> dummy >> posMin_
       >> dummy >> posMax_
       >> dummy >> specificHeatStrength_;
}

/*!
 * \details Writes the boundary properties to an std::ostream. 
 * \param[out] os   the ostream the properties are to be written to.
 */
void HeaterBoundary::write(std::ostream& os) const
{
    BaseBoundary::write(os);
    os << " posMin " << posMin_
       << " posMax " << posMax_
       << " specificHeatStrength " << specificHeatStrength_;
}

/*!
 * \details Returns the object's class name (i.e. 'HeaterBoundary').
 * \return the object's class name
 */
std::string HeaterBoundary::getName() const
{
    return "HeaterBoundary";
}

