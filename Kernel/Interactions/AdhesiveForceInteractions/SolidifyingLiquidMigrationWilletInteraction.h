//Copyright (c) 2013-2023, The MercuryDPM Developers Team. All rights reserved.
//For the list of developers, see <http://www.MercuryDPM.org/Team>.
//
//Redistribution and use in source and binary forms, with or without
//modification, are permitted provided that the following conditions are met:
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//  * Neither the name MercuryDPM nor the
//    names of its contributors may be used to endorse or promote products
//    derived from this software without specific prior written permission.
//
//THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
//ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
//WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
//DISCLAIMED. IN NO EVENT SHALL THE MERCURYDPM DEVELOPERS TEAM BE LIABLE FOR ANY
//DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
//(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
//LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
//ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
//(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
//SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

#ifndef SolidifyingLiquidMigrationWilletInteraction_h
#define SolidifyingLiquidMigrationWilletInteraction_h

#include "LiquidMigrationWilletInteraction.h"
#include "Math/Vector.h"

class BaseParticle;

class SolidifyingLiquidMigrationWilletSpecies;

class BaseInteractable;

/*!
 * \class SolidifyingLiquidMigrationWilletInteraction
 */
class SolidifyingLiquidMigrationWilletInteraction : public virtual LiquidMigrationWilletInteraction
{
public:
    /*!
     * \brief An alias name for SolidifyingLiquidMigrationWilletSpecies data type.
     */
    typedef SolidifyingLiquidMigrationWilletSpecies SpeciesType;

    /*!
     * \brief Constructor.
     */
    SolidifyingLiquidMigrationWilletInteraction(BaseInteractable* P, BaseInteractable* I, unsigned timeStamp);

    /*!
     * \brief Copy constructor.
     */
    SolidifyingLiquidMigrationWilletInteraction(const SolidifyingLiquidMigrationWilletInteraction& p);

    /// \todo MX: add to MPI
    SolidifyingLiquidMigrationWilletInteraction();

    /*!
     * \brief Destructor.
     */
    ~SolidifyingLiquidMigrationWilletInteraction() override;

    /*!
     * \brief Computes the Adhesive force.
     */
    void computeAdhesionForce();

    /*!
     * \brief Interaction read function, which accepts an std::istream as input.
     */
    void read(std::istream& is) override;

    /*!
     * \brief Interaction print function, which accepts an std::ostream as input.
     */
    void write(std::ostream& os) const override;

    /*!
     * \brief Returns a pointer to the adhesive force species SolidifyingLiquidMigrationWilletSpecies.
     */
    const SolidifyingLiquidMigrationWilletSpecies* getSpecies() const;

    /*!
     * \brief Returns the name of the interaction, see Interaction.h.
     */
    std::string getBaseName() const;

    /*!
     * \brief Returns the elastic energy stored in the adhesive spring. 
     */
    Mdouble getElasticEnergy() const override;

    /*!
     * \brief Get function.
     */
    bool getBonded() const;

    /*!
     * \brief Set function.
     */
    void setBonded(bool bonded);

    /*!
     * This access function create a bond between the particles with the stength of the maximum bondForce
     */
    void bond();

    /*!
     * This access function create a bond between the particles with the stength of the current force.
     * So the current distance should be come the equilibrum distance.
     * \bug Make have to things about dispative forces
     */
    void bondInPlace();

    void unbond();

    /*!
     * \brief Simple access function to return the force of the bond
     */
    Mdouble getBondForce() {
        return bondForce_;
    }

    void setBondForce(Mdouble bondForce) {
        bondForce_ = bondForce;
    }

    ///Different from bondInteraction:

    double getSolidVolume() { return solidVolume_; }

    void setSolidVolume(Mdouble solidVolume) {
        solidVolume_ = solidVolume;
    }

    void setLiquidBridgeVolume(Mdouble liquidBridgeVolume) override;

    void addLiquidBridgeVolume(Mdouble liquidBridgeVolume) override;

    /**
     * writes extra information to the VTK output
     */
    unsigned getNumberOfFieldsVTK() const override;

    std::string getTypeVTK(unsigned i) const override;

    std::string getNameVTK(unsigned i) const override;

    std::vector<Mdouble> getFieldVTK(unsigned i) const override;

    void rupture() {
        if (!getBonded()) LiquidMigrationWilletInteraction::rupture();
    }

private:
    /*!
     * \brief A history parameter to store if the particles were in contact or not. Useful
     *        to compute adhesive forces.
     */
    bool bonded_;

    /*!
     * \brief This stores the extra force created by the bond
     */
    Mdouble bondForce_;


    Mdouble solidVolume_;


};

#endif
