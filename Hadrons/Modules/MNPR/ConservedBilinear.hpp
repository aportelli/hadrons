/*
 * ConservedBilinear.hpp, part of Hadrons (https://github.com/aportelli/Hadrons)
 *
 * Copyright (C) 2015 - 2022
 *
 * Author: Antonin Portelli <antonin.portelli@me.com>
 * Author: Julia Kettle J.R.Kettle-2@sms.ed.ac.uk
 * Author: Peter Boyle <paboyle@ph.ed.ac.uk>
 * Author: Fabian Joswig <fabian.joswig@wwu.de>
 * Author: Felix Erben <felix.erben@ed.ac.uk>
 *
 * Hadrons is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 2 of the License, or
 * (at your option) any later version.
 *
 * Hadrons is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Hadrons.  If not, see <http://www.gnu.org/licenses/>.
 *
 * See the full license in the file "LICENSE" in the top level distribution
 * directory.
 */

/*  END LEGAL */

#ifndef Hadrons_ConservedBilinear_hpp_
#define Hadrons_ConservedBilinear_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <Hadrons/Modules/MNPR/NPRUtils.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                                TConservedBilinear                                       *
        Performs bilinear contractions of the type tr[g5*adj(qOut')*g5*G*qIn']
        Suitable for non exceptional momenta in Rome-Southampton NPR where G is
        the non-local conserved vector current which depends on the action.

Returns a spin-colour matrix, with indices si,sj, ci,cj
**************************************************************************/
BEGIN_MODULE_NAMESPACE(MNPR)

class ConservedBilinearPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(ConservedBilinearPar,
                                    std::string,    action,
                                    std::string,    qIn,
                                    std::string,    qOut,
                                    std::string,    pIn,
                                    std::string,    pOut,
                                    std::string,    output);
};

template <typename FImpl>
class TConservedBilinear: public Module<ConservedBilinearPar>
{
public:
    FERM_TYPE_ALIASES(FImpl,);
    class Metadata: Serializable
    {
    public:
        GRID_SERIALIZABLE_CLASS_MEMBERS(Metadata,
                                        Gamma::Algebra, gamma,
                                        std::string,  pIn,
                                        std::string,  pOut);
    };
    typedef Correlator<Metadata, SpinColourMatrix> Result;
public:
    // constructor
    TConservedBilinear(const std::string name);
    // destructor
    virtual ~TConservedBilinear(void) {};
    // dependencies/products
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
private:
};

MODULE_REGISTER_TMP(ConservedBilinear, ARG(TConservedBilinear<FIMPL>), MNPR);

/******************************************************************************
 *                           TConservedBilinear implementation                            *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl>
TConservedBilinear<FImpl>::TConservedBilinear(const std::string name)
: Module<ConservedBilinearPar>(name)
{}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl>
void TConservedBilinear<FImpl>::setup(void)
{
    LOG(Message) << "Running setup for ConservedBilinear" << std::endl;

    // The propagator can be 4d or 5d, but must match the action
    const unsigned int ActionLs_{ env().getObjectLs(par().action) };
    if (ActionLs_ != 1)
    {
        std::string sError{ "ConservedBilinear currently only implemented for 4d actions."};
        HADRONS_ERROR(Size, sError);
    }

    envTmpLat(PropagatorField, "qIn_phased");
    envTmpLat(PropagatorField, "qOut_phased");
    envTmpLat(PropagatorField, "lret");
    envTmpLat(PropagatorField, "dummy_phys_source");
    envTmpLat(ComplexField, "pDotXIn");
    envTmpLat(ComplexField, "pDotXOut");
    envTmpLat(ComplexField, "xMu");
}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl>
std::vector<std::string> TConservedBilinear<FImpl>::getInput(void)
{
    std::vector<std::string> input = {par().action, par().qIn, par().qOut};

    return input;
}

template <typename FImpl>
std::vector<std::string> TConservedBilinear<FImpl>::getOutput(void)
{
    std::vector<std::string> out = {};

    return out;
}

template <typename FImpl>
void TConservedBilinear<FImpl>::execute(void)
{

    LOG(Message) << "Computing conserved bilinear contractions '" << getName() << "' using"
                 << " propagators '" << par().qIn << "' and '" << par().qOut << "'"
                 << "with action " << par().action << std::endl;

    auto &act = envGet(FMat, par().action);

    // Propagators
    auto  &qIn    = envGet(PropagatorField, par().qIn);
    auto  &qOut   = envGet(PropagatorField, par().qOut);
    envGetTmp(PropagatorField, qIn_phased);
    envGetTmp(PropagatorField, qOut_phased);
    envGetTmp(PropagatorField, lret);
    envGetTmp(PropagatorField, dummy_phys_source);
    envGetTmp(ComplexField, pDotXIn);
    envGetTmp(ComplexField, pDotXOut);
    envGetTmp(ComplexField, xMu);

    // momentum on legs
    std::vector<Real>           pIn  = strToVec<Real>(par().pIn),
                                pOut = strToVec<Real>(par().pOut);
    Coordinate                  latt_size = GridDefaultLatt();
    Complex                     Ci(0.0,1.0);
    std::vector<Result>         result;
    Result                      r;

    Real volume = 1.0;
    for (int mu = 0; mu < Nd; mu++) {
        volume *= latt_size[mu];
    }

    NPRUtils<FImpl>::dot(pDotXIn,pIn);
    qIn_phased  = qIn  * exp(-Ci * pDotXIn);
    NPRUtils<FImpl>::dot(pDotXOut,pOut);
    qOut_phased = qOut * exp(-Ci * pDotXOut);


    r.info.pIn  = par().pIn;
    r.info.pOut = par().pOut;

    Gamma::Algebra Gmu [] = {
        Gamma::Algebra::GammaX,
        Gamma::Algebra::GammaY,
        Gamma::Algebra::GammaZ,
        Gamma::Algebra::GammaT,
    };

    for (int mu = 0; mu < Nd; mu++)
    {
        r.info.gamma = Gmu[mu];
        act.ContractConservedCurrent(qOut_phased, qIn_phased, lret, dummy_phys_source, Current::Vector, mu);
        r.corr.push_back( (1.0 / volume) * sum_large(lret) );
        result.push_back(r);
        r.corr.erase(r.corr.begin());
    }

    //////////////////////////////////////////////////
    saveResult(par().output, "ConservedBilinear", result);
    LOG(Message) << "Complete. Writing results to " << par().output << std::endl;
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_ConservedBilinear_hpp_