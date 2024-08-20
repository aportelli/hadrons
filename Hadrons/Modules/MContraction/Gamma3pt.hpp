/*
 * Gamma3pt.hpp, part of Hadrons (https://github.com/aportelli/Hadrons)
 *
 * Copyright (C) 2015 - 2023
 *
 * Author: Antonin Portelli <antonin.portelli@me.com>
 * Author: Fionn O hOgain <fionn.o.hogain@ed.ac.uk>
 * Author: Lanny91 <andrew.lawson@gmail.com>
 * Author: Ryan Hill <rchrys.hill@gmail.com>
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

#ifndef Hadrons_MContraction_Gamma3pt_hpp_
#define Hadrons_MContraction_Gamma3pt_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <Hadrons/Serialization.hpp>

BEGIN_HADRONS_NAMESPACE

/*
 * 3pt contraction with gamma matrix insertion.
 *
 * Schematic:
 *
 *             q2           q3
 *        /----<------*------<----¬
 *       /          gamma          \
 *      /                           \
 *   i *                            * f
 *      \                          /
 *       \                        /
 *        \----------->----------/
 *                   q1
 *
 *      trace(g5*q1*adj(q2)*g5*gamma*q3)
 * 
 *  options:
 *   - q1: sink smeared propagator, source at i
 *   - q2: propagator, source at i
 *   - q3: propagator, source at f
 *   - gammas: gamma matrices to insert
 *             (space-separated strings e.g. "GammaT GammaX GammaY") 
 *   - tSnk: sink position for propagator q1.
 *
 */

/* 11/07/2024 (AP): This module has been extended to allow different gamma
 * matrices at the source and the sink
 * ----------------------------------------------------
 * Previous version input parameters
 * <options>
 *  <q1>q1</q1>
 *  <q2>q2</q2>
 *  <q3>q3</q3>
 *  <gamma>all</gamma>
 *  <tSnk>16</tSnk>
 *  <output>3pt</output>
 * </options>
 * 
 * is equivalent to the code below in the new version
 * 
 * <options>
 *  <q1>q1</q1>
 *  <q2>q2</q2>
 *  <q3>q3</q3>
 *  <gamma>
 *    <vertex>all</vertex>
 *    <source>Gamma5</source>
 *    <sink>Gamma5</sink>
 *  </gamma>
 *  <tSnk>16</tSnk>
 *  <save4d>false</save4d>
 *  <output>3pt</output>
 * </options>
 * 
 * HDF5 metadata is changed as well to reflect that, hopefully that part
 * is self-explanatory.
 */

/******************************************************************************
 *                               Gamma3pt                                     *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MContraction)

class Gamma3ptPar: Serializable
{
public:
    class GammaPar: Serializable
    {
    public:
        GRID_SERIALIZABLE_CLASS_MEMBERS(GammaPar,
                                        std::string, vertex,
                                        std::string, source,
                                        std::string, sink);
    };
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(Gamma3ptPar,
                                    std::string,  q1,
                                    std::string,  q2,
                                    std::string,  q3,
                                    GammaPar,  gamma,
                                    unsigned int, tSnk,
                                    bool, save4d,
                                    std::string,  output);
};

template <typename FImpl1, typename FImpl2, typename FImpl3>
class TGamma3pt: public Module<Gamma3ptPar>
{
    FERM_TYPE_ALIASES(FImpl1, 1);
    FERM_TYPE_ALIASES(FImpl2, 2);
    FERM_TYPE_ALIASES(FImpl3, 3); 
public:
    class GammaTriad: Serializable
    {
    public:
        GRID_SERIALIZABLE_CLASS_MEMBERS(GammaTriad,
                                        Gamma::Algebra, vertex, 
                                        Gamma::Algebra, source, 
                                        Gamma::Algebra, sink);
    };
    class Result: Serializable
    {
    public:
        GRID_SERIALIZABLE_CLASS_MEMBERS(Result,
                                        GammaTriad, gamma,
                                        std::vector<Complex>, corr);
    };
public:
    // constructor
    TGamma3pt(const std::string name);
    // destructor
    virtual ~TGamma3pt(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    virtual std::vector<std::string> getOutputFiles(void);
    static std::vector<Gamma::Algebra> parseGammaString(const std::string &gammaString);
protected:
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
private:
    std::vector<GammaTriad> gammaList_;
};

MODULE_REGISTER_TMP(Gamma3pt, ARG(TGamma3pt<FIMPL, FIMPL, FIMPL>), MContraction);

/******************************************************************************
 *                       TGamma3pt implementation                             *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl1, typename FImpl2, typename FImpl3>
TGamma3pt<FImpl1, FImpl2, FImpl3>::TGamma3pt(const std::string name)
: Module<Gamma3ptPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl1, typename FImpl2, typename FImpl3>
std::vector<std::string> TGamma3pt<FImpl1, FImpl2, FImpl3>::getInput(void)
{
    std::vector<std::string> in = {par().q1, par().q2, par().q3};
    
    return in;
}

template <typename FImpl1, typename FImpl2, typename FImpl3>
std::vector<std::string> TGamma3pt<FImpl1, FImpl2, FImpl3>::getOutput(void)
{
    std::vector<std::string> out = {getName()};

    if (par().save4d)
    {
        out.push_back(getName() + "_4d");
    }
    
    return out;
}

template <typename FImpl1, typename FImpl2, typename FImpl3>
std::vector<std::string> TGamma3pt<FImpl1, FImpl2, FImpl3>::getOutputFiles(void)
{
    std::vector<std::string> output;
    
    if (!par().output.empty())
        output.push_back(resultFilename(par().output));
    
    return output;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl1, typename FImpl2, typename FImpl3>
void TGamma3pt<FImpl1, FImpl2, FImpl3>::setup(void)
{
    std::vector<Gamma::Algebra> vGammaList, srcGammaList, snkGammaList;

    vGammaList = parseGammaString(par().gamma.vertex);
    srcGammaList = parseGammaString(par().gamma.source);
    snkGammaList = parseGammaString(par().gamma.sink);
    gammaList_.clear();
    for (auto gV: vGammaList)
    for (auto gSrc: srcGammaList)
    for (auto gSnk: snkGammaList)
    {
        GammaTriad t;

        t.vertex = gV;
        t.source = gSrc;
        t.sink = gSnk;
        gammaList_.push_back(t);
    }
    if (par().save4d)
    {
        envCreate(std::vector<ComplexField1>, getName() + "_4d", 1, 
                  gammaList_.size(), envGetGrid(ComplexField1));
    }
    else
    {
        envTmpLat(LatticeComplex, "c");
    }
    envCreate(HadronsSerializable, getName(), 1, 0);
}

template <typename FImpl1, typename FImpl2, typename FImpl3>
std::vector<Gamma::Algebra> 
TGamma3pt<FImpl1, FImpl2, FImpl3>::parseGammaString(const std::string &gammaString)
{
    std::vector<Gamma::Algebra> gammaList;
    // Determine gamma matrices to insert
    if (gammaString == "all")
    {
        // Do all contractions.
        for (unsigned int i = 1; i < Gamma::nGamma; i += 2)
        {
            gammaList.push_back((Gamma::Algebra)i);
        }
    }
    else
    {
        // Parse individual contractions from input string.
        gammaList = strToVec<Gamma::Algebra>(gammaString);
    }

    return gammaList;
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl1, typename FImpl2, typename FImpl3>
void TGamma3pt<FImpl1, FImpl2, FImpl3>::execute(void)
{
    LOG(Message) << "Computing 3pt contractions '" << getName() << "' using"
                 << " quarks '" << par().q1 << "' (spectator), '" << par().q2 << "', and '"
                 << par().q3 << "'" << std::endl;
    LOG(Message) << "Vertices: " << par().gamma.vertex << std::endl;
    LOG(Message) << " Sources: " << par().gamma.source << std::endl;
    LOG(Message) << "   Sinks: " << par().gamma.sink << std::endl;
    LOG(Message) << "   Total: " << gammaList_.size() << " correlator(s)" << std::endl;
    if (par().save4d)
    {
        LOG(Message) << "4D correlators are kept in memory" << std::endl;
    }

    auto                        &q1 = envGet(SlicedPropagator1, par().q1);
    auto                        &q2 = envGet(PropagatorField2, par().q2);
    auto                        &q3 = envGet(PropagatorField2, par().q3);
    Gamma                       g5(Gamma::Algebra::Gamma5);
    ComplexField1               *pt;
    std::vector<TComplex>       buf;
    std::vector<Result>         result;
    int                         nt = env().getDim(Tp);
    SitePropagator1             q1Snk = q1[par().tSnk];

    for (unsigned int i = 0; i < gammaList_.size(); ++i)
    {
        if (par().save4d)
        {
            auto &vec = envGet(std::vector<ComplexField1>, getName() + "_4d");
            pt = &(vec[i]);
        }
        else
        {
            envGetTmp(ComplexField1, c);
            pt = &c;
        }

        Gamma gV(gammaList_[i].vertex), gSrc(gammaList_[i].source), 
              gSnk(gammaList_[i].sink);
        Result r;
        
        (*pt) = trace(gSnk*q1Snk*(adj(gSrc)*g5)*adj(q2)*(g5*gV)*q3);
        r.gamma = gammaList_[i];
        r.corr.resize(nt);
        sliceSum(*pt, buf, Tp);
        for (unsigned int t = 0; t < buf.size(); ++t)
        {
            r.corr[t] = TensorRemove(buf[t]);
        }
        result.push_back(r);
    }
    saveResult(par().output, "gamma3pt", result);
    auto &out = envGet(HadronsSerializable, getName());
    out = result;
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MContraction_Gamma3pt_hpp_
