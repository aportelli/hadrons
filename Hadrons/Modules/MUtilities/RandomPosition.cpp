#include <Hadrons/Modules/MUtilities/RandomPosition.hpp>

using namespace Grid;
using namespace Hadrons;
using namespace MUtilities;

/******************************************************************************
*                  TRandomPosition implementation                             *
******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
TRandomPosition::TRandomPosition(const std::string name)
: Module<RandomPositionPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
std::vector<std::string> TRandomPosition::getInput(void)
{
    std::vector<std::string> in;
    
    return in;
}

std::vector<std::string> TRandomPosition::getOutput(void)
{
    std::vector<std::string> out = {getName()};
    
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
void TRandomPosition::setup(void)
{
    envCreate(HadronsSerializable, getName(), 1, 0);
}

uint32_t TRandomPosition::randInt(int max)
{
    // 1) Generate a double in the interval [0,1].
    // 2) Convert to an int, rounding down, using 48 / 52 available bits from a double-precision number.
    // 3) This has 0 bias for power-of-2 dimension size. For a size of 96, it has a modulo bias of ~2.27e-11%.
    auto &rng = rngSerial();
    RealD randfloat;
    random(rng, randfloat);
    size_t max_size = ((size_t)1 << 48);
    size_t randint = static_cast<size_t>(randfloat * max_size);
    return randint % abs(max);
}

// execution ///////////////////////////////////////////////////////////////////
void TRandomPosition::execute(void)
{
    std::vector<int> dims = env().getDim();

    std::vector<int> pos(env().getNd());
    for (int i=0; i < dims.size(); ++i)
        pos[i] = this->randInt(dims[i]);

    LOG(Message) << "Created random position vector " << pos << std::endl;
    saveResult(par().output, "pos", pos);
    auto &out = envGet(HadronsSerializable, getName());
    out = pos;
}
