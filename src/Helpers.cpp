#include "Helpers.h"

#include <filesystem>
#include <regex>

namespace MarlinACTS {

std::string findFile(const std::string& inpath)
{
    if (inpath.empty()) return inpath;

    if (inpath[0] == '/') return inpath;

    if (std::filesystem::exists(inpath)) return inpath;

    return inpath;
}

ZBinSchema parseZBinSchema(std::string schema_str)
{
    float pos = 0.0;
    int t1 = 1;
    int t2 = 1;
    int b1 = 1;
    int b2 = 1;
    bool err = true;

    try
    {
      std::regex separator { "," };
      std::sregex_token_iterator item(schema_str.begin(), schema_str.end(), separator, -1);
      std::sregex_token_iterator end_it;
      for(int k = 0; item != end_it; k++)
      {
        switch(k)
        {
        case 0:
          pos = std::stof(*item); break;
        case 1:
          t1 = std::stoi(*item); break;
        case 2:
          t2 = std::stoi(*item); break;
        case 3:
          b1 = std::stoi(*item); break;
        case 4:
          b2 = std::stoi(*item); err = false;
        }
        item++;
      }
    }
    catch (std::invalid_argument const& ex)
    {}
    catch (std::out_of_range const& ex)
    {}

    ZBinSchema result { pos, t1, t2, b1, b2, err };
    return result;
}

Acts::ParticleHypothesis getParticleHypothesis(const EVENT::MCParticle* mcParticle)
{
    switch (std::abs(mcParticle->getPDG())) {
    case 11:
        return Acts::ParticleHypothesis {Acts::PdgParticle::eElectron};
    case 13:
        return Acts::ParticleHypothesis {Acts::PdgParticle::eMuon};
    case 15:
        return Acts::ParticleHypothesis {Acts::PdgParticle::eTau};
    case 22:
        return Acts::ParticleHypothesis {Acts::PdgParticle::eGamma};
    case 111:
        return Acts::ParticleHypothesis {Acts::PdgParticle::ePionZero};
    case 211:
        return Acts::ParticleHypothesis {Acts::PdgParticle::ePionPlus};
    case 2112:
        return Acts::ParticleHypothesis {Acts::PdgParticle::eNeutron};
    case 2212:
        return Acts::ParticleHypothesis {Acts::PdgParticle::eProton};
    }

    Acts::PdgParticle pdg = Acts::PdgParticle::eInvalid;
    float mass = 0.0f;
    Acts::AnyCharge charge_type { 0.0f };
    return Acts::ParticleHypothesis { pdg, mass, charge_type };
}

} //namespace MarlinACTS

