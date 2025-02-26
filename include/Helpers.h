#ifndef HELPERS_H
#define HELPERS_H 1

#include <string>
#include <tuple>

#include <EVENT/MCParticle.h>
#include <Acts/EventData/ParticleHypothesis.hpp>

namespace MarlinACTS {

std::string findFile(const std::string& inpath);

using ZBinSchema = std::tuple<float, int, int, int, int, bool>;

ZBinSchema parseZBinSchema(std::string schema_str);

Acts::ParticleHypothesis getParticleHypothesis(const EVENT::MCParticle* mcParticle);

} //namespace MarlinACTS

#endif //HELPERS_H
