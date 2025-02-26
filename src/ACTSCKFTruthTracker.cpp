#include "ACTSCKFTruthTracker.h"
#include "Helpers.h"

#include <EVENT/MCParticle.h>

ACTSCKFTruthTracker aACTSCKFTruthTracker;

ACTSCKFTruthTracker::ACTSCKFTruthTracker() :
    ACTSCKFBaseTracker("ACTSCKFTruthTracker")
{
    registerInputCollection(LCIO::MCPARTICLE, "MCParticleCollectionName",
                            "Name of the MCParticle input collection (used for seeding).",
                            _inputParticleCollection, std::string("MCParticle"));
}

void ACTSCKFTruthTracker::init()
{
    ACTSCKFBaseTracker::init();
}

ACTSCKFBaseTracker::SeedParamList
ACTSCKFTruthTracker::getSeeds(const MarlinACTS::MeasurementContainer& m_list, LCEvent* evt)
{
    SeedParamList seeds;

    LCCollection* particleCollection = getCollection(_inputParticleCollection, evt);
    if (particleCollection == nullptr) return seeds;

    for (uint32_t idxP = 0; idxP < particleCollection->getNumberOfElements(); ++idxP)
    {
        const MCParticle* mcParticle =
            static_cast<const MCParticle*>(particleCollection->getElementAt(idxP));

        // Tracks are made by stable charged particles from generation
        if (mcParticle->isCreatedInSimulation() ||
            mcParticle->getGeneratorStatus() != 1 || mcParticle->getCharge() == 0) continue;

        // Create initial parameters
        double px = mcParticle->getMomentum()[0];
        double py = mcParticle->getMomentum()[1];
        double pz = mcParticle->getMomentum()[2];
        double pt = sqrt(px * px + py * py);
        double p = sqrt(px * px + py * py + pz * pz);

        Acts::BoundVector params = Acts::BoundVector::Zero();
        // position/time
        params[Acts::eBoundLoc0] = 0;
        params[Acts::eBoundLoc1] = 0;
        params[Acts::eBoundTime] = mcParticle->getTime();
        // direction angles phi,theta
        params[Acts::eBoundPhi] = atan2(py, px);
        params[Acts::eBoundTheta] = atan2(pt, pz);
        // absolute momentum vector
        params[Acts::eBoundQOverP] = mcParticle->getCharge() / p;

        // build the track covariance matrix using the smearing sigmas
        Acts::BoundSquareMatrix cov = Acts::BoundSquareMatrix::Zero();
        cov(Acts::eBoundLoc0, Acts::eBoundLoc0) = std::pow(_initialTrackError_d0, 2);
        cov(Acts::eBoundLoc1, Acts::eBoundLoc1) = std::pow(_initialTrackError_z0, 2);
        cov(Acts::eBoundTime, Acts::eBoundTime) = std::pow(_initialTrackError_time, 2);
        cov(Acts::eBoundPhi, Acts::eBoundPhi) = std::pow(_initialTrackError_phi, 2);
        cov(Acts::eBoundTheta, Acts::eBoundTheta) = std::pow(_initialTrackError_lambda, 2);
        cov(Acts::eBoundQOverP, Acts::eBoundQOverP) = std::pow(_initialTrackError_relP * p / (p * p), 2);

        Acts::BoundTrackParameters seed(perigeeSurface, params, cov,
                                        MarlinACTS::getParticleHypothesis(mcParticle));
        seeds.push_back(seed);
    }
    return seeds;
}
