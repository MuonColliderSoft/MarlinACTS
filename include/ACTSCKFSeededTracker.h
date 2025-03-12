#ifndef ACTSCKFSEEDEDTRACKER_H
#define ACTSCKFSEEDEDTRACKER_H 1

#include "ACTSCKFBaseTracker.h"
#include "SeedSpacePoint.h"
#include "SpacePointContainer.h"
#include "GeometryIdSelector.h"

#include <Acts/Seeding/SeedFinder.hpp>
#include <Acts/Seeding/SpacePointGrid.hpp>

class ACTSCKFSeededTracker : public ACTSCKFBaseTracker
{
public:
    ACTSCKFSeededTracker(const ACTSCKFSeededTracker&) = delete;
    ACTSCKFSeededTracker& operator=(const ACTSCKFSeededTracker&) = delete;
    ACTSCKFSeededTracker();

    virtual marlin::Processor *newProcessor() { return new ACTSCKFSeededTracker; }

    virtual void init() override;

    virtual void processEvent(LCEvent* evt) override;

protected:
    using SSPoint = typename Acts::SpacePointContainer<
            MarlinACTS::SpacePointContainer<std::vector<const MarlinACTS::SeedSpacePoint*>>,
            Acts::detail::RefHolder>::SpacePointProxyType;
    using SSPointGrid = Acts::CylindricalSpacePointGrid<SSPoint>;

    virtual SeedParamList getSeeds(const MarlinACTS::MeasurementContainer& m_list,
                                   LCEvent* evt) override;

private:
    float _seedFinding_rMax;
    float _seedFinding_deltaRMin;
    float _seedFinding_deltaRMax;
    float _seedFinding_deltaRMinTop;
    float _seedFinding_deltaRMaxTop;
    float _seedFinding_deltaRMinBottom;
    float _seedFinding_deltaRMaxBottom;
    float _seedFinding_collisionRegion;
    float _seedFinding_zMax;
    float _seedFinding_sigmaScattering;
    float _seedFinding_radLengthPerSeed;
    float _seedFinding_minPt;
    float _seedFinding_impactMax;
    float _seedFinding_cotThetaMax;

    std::string _outputSeedCollection;

    int _zTopBinLen;
    int _zBottomBinLen;
    int _phiTopBinLen;
    int _phiBottomBinLen;

    StringVec _seedFinding_zBinSchema;
    std::vector<float> _zBinEdges;
    std::vector<std::pair<int, int>> _ZTopBinSchema;
    std::vector<std::pair<int, int>> _ZBottomBinSchema;

    std::vector<std::string> _seedingLayers;

    Acts::SeedFinderConfig<SSPoint> finderCfg;
    Acts::CylindricalSpacePointGridConfig gridCfg;

    MarlinACTS::GeometryIdSelector _seedGeometrySelection;

    LCCollectionVec* seedCollection;
};

#endif //ACTSCKFSEEDEDTRACKER_H
