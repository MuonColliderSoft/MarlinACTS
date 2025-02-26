#include "ACTSCKFSeededTracker.h"
#include "Helpers.h"

#include <Acts/Definitions/Units.hpp>
#include <Acts/Seeding/EstimateTrackParamsFromSeed.hpp>

using SSPoint = MarlinACTS::SeedSpacePoint;
using SSPointGrid = Acts::CylindricalSpacePointGrid<SSPoint>;

ACTSCKFSeededTracker aACTSCKFSeededTracker;

ACTSCKFSeededTracker::ACTSCKFSeededTracker() :
    ACTSCKFBaseTracker("ACTSCKFSeededTracker")
{
    registerProcessorParameter("SeedFinding_RMax", "Maximum radius of hits to consider.",
                               _seedFinding_rMax, 150.f);

    registerProcessorParameter("SeedFinding_DeltaRMin", "Minimum dR between hits in a seed.",
                               _seedFinding_deltaRMin, 5.f);

    registerProcessorParameter("SeedFinding_DeltaRMax", "Maximum dR between hits in a seed.",
                               _seedFinding_deltaRMax, 80.f);

    registerProcessorParameter("SeedFinding_DeltaRMinTop",
                               "Minimum dR between the reference hit and outer ones in a seed.",
                               _seedFinding_deltaRMinTop, 0.f);

    registerProcessorParameter("SeedFinding_DeltaRMaxTop",
                               "Maximum dR between the reference hit and outer ones in a seed.",
                               _seedFinding_deltaRMaxTop, 0.f);

    registerProcessorParameter("SeedFinding_DeltaRMinBottom",
                               "Minimum dR between the reference hit and inner ones in a seed.",
                               _seedFinding_deltaRMinBottom, 0.f);

    registerProcessorParameter("SeedFinding_DeltaRMaxBottom",
                               "Maximum dR between the reference hit and inner ones in a seed.",
                               _seedFinding_deltaRMaxBottom, 0.f);

    registerProcessorParameter("SeedFinding_zTopBinLen", "Number of top bins along Z for seeding",
                               _zTopBinLen, 1);

    registerProcessorParameter("SeedFinding_zBottomBinLen", "Number of bottom bins along Z for seeding",
                               _zBottomBinLen, 1);

    registerProcessorParameter("SeedFinding_phiTopBinLen", "Number of top bins along phi for seeding",
                               _phiTopBinLen, 1);

    registerProcessorParameter("SeedFinding_phiBottomBinLen", "Number of bottom bins along phi for seeding",
                               _phiBottomBinLen, 1);

    registerProcessorParameter("SeedFinding_zBinSchema", "Bins schema along Z for seeding.",
                               _seedFinding_zBinSchema, StringVec(0));

    registerProcessorParameter("SeedFinding_CollisionRegion",
                               "Size of the collision region in one direction (assumed symmetric).",
                               _seedFinding_collisionRegion, 75.f);

    registerProcessorParameter("SeedFinding_ZMax", "Maximum z of hits hits to consider.",
                               _seedFinding_zMax, 600.f);

    registerProcessorParameter("SeedFinding_RadLengthPerSeed", "Average radiation length per seed.",
                               _seedFinding_radLengthPerSeed, 0.1f);

    registerProcessorParameter("SeedFinding_SigmaScattering", "Number of sigmas to allow in scattering angle.",
                               _seedFinding_sigmaScattering, 50.f);

    registerProcessorParameter("SeedFinding_MinPt", "Minimum pT of tracks to seed.",
                               _seedFinding_minPt, 500.f);

    registerProcessorParameter("SeedFinding_ImpactMax", "Maximum d0 of tracks to seed.",
                               _seedFinding_impactMax, static_cast<float>(3.f * Acts::UnitConstants::mm));

    registerProcessorParameter("SeedFinding_cotThetaMax", "Maximum cot theta for seed.",
                               _seedFinding_cotThetaMax, 7.40627f);

    registerProcessorParameter("SeedingLayers",
                               "Layers to use for seeding in format \"VolID LayID\", one per line. ID's "
                               "are ACTS GeometryID's. * can be used to wildcard.",
                               _seedingLayers, {"*", "*"});

}

void ACTSCKFSeededTracker::init()
{
    ACTSCKFBaseTracker::init();

    // Initialize seeding layers
    std::vector<std::string> seedingLayers;
    std::copy_if(_seedingLayers.begin(), _seedingLayers.end(),
                 std::back_inserter(seedingLayers),
                 [](const std::string &s) { return !s.empty(); });

    if (seedingLayers.size() % 2 != 0)
    {
        throw std::runtime_error("SeedingLayers needs an even number of entries");
    }

    std::vector<Acts::GeometryIdentifier> geoSelection;
    for (uint32_t i = 0; i < seedingLayers.size(); i += 2)
    {
        Acts::GeometryIdentifier geoid;
        if (_seedingLayers[i + 0] != "*")
            // volume
            geoid = geoid.setVolume(std::stoi(_seedingLayers[i + 0]));
        if (_seedingLayers[i + 1] != "*")
            // layer
            geoid = geoid.setLayer(std::stoi(_seedingLayers[i + 1]));

        geoSelection.push_back(geoid);
    }

    _seedGeometrySelection = MarlinACTS::GeometryIdSelector(geoSelection);

    if (_seedFinding_deltaRMinTop == 0.f) _seedFinding_deltaRMinTop = _seedFinding_deltaRMin;
    if (_seedFinding_deltaRMaxTop == 0.f) _seedFinding_deltaRMaxTop = _seedFinding_deltaRMax;
    if (_seedFinding_deltaRMinBottom == 0.f) _seedFinding_deltaRMinBottom = _seedFinding_deltaRMin;
    if (_seedFinding_deltaRMaxBottom == 0.f) _seedFinding_deltaRMaxBottom = _seedFinding_deltaRMax;

    if (_seedFinding_zBinSchema.size() > 0)
    {
        for (std::string token : _seedFinding_zBinSchema)
        {
            auto [ pos, ztop_sx, ztop_dx, zbottom_sx, zbottom_dx, err ] = MarlinACTS::parseZBinSchema(token);
            if (err || pos < -_seedFinding_zMax || pos >= _seedFinding_zMax)
            {
                streamlog_out(WARNING) << "Wrong parameter SeedFinding_zBinSchema; default used" << std::endl;
                _zBinEdges.clear();
                break;
            }
            _zBinEdges.push_back(pos);
            _ZTopBinSchema.emplace_back(-ztop_sx, ztop_dx);
            _ZBottomBinSchema.emplace_back(-zbottom_sx, zbottom_dx);
        }
    }

    if (_zBinEdges.empty())
    {
        // taken from Acts::CylindricalSpacePointGridCreator::createGrid
        float tmpf = 2 * _seedFinding_zMax / (_seedFinding_cotThetaMax * _seedFinding_deltaRMax);
        int num_bins = std::max(2, static_cast<int>(std::max(1.f, std::floor(tmpf))));

        for (int bin = 0; bin < num_bins; bin++)
        {
            _zBinEdges.push_back(-_seedFinding_zMax + bin * 2 * _seedFinding_zMax / num_bins);
            _ZTopBinSchema.emplace_back(-_zTopBinLen, _zTopBinLen);
            _ZBottomBinSchema.emplace_back(-_zBottomBinLen, _zBottomBinLen);
        }
    }

    finderCfg.rMax = _seedFinding_rMax;
    finderCfg.deltaRMin = _seedFinding_deltaRMin;
    finderCfg.deltaRMax = _seedFinding_deltaRMax;
    finderCfg.deltaRMinTopSP = _seedFinding_deltaRMinTop;
    finderCfg.deltaRMaxTopSP = _seedFinding_deltaRMaxTop;
    finderCfg.deltaRMinBottomSP = _seedFinding_deltaRMinBottom;
    finderCfg.deltaRMaxBottomSP = _seedFinding_deltaRMaxBottom;
    finderCfg.collisionRegionMin = -_seedFinding_collisionRegion;
    finderCfg.collisionRegionMax = _seedFinding_collisionRegion;
    finderCfg.zMin = -_seedFinding_zMax;
    finderCfg.zMax = _seedFinding_zMax;
    finderCfg.maxSeedsPerSpM = 1;
    finderCfg.cotThetaMax = 7.40627;  // 2.7 eta;
    finderCfg.sigmaScattering = _seedFinding_sigmaScattering;
    finderCfg.radLengthPerSeed = _seedFinding_radLengthPerSeed;
    finderCfg.minPt = _seedFinding_minPt * Acts::UnitConstants::MeV;
    finderCfg.impactMax = _seedFinding_impactMax * Acts::UnitConstants::mm;
    finderCfg.useVariableMiddleSPRange = true;

    Acts::SeedFilterConfig filterCfg;
    filterCfg.maxSeedsPerSpM = finderCfg.maxSeedsPerSpM;

    finderCfg.seedFilter = std::make_unique<Acts::SeedFilter<MarlinACTS::SeedSpacePoint>>(
        Acts::SeedFilter<MarlinACTS::SeedSpacePoint>(filterCfg.toInternalUnits()));
    finderCfg = finderCfg.toInternalUnits().calculateDerivedQuantities();

    gridCfg.cotThetaMax = finderCfg.cotThetaMax;
    gridCfg.deltaRMax = finderCfg.deltaRMax;
    gridCfg.minPt = finderCfg.minPt;
    gridCfg.rMax = finderCfg.rMax;
    gridCfg.zMax = finderCfg.zMax;
    gridCfg.zMin = finderCfg.zMin;
    gridCfg.impactMax = finderCfg.impactMax;

    gridCfg.zBinEdges.resize(_zBinEdges.size());
    for (int k = 0; k < _zBinEdges.size(); k++) gridCfg.zBinEdges[k] = _zBinEdges[k];

}

ACTSCKFBaseTracker::SeedParamList ACTSCKFSeededTracker::getSeeds(const MarlinACTS::MeasurementContainer& m_list)
{
    /* ********************************************************************************************
     *  Space points
     ******************************************************************************************* */

    MarlinACTS::SeedSpacePointContainer spacePoints;

    for (auto meas : m_list)
    {
        const MarlinACTS::SourceLink& s_link = meas.sourceLink().get<MarlinACTS::SourceLink>();
        if (!_seedGeometrySelection.check(s_link.geometryId()))
        {
            std::cout << "Check failure for geometry ID " << s_link.geometryId() << std::endl;
            continue;
        }

        const EVENT::TrackerHit* lciohit = s_link.lciohit();
        const Acts::Surface *surface = findSurface(lciohit);

        const double* lcioglobalpos = lciohit->getPosition();
        Acts::Vector3 globalPos = { lcioglobalpos[0], lcioglobalpos[1], lcioglobalpos[2] };
        Acts::RotationMatrix3 rotLocalToGlobal = surface->referenceFrame(geometryContext(), globalPos, {0, 0, 0});

        // Convert to a seed space point
        // the space point requires only the variance of the transverse and
        // longitudinal position. reduce computations by transforming the
        // covariance directly from local to rho/z.
        //
        // compute Jacobian from global coordinates to rho/z
        //
        //         rho = sqrt(x² + y²)
        // drho/d{x,y} = (1 / sqrt(x² + y²)) * 2 * {x,y}
        //             = 2 * {x,y} / r
        //       dz/dz = 1 (duuh!)
        //
        double x = globalPos[Acts::ePos0];
        double y = globalPos[Acts::ePos1];
        double scale = 2 / std::hypot(x, y);

        Acts::ActsMatrix<2, 3> jacXyzToRhoZ = Acts::ActsMatrix<2, 3>::Zero();
        jacXyzToRhoZ(0, Acts::ePos0) = scale * x;
        jacXyzToRhoZ(0, Acts::ePos1) = scale * y;
        jacXyzToRhoZ(1, Acts::ePos2) = 1;

        // compute Jacobian from local coordinates to rho/z
        Acts::ActsMatrix<2, 2> jac =
          jacXyzToRhoZ * rotLocalToGlobal.block<3, 2>(Acts::ePos0, Acts::ePos0);
        // compute rho/z variance
        Acts::ActsVector<2> var = (jac * meas.covariance() * jac.transpose()).diagonal();  // TODO check variance

        spacePoints.push_back(MarlinACTS::SeedSpacePoint(globalPos, var[0], var[1], s_link));

    }

    /* ********************************************************************************************
     *  Seed finder setup
     ******************************************************************************************* */

    const Acts::Vector3 zeropos(0, 0, 0);
    double B_zero = magneticFieldValue(zeropos)[2];

    Acts::SeedFinderOptions finderOpts;
    finderOpts.bFieldInZ = B_zero;
    finderOpts.beamPos = { 0, 0 };
    finderOpts = finderOpts.toInternalUnits();
    finderOpts = finderOpts.calculateDerivedQuantities(finderCfg);

    Acts::CylindricalSpacePointGridOptions gridOpts;
    gridOpts.bFieldInZ = B_zero;

    auto extractGlobalQuantities = [](const SSPoint& sp, float, float, float) {
        Acts::Vector3 position { sp.x(), sp.y(), sp.z() };
        Acts::Vector2 covariance { sp.varianceR(), sp.varianceZ() };
        return std::make_tuple(position, covariance, sp.t());
    };

    std::vector<const MarlinACTS::SeedSpacePoint *> spacePointPtrs(
        spacePoints.size(), nullptr);
    std::transform(spacePoints.begin(), spacePoints.end(), spacePointPtrs.begin(),
                 [](const MarlinACTS::SeedSpacePoint &sp) { return &sp; });

    Acts::Extent rRangeSPExtent;

    SSPointGrid grid = Acts::CylindricalSpacePointGridCreator::createGrid<SSPoint>(
        gridCfg.toInternalUnits(), gridOpts.toInternalUnits());

    Acts::CylindricalSpacePointGridCreator::fillGrid(finderCfg, finderOpts, grid,
        spacePointPtrs.begin(), spacePointPtrs.end(), extractGlobalQuantities,
        rRangeSPExtent);

    const Acts::GridBinFinder<2ul> bottomBinFinder(_phiBottomBinLen, _ZBottomBinSchema);
    const Acts::GridBinFinder<2ul> topBinFinder(_phiTopBinLen, _ZTopBinSchema);

    auto spacePointsGrouping = Acts::CylindricalBinnedGroup<SSPoint>(std::move(grid),
        bottomBinFinder, topBinFinder);

    Acts::SeedFinder<SSPoint, SSPointGrid> finder(finderCfg);
    decltype(finder)::SeedingState state;

    state.spacePointData.resize(spacePointPtrs.size(),
        finderCfg.useDetailedDoubleMeasurementInfo);

    float up = Acts::clampValue<float>(
      std::floor(rRangeSPExtent.max(Acts::BinningValue::binR) / 2) * 2);
    const Acts::Range1D<float> rMiddleSPRange(
      std::floor(rRangeSPExtent.min(Acts::BinningValue::binR) / 2) * 2 +
          finderCfg.deltaRMiddleMinSPRange,
      up - finderCfg.deltaRMiddleMaxSPRange);

    /* ********************************************************************************************
     *  Seeding
     ******************************************************************************************* */

    SeedParamList paramseeds;

    for (const auto [bottom, middle, top] : spacePointsGrouping)
    {
        std::vector<Acts::Seed<SSPoint>> seeds;
        finder.createSeedsForGroup(finderOpts, state, spacePointsGrouping.grid(),
                                   seeds, bottom, middle, top, rMiddleSPRange);

        for (const Acts::Seed<SSPoint> &seed : seeds)
        {
            const SSPoint* bottomSP = seed.sp().front();

            const auto& sourceLink = bottomSP->sourceLink();
            const Acts::GeometryIdentifier& geoId = sourceLink.geometryId();
            const Acts::Surface* surface = trackingGeometry()->findSurface(geoId);
            if (surface == nullptr)
            {
                std::cout << "surface with geoID " << geoId
                          << " is not found in the tracking geometry";
                continue;
            }

            const Acts::Vector3 seedPos(bottomSP->x(), bottomSP->y(), bottomSP->z());
            Acts::Vector3 seedField = magneticFieldValue(seedPos);

            std::optional<Acts::BoundVector> optParams =
              Acts::estimateTrackParamsFromSeed(geometryContext(), seed.sp().begin(), seed.sp().end(),
                                                *surface, seedField, 0.1 * Acts::UnitConstants::T);
            if (!optParams.has_value())
            {
                std::cout << "Failed estimation of track parameters for seed." << std::endl;
                continue;
            }

            const Acts::BoundVector &params = optParams.value();

            //float charge = std::copysign(1, params[Acts::eBoundQOverP]);
            float p = std::abs(1 / params[Acts::eBoundQOverP]);

            Acts::BoundSquareMatrix cov = Acts::BoundSquareMatrix::Zero();
            cov(Acts::eBoundLoc0, Acts::eBoundLoc0) = std::pow(_initialTrackError_d0, 2);
            cov(Acts::eBoundLoc1, Acts::eBoundLoc1) = std::pow(_initialTrackError_z0, 2);
            cov(Acts::eBoundTime, Acts::eBoundTime) = std::pow(_initialTrackError_time, 2);
            cov(Acts::eBoundPhi, Acts::eBoundPhi) = std::pow(_initialTrackError_phi, 2);
            cov(Acts::eBoundTheta, Acts::eBoundTheta) = std::pow(_initialTrackError_lambda, 2);
            cov(Acts::eBoundQOverP, Acts::eBoundQOverP) = std::pow(_initialTrackError_relP * p / (p * p), 2);

            Acts::BoundTrackParameters paramseed(surface->getSharedPtr(), params,
                                                 cov, Acts::ParticleHypothesis::pion());
            paramseeds.push_back(paramseed);
        }
    }

    return paramseeds;
}
