#include "ACTSCKFBaseTracker.h"

#include <Acts/Definitions/Units.hpp>
#include <Acts/TrackFinding/MeasurementSelector.hpp>
#include "Acts/TrackFinding/TrackStateCreator.hpp"

#include <IMPL/TrackImpl.h>
#include <IMPL/LCFlagImpl.h>
#include <IMPL/TrackerHitPlaneImpl.h>

#include "Measurement.h"
#include "MeasurementCalibrator.h"

#include <cmath>

ACTSCKFBaseTracker::ACTSCKFBaseTracker(const string& procname) :
    ACTSBaseTracker(procname)
{
    registerProcessorParameter("CKF_Chi2CutOff", "Maximum local chi2 contribution.",
                               _CKF_chi2CutOff, _CKF_chi2CutOff);

    registerProcessorParameter("CKF_NumMeasurementsCutOff",
                               "Maximum number of associated measurements on a single surface.",
                               _CKF_numMeasurementsCutOff, _CKF_numMeasurementsCutOff);

    registerProcessorParameter("InitialTrackError_RelP", 
                               "Track error estimate, momentum component (relative).",
                               _initialTrackError_relP, 0.25);

    registerProcessorParameter("InitialTrackError_Phi", "Track error estimate, phi (degree).",
                               _initialTrackError_phi, 1.0);

    registerProcessorParameter("InitialTrackError_Lambda", "Track error estimate, lambda (degree).",
                             _initialTrackError_lambda, 1.0);

    registerProcessorParameter("InitialTrackError_D0", "Track error estimate, local position for D0 (mm).",
                               _initialTrackError_d0, 0.01);

    registerProcessorParameter("InitialTrackError_Z0", "Track error estimate, local position for Z0 (mm).",
                               _initialTrackError_z0, 0.01);

    registerProcessorParameter("InitialTrackError_Time ", "Track error estimate, time (ns).",
                               _initialTrackError_time, 1.0);

    registerProcessorParameter("ThetaTolerance", "Tolerance for theta in rad.",
                               theta_tolerance, 0.01f * static_cast<float>(M_PI));

    registerProcessorParameter("MTCKFMode", "Enable MT mode for CKF", mtckf_mode, false);
}

void ACTSCKFBaseTracker::init()
{
    ACTSBaseTracker::init();

    _initialTrackError_phi *= Acts::UnitConstants::degree;
    _initialTrackError_lambda *= Acts::UnitConstants::degree;
    _initialTrackError_d0 *= Acts::UnitConstants::mm;
    _initialTrackError_z0 *= Acts::UnitConstants::mm;
    _initialTrackError_time *= Acts::UnitConstants::ns;

    trackFinder.reset(new CKF(*propagator));
}

void ACTSCKFBaseTracker::processEvent(LCEvent* evt)
{
    ACTSBaseTracker::processEvent(evt);

    vector<std::pair<Acts::GeometryIdentifier, EVENT::TrackerHit*>> sortedHits;

    for (const std::string& collection : _inputTrackerHitCollections)
    {
        LCCollection* trackerHitCollection = getCollection(collection, evt);
        if (trackerHitCollection == nullptr) continue;

        for (int idxHit = 0; idxHit < trackerHitCollection->getNumberOfElements(); idxHit++)
        {
            EVENT::TrackerHit* hit = static_cast<EVENT::TrackerHit*>(
                trackerHitCollection->getElementAt(idxHit));

            sortedHits.push_back(std::make_pair(geoIDMappingTool()->getGeometryID(hit), hit));
        }
    }

    // Sort by GeoID
    std::sort( sortedHits.begin(), sortedHits.end(),
        [](const std::pair<Acts::GeometryIdentifier, EVENT::TrackerHit*>& hit0,
           const std::pair<Acts::GeometryIdentifier, EVENT::TrackerHit*>& hit1)
           -> bool { return hit0.first < hit1.first; });

    /* ********************************************************************************************
     *  Measurements setup
     ******************************************************************************************* */

    MarlinACTS::SourceLinkContainer sourceLinks;
    MarlinACTS::MeasurementContainer measurements;

    sourceLinks.reserve(sortedHits.size());
    for (std::pair<Acts::GeometryIdentifier, EVENT::TrackerHit*>& hit : sortedHits)
    {
        const Acts::Surface* surface = trackingGeometry()->findSurface(hit.first);
        if (surface == nullptr) throw std::runtime_error("Surface not found");

        const double* lcioglobalpos = hit.second->getPosition();
        Acts::Vector3 globalPos = {
            lcioglobalpos[0], lcioglobalpos[1], lcioglobalpos[2] };

        Acts::Result<Acts::Vector2> lpResult =
            surface->globalToLocal(geometryContext(), globalPos,
                                   {0, 0, 0}, 0.5 * Acts::UnitConstants::um);
        if (!lpResult.ok())
          throw std::runtime_error("Global to local transformation did not succeed.");

        Acts::Vector2 loc = lpResult.value();

        Acts::SquareMatrix2 localCov = Acts::SquareMatrix2::Zero();
        const EVENT::TrackerHitPlane* hitplane =
            dynamic_cast<const EVENT::TrackerHitPlane*>(hit.second);
        if (hitplane) {
            localCov(0, 0) = std::pow(hitplane->getdU() * Acts::UnitConstants::mm, 2);
            localCov(1, 1) = std::pow(hitplane->getdV() * Acts::UnitConstants::mm, 2);
        } else {
            throw std::runtime_error("Currently only support TrackerHitPlane.");
        }

        MarlinACTS::SourceLink sourceLink(surface->geometryId(),
                                          measurements.size(), hit.second);
        Acts::SourceLink src_wrap { sourceLink };
        MarlinACTS::Measurement meas = MarlinACTS::makeMeasurement(
            src_wrap, loc, localCov, Acts::eBoundLoc0, Acts::eBoundLoc1);

        measurements.push_back(meas);
        sourceLinks.emplace_hint(sourceLinks.end(), sourceLink);
    }

    Acts::MeasurementSelector::Config measurementSelectorCfg = {
        { Acts::GeometryIdentifier(),
          { {}, { _CKF_chi2CutOff }, { (std::size_t)(_CKF_numMeasurementsCutOff) } } }
    };

    Acts::MeasurementSelector measSel { measurementSelectorCfg };
    MarlinACTS::MeasurementCalibrator measCal { measurements };

    MarlinACTS::SourceLinkAccessor slAccessor;
    slAccessor.container = &sourceLinks;

    /* ********************************************************************************************
     *  Seed setup
     ******************************************************************************************* */

    auto seeds = getSeeds(measurements, evt);

    /* ********************************************************************************************
     *  CKF setup
     ******************************************************************************************* */

    Acts::PropagatorPlainOptions pOptions { geometryContext(), magneticFieldContext() };
    pOptions.maxSteps = 10000;

    Acts::GainMatrixUpdater kfUpdater;

    using TrackStateCreatorType = Acts::TrackStateCreator<MarlinACTS::SourceLinkAccessor::Iterator,
            TrackContainer>;
    TrackStateCreatorType trackStateCreator;
    trackStateCreator.sourceLinkAccessor
        .template connect<&MarlinACTS::SourceLinkAccessor::range>(&slAccessor);
    trackStateCreator.calibrator
        .template connect<&MarlinACTS::MeasurementCalibrator::calibrate>(&measCal);
    trackStateCreator.measurementSelector
        .template connect<&Acts::MeasurementSelector::select<Acts::VectorMultiTrajectory>>(&measSel);

    Acts::CombinatorialKalmanFilterExtensions<TrackContainer> extensions;
    extensions.updater.connect<
        &Acts::GainMatrixUpdater::operator()<Acts::VectorMultiTrajectory>>(&kfUpdater);
    extensions.createTrackStates
          .template connect<&TrackStateCreatorType::createTrackStates>( &trackStateCreator);

    TrackFinderOptions ckfOptions = TrackFinderOptions(geometryContext(), magneticFieldContext(),
        calibrationContext(), extensions, pOptions);

    /* ********************************************************************************************
     *  Track finding
     ******************************************************************************************* */

    MarlinACTS::CollectionHandler coll_handler { _outputTrackCollection, true, theta_tolerance };

    ExtrapolationOpts extrapolationOptions(geometryContext(), magneticFieldContext());

    if (mtckf_mode)
    {
        tbb::parallel_for(tbb::blocked_range<std::size_t>(0, seeds.size()),
                          ACTSTrackFinderThread(*this, seeds, coll_handler,
                                                extrapolationOptions, ckfOptions));
    }
    else
    {
        for (std::size_t iseed = 0; iseed < seeds.size(); ++iseed)
        {
            try
            {
                runCKFonSeed(seeds.at(iseed), ckfOptions, extrapolationOptions, coll_handler);
            }
            catch (std::error_code error)
            {
                streamlog_out(DEBUG) << "Track finding failure: " << error << std::endl;
            }

        }
    }

    coll_handler.process();

    coll_handler.saveCollection(evt);
}


/* ********************************************************************************************
 * The common track finder
 ******************************************************************************************* */
void ACTSCKFBaseTracker::runCKFonSeed(Acts::BoundTrackParameters& seed,
                                      TrackFinderOptions& ckfOptions,
                                      ExtrapolationOpts& extrapolationOptions,
                                      MarlinACTS::CollectionHandler& coll_handler)
{
    auto trackContainer = std::make_shared<Acts::VectorTrackContainer>();
    auto trackStateContainer = std::make_shared<Acts::VectorMultiTrajectory>();
    TrackContainer tracks(trackContainer, trackStateContainer);

    auto result = trackFinder->findTracks(seed, ckfOptions, tracks);
    if (!result.ok())
    {
        if (mtckf_mode)
        {
            tbb::queuing_mutex::scoped_lock lock(p_mutex);
            _fitFails++;
        }
        else _fitFails++;
        throw result.error();
    }

    const auto& fitOutput = result.value();
    for (const auto& trackItem : fitOutput)
    {
        auto trackTip = tracks.makeTrack();
        trackTip.copyFrom(trackItem, true);
        auto smoothResult = Acts::smoothTrack(geometryContext(), trackTip);
        if (!smoothResult.ok())
        {
            if (mtckf_mode)
            {
                tbb::queuing_mutex::scoped_lock lock(p_mutex);
                _fitFails++;
            }
            else _fitFails++;
            throw smoothResult.error();
        }

        auto extrapolationResult = Acts::extrapolateTrackToReferenceSurface(
            trackTip, *perigeeSurface, *propagator, extrapolationOptions,
            Acts::TrackExtrapolationStrategy::firstOrLast);
        if (!extrapolationResult.ok())
        {
            if (mtckf_mode)
            {
                tbb::queuing_mutex::scoped_lock lock(p_mutex);
                _fitFails++;
            }
            else _fitFails++;
            throw extrapolationResult.error();
        }

        if (mtckf_mode)
        {
            tbb::queuing_mutex::scoped_lock lock(p_mutex);
            coll_handler.addTrack(convert_track(trackTip));
        }
        else coll_handler.addTrack(convert_track(trackTip));
    }
}

/* ********************************************************************************************
 * Functor for MT CKF
 ******************************************************************************************* */
void ACTSTrackFinderThread::operator()
    (const tbb::blocked_range<std::size_t>& prange) const
{
    for (std::size_t iseed = prange.begin(); iseed != prange.end(); ++iseed)
    {
        try
        {
            processor.runCKFonSeed(paramseeds.at(iseed), ckfFinderOptions,
                                   extrapolOptions, trackCollection);
        }
        catch (std::error_code error)
        {
            // TODO synchronized log
        }
    }
}
