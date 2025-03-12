#include "ACTSCKFBaseTracker.h"

#include <Acts/Definitions/Units.hpp>
#include <Acts/TrackFinding/MeasurementSelector.hpp>

#include <IMPL/TrackImpl.h>
#include <IMPL/LCFlagImpl.h>
#include <IMPL/TrackerHitPlaneImpl.h>

#include "CollectionHandler.h"
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

    registerProcessorParameter("InitialTrackError_Phi", "Track error estimate, phi (radians).",
                               _initialTrackError_phi, 1 * Acts::UnitConstants::degree);

    registerProcessorParameter("InitialTrackError_Lambda", "Track error estimate, lambda (radians).",
                             _initialTrackError_lambda, 1 * Acts::UnitConstants::degree);

    registerProcessorParameter("InitialTrackError_D0", "Track error estimate, local position for D0 (mm).",
                               _initialTrackError_d0, 10 * Acts::UnitConstants::um);

    registerProcessorParameter("InitialTrackError_Z0", "Track error estimate, local position for Z0 (mm).",
                               _initialTrackError_z0, 10 * Acts::UnitConstants::um);

    registerProcessorParameter("InitialTrackError_Time ", "Track error estimate, time (ns).",
                               _initialTrackError_z0, 1 * Acts::UnitConstants::ns);

    registerProcessorParameter("ThetaTolerance", "Tolerance for theta in rad.",
                               theta_tolerance, 0.01f * static_cast<float>(M_PI));
}

void ACTSCKFBaseTracker::init()
{
    ACTSBaseTracker::init();

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

    Acts::SourceLinkAccessorDelegate<MarlinACTS::SourceLinkAccessor::Iterator> slAccessorDelegate;
    slAccessorDelegate.connect<&MarlinACTS::SourceLinkAccessor::range>(&slAccessor);

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

    Acts::CombinatorialKalmanFilterExtensions<TrackContainer> extensions;
    extensions.calibrator.connect<&MarlinACTS::MeasurementCalibrator::calibrate>(&measCal);
    extensions.updater.connect<
        &Acts::GainMatrixUpdater::operator()<Acts::VectorMultiTrajectory>>(&kfUpdater);
    extensions.measurementSelector.connect<
        &Acts::MeasurementSelector::select<Acts::VectorMultiTrajectory>>(&measSel);


    TrackFinderOptions ckfOptions = TrackFinderOptions(geometryContext(), magneticFieldContext(),
        calibrationContext(), slAccessorDelegate, extensions, pOptions);

    /* ********************************************************************************************
     *  Track finding
     ******************************************************************************************* */

    MarlinACTS::CollectionHandler coll_handler { _outputTrackCollection, true, theta_tolerance };

    auto trackContainer = std::make_shared<Acts::VectorTrackContainer>();
    auto trackStateContainer = std::make_shared<Acts::VectorMultiTrajectory>();
    TrackContainer tracks(trackContainer, trackStateContainer);

    Propagator::template Options<Acts::ActorList<Acts::MaterialInteractor, Acts::EndOfWorldReached>>
        extrapolationOptions(geometryContext(), magneticFieldContext());

    for (std::size_t iseed = 0; iseed < seeds.size(); ++iseed)
    {
        tracks.clear();

        auto result = trackFinder->findTracks(seeds.at(iseed), ckfOptions, tracks);
        if (!result.ok())
        {
            streamlog_out(WARNING) << "Track fit error: " << result.error() << std::endl;
            _fitFails++;
            continue;
        }

        const auto& fitOutput = result.value();
        for (const auto& trackItem : fitOutput)
        {
            auto trackTip = tracks.makeTrack();
            trackTip.copyFrom(trackItem, true);
            auto smoothResult = Acts::smoothTrack(geometryContext(), trackTip);
            if (!smoothResult.ok())
            {
                streamlog_out(DEBUG) << "Smooth failure: " << smoothResult.error() << std::endl;
                continue;
            }

            auto extrapolationResult = Acts::extrapolateTrackToReferenceSurface(
                trackTip, *perigeeSurface, *propagator, extrapolationOptions,
                Acts::TrackExtrapolationStrategy::firstOrLast);
            if (!extrapolationResult.ok())
            {
                streamlog_out(DEBUG) << "Extrapolation failure: " 
                                     << extrapolationResult.error() << std::endl;
                continue;
            }

            coll_handler.addTrack(convert_track(trackTip));
        }
    }

    coll_handler.saveCollection(evt);
}
