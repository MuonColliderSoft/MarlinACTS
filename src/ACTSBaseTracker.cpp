#include "ACTSBaseTracker.h"

#include <DD4hep/DD4hepUnits.h>
#include <DD4hep/Detector.h>

#include <UTIL/LCTrackerConf.h>
#include <IMPL/LCFlagImpl.h>
#include <IMPL/TrackImpl.h>
#include <IMPL/TrackStateImpl.h>
#include <UTIL/CellIDDecoder.h>

#include <Acts/Definitions/Units.hpp>
#include <Acts/Geometry/CylinderVolumeBuilder.hpp>
#include <Acts/Geometry/CylinderVolumeHelper.hpp>
#include <Acts/Geometry/ITrackingVolumeBuilder.hpp>
#include <Acts/Geometry/LayerArrayCreator.hpp>
#include <Acts/Geometry/LayerCreator.hpp>
#include <Acts/Geometry/ProtoLayerHelper.hpp>
#include <Acts/Geometry/SurfaceArrayCreator.hpp>
#include <Acts/Geometry/TrackingGeometryBuilder.hpp>
#include <Acts/Geometry/TrackingVolumeArrayCreator.hpp>
#include <Acts/MagneticField/ConstantBField.hpp>
#include <Acts/Plugins/DD4hep/DD4hepFieldAdapter.hpp>
#include <Acts/Plugins/Json/JsonMaterialDecorator.hpp>
#include <Acts/Plugins/TGeo/TGeoDetectorElement.hpp>
#include <Acts/Plugins/TGeo/TGeoLayerBuilder.hpp>
#include "Acts/Propagator/MaterialInteractor.hpp"
#include <Acts/Utilities/BinningType.hpp>

#include "Helpers.h"
#include "SourceLink.h"

#include <string>

using MarlinACTS::GeometryIdMappingTool;
using DetSchema = GeometryIdMappingTool::DetSchema;

ACTSBaseTracker::ACTSBaseTracker(const string& procname) :
    Processor(procname)
{
    registerProcessorParameter("MatFile", "Path to the material description JSON file. Can be empty.",
                               _matFile, _matFile);

    registerProcessorParameter("TGeoFile", "Path to the tracker geometry file.",
                               _tgeoFile, _tgeoFile);

    registerProcessorParameter("TGeoDescFile", "Path to the JSON file describing the subdetectors.",
                               _tgeodescFile, _tgeodescFile);

    registerProcessorParameter("DetectorSchema",
                               "Detector schema name (MuColl_v1, MAIA_v0, MuSIC_v1, MuSIC_v2).",
                               _detSchema, _detSchema);

    registerInputCollections(LCIO::TRACKERHITPLANE, "TrackerHitCollectionNames",
                             "Name of the TrackerHit input collections.",_inputTrackerHitCollections, {});

    registerOutputCollection(LCIO::TRACK, "TrackCollectionName", "Name of track output collection.",
                             _outputTrackCollection, string("Tracks"));

    registerProcessorParameter("CaloFace_Radius", "ECAL Inner Radius (mm).", _caloFaceR, 1857.f);

    registerProcessorParameter("CaloFace_Z", "ECAL half length (mm).", _caloFaceZ, 2307.f);

}

const Acts::Surface* ACTSBaseTracker::findSurface(const EVENT::TrackerHit* hit) const
{
    uint64_t moduleGeoId = _geoIDMappingTool->getGeometryID(hit);
    return _trackingGeometry->findSurface(moduleGeoId);
}

/* ********************************************************************************************
 *
 * Processor initialization
 *
 * ******************************************************************************************** */
void ACTSBaseTracker::init()
{
    _matFile = MarlinACTS::findFile(_matFile);
    _tgeoFile = MarlinACTS::findFile(_tgeoFile);
    _tgeodescFile = MarlinACTS::findFile(_tgeodescFile);

    _runNumber = 0;
    _eventNumber = 0;
    _fitFails = 0;

    printParameters();

    streamlog_out(MESSAGE) << " -------------------------------------" << std::endl;

    streamlog_out(MESSAGE) << " -- Building magnetic field" << std::endl;
    buildBfield();
    streamlog_out(MESSAGE) << " -- Building tracking detector" << std::endl;
    buildDetector();

    streamlog_out(MESSAGE) << " -------------------------------------" << std::endl;

    _magCache = _magneticField->makeCache(_magneticFieldContext);

    DetSchema dSchema;
    if (_detSchema == "MuSIC_v1") dSchema = DetSchema::MuSIC_v1;
    else if (_detSchema == "MuSIC_v2") dSchema = DetSchema::MuSIC_v2;
    else if (_detSchema == "MAIA_v0") dSchema = DetSchema::MAIA_v0;
    else if (_detSchema == "MuColl_v1") dSchema = DetSchema::MuColl_v1;
    else dSchema = DetSchema::MuColl_v1;

    _geoIDMappingTool = std::make_shared<GeometryIdMappingTool>(
        lcio::LCTrackerCellID::encoding_string(), dSchema);

    Navigator::Config navigatorCfg { trackingGeometry() };
    navigatorCfg.resolvePassive = false;
    navigatorCfg.resolveMaterial = true;
    navigatorCfg.resolveSensitive = true;

    Stepper stepper(magneticField());
    Navigator navigator(navigatorCfg);

    propagator.reset(new Propagator(std::move(stepper), std::move(navigator)));

    perigeeSurface = Acts::Surface::makeShared<Acts::PerigeeSurface>(Acts::Vector3 {0., 0., 0.});

    auto caloCylinder = std::make_shared<Acts::CylinderBounds>(_caloFaceR, _caloFaceZ);
    caloSurface = Acts::Surface::makeShared<Acts::CylinderSurface>(Acts::Transform3::Identity(), caloCylinder);

    Acts::Translation3 circlePositionL(0, 0, -_caloFaceZ);
    Acts::Translation3 circlePositionR(0, 0, _caloFaceZ);
    caloCapL = Acts::Surface::makeShared<Acts::DiscSurface>(Acts::Transform3(circlePositionL), 0. , _caloFaceR);
    caloCapR = Acts::Surface::makeShared<Acts::DiscSurface>(Acts::Transform3(circlePositionR), 0. , _caloFaceR);
    
    trackCollection.reset(new LCCollectionVec(LCIO::TRACK));
    LCFlagImpl trkFlag(0);
    trkFlag.setBit(LCIO::TRBIT_HITS);
    trackCollection->setFlag(trkFlag.getFlag());
}

/* ********************************************************************************************
 *
 * Detector geometry conversion
 *
 * ******************************************************************************************** */
void ACTSBaseTracker::buildDetector()
{
    Acts::Logging::Level surfaceLogLevel = Acts::Logging::INFO;
    Acts::Logging::Level layerLogLevel = Acts::Logging::INFO;
    Acts::Logging::Level volumeLogLevel = Acts::Logging::INFO;

    // Material description
    std::shared_ptr<const Acts::IMaterialDecorator> matDeco = nullptr;
    if (!_matFile.empty())
    {
        Acts::MaterialMapJsonConverter::Config jsonGeoConvConfig;
        matDeco = std::make_shared<const Acts::JsonMaterialDecorator>(
            jsonGeoConvConfig, _matFile, Acts::Logging::INFO);
    }

    // Geometry
    TGeoManager* gGeoManagerOld = nullptr;
    if (!_tgeoFile.empty())
    {
        // Save current geometry. This is needed by all the other Processors
        gGeoManagerOld = gGeoManager;
        gGeoManager = nullptr;  // prevents it from being deleted

        TGeoManager::Import(_tgeoFile.c_str());
    }

    // configure surface array creator
    Acts::SurfaceArrayCreator::Config sacConfig;
    auto surfaceArrayCreator = std::make_shared<const Acts::SurfaceArrayCreator>(
        sacConfig, Acts::getDefaultLogger("SurfaceArrayCreator", surfaceLogLevel));

    // configure the proto layer helper
    Acts::ProtoLayerHelper::Config plhConfig;
    auto protoLayerHelper = std::make_shared<const Acts::ProtoLayerHelper>(
        plhConfig, Acts::getDefaultLogger("ProtoLayerHelper", layerLogLevel));

    // configure the layer creator that uses the surface array creator
    Acts::LayerCreator::Config lcConfig;
    lcConfig.surfaceArrayCreator = surfaceArrayCreator;
    auto layerCreator = std::make_shared<const Acts::LayerCreator>(
        lcConfig, Acts::getDefaultLogger("LayerCreator", layerLogLevel));

    // configure the layer array creator
    Acts::LayerArrayCreator::Config lacConfig;
    auto layerArrayCreator = std::make_shared<const Acts::LayerArrayCreator>(
        lacConfig, Acts::getDefaultLogger("LayerArrayCreator", layerLogLevel));

    // tracking volume array creator
    Acts::TrackingVolumeArrayCreator::Config tvacConfig;
    auto tVolumeArrayCreator =
        std::make_shared<const Acts::TrackingVolumeArrayCreator>(tvacConfig,
            Acts::getDefaultLogger("TrackingVolumeArrayCreator", volumeLogLevel));

    // configure the cylinder volume helper
    Acts::CylinderVolumeHelper::Config cvhConfig;
    cvhConfig.layerArrayCreator = layerArrayCreator;
    cvhConfig.trackingVolumeArrayCreator = tVolumeArrayCreator;
    auto cylinderVolumeHelper = std::make_shared<const Acts::CylinderVolumeHelper>(
        cvhConfig, Acts::getDefaultLogger("CylinderVolumeHelper", volumeLogLevel));

    // list the volume builders
    std::list<std::shared_ptr<const Acts::ITrackingVolumeBuilder>> volumeBuilders;

    // Detector definition
    std::vector<Acts::TGeoLayerBuilder::Config> layerBuilderConfigs;

    if (_tgeodescFile.empty())
        throw std::runtime_error("Required geometry description file (TGeoDescFile) missing.");

    // Open the description
    nlohmann::json tgeodesc;
    std::ifstream tgeodescFile(_tgeodescFile, std::ifstream::in | std::ifstream::binary);
    if (!tgeodescFile.is_open())
        throw std::runtime_error("Unable to open TGeo description file: "+_tgeodescFile);

    tgeodescFile >> tgeodesc;

    // Helper parsing functions
    auto range_from_json = [](const nlohmann::json& jsonpair) -> std::pair<double, double>
    {
        return { jsonpair["lower"], jsonpair["upper"] };
    };

    // Loop over volumes to define sub-detectors
    for(const auto& volume : tgeodesc["Volumes"])
    {
        // Volume information
        Acts::TGeoLayerBuilder::Config layerBuilderConfig;
        layerBuilderConfig.configurationName = volume["geo-tgeo-volume-name"];
        layerBuilderConfig.unit = 1 * Acts::UnitConstants::cm;
        layerBuilderConfig.autoSurfaceBinning = true;

        // AutoBinning
        std::vector<std::pair<double, double>> binTolerances { Acts::numBinningValues(), { 0., 0. } };
        binTolerances[static_cast<int>(Acts::BinningValue::binR)] = 
            range_from_json(volume["geo-tgeo-sfbin-r-tolerance"]);
        binTolerances[static_cast<int>(Acts::BinningValue::binZ)] =
            range_from_json(volume["geo-tgeo-sfbin-z-tolerance"]);
        binTolerances[static_cast<int>(Acts::BinningValue::binPhi)] =
            range_from_json(volume["geo-tgeo-sfbin-phi-tolerance"]);
        layerBuilderConfig.surfaceBinMatcher = Acts::SurfaceBinningMatcher(binTolerances);

        // Loop over subvolumes (two endcaps and one barrel)
        // List of possible subvolume names. Order corresponds to layerConfigurations.
        std::array<string, 3> subvolumeNames = {"negative","central","positive"};
        for(std::size_t idx = 0; idx < 3; idx++) {
            const string& subvolumeName = subvolumeNames[idx];
            if (!volume["geo-tgeo-volume-layers"][subvolumeName]) continue;

            // Create the layer config object and fill it
            Acts::TGeoLayerBuilder::LayerConfig lConfig;
            lConfig.volumeName = volume["geo-tgeo-subvolume-names"][subvolumeName];
            lConfig.sensorNames = volume["geo-tgeo-sensitive-names"][subvolumeName];
            lConfig.localAxes = volume["geo-tgeo-sensitive-axes"][subvolumeName];
            lConfig.envelope = std::pair<double, double>(0.1 * Acts::UnitConstants::mm, 
                                                         0.1 * Acts::UnitConstants::mm);

            // Fill the parsing restrictions in r
            lConfig.parseRanges.push_back({Acts::BinningValue::binR,
                range_from_json(volume["geo-tgeo-layer-r-ranges"][subvolumeName])});

            // Fill the parsing restrictions in z
            lConfig.parseRanges.push_back({Acts::BinningValue::binZ,
                range_from_json(volume["geo-tgeo-layer-z-ranges"][subvolumeName])});

            // Fill the layer splitting parameters in z
            float rsplit = volume["geo-tgeo-layer-r-split"][subvolumeName];
            if(rsplit > 0) {
            lConfig.splitConfigs.push_back({Acts::BinningValue::binR, rsplit});
            }

            // Fill the layer splitting parameters in z
            float zsplit = volume["geo-tgeo-layer-z-split"][subvolumeName];
            if(zsplit > 0) {
            lConfig.splitConfigs.push_back({Acts::BinningValue::binZ, zsplit});
            }

            layerBuilderConfig.layerConfigurations[idx].push_back(lConfig);
        }

        layerBuilderConfigs.push_back(layerBuilderConfig);
    }

    // remember the layer builders to collect the detector elements
    std::vector<std::shared_ptr<const Acts::TGeoLayerBuilder>> tgLayerBuilders;

    for (auto& lbc : layerBuilderConfigs)
    {
        std::shared_ptr<const Acts::LayerCreator> layerCreatorLB = nullptr;

        if (lbc.autoSurfaceBinning)
        {
            // Configure surface array creator (optionally) per layer builder
            // (in order to configure them to work appropriately)
            Acts::SurfaceArrayCreator::Config sacConfigLB;
            sacConfigLB.surfaceMatcher = lbc.surfaceBinMatcher;
            auto surfaceArrayCreatorLB =
            std::make_shared<const Acts::SurfaceArrayCreator>(sacConfigLB,
                Acts::getDefaultLogger(lbc.configurationName + "SurfaceArrayCreator",
                                       surfaceLogLevel));

            // configure the layer creator that uses the surface array creator
            Acts::LayerCreator::Config lcConfigLB;
            lcConfigLB.surfaceArrayCreator = surfaceArrayCreatorLB;
            layerCreatorLB = std::make_shared<const Acts::LayerCreator>(lcConfigLB,
                Acts::getDefaultLogger(lbc.configurationName + "LayerCreator",
                                       layerLogLevel));
        }

        // Configure the proto layer helper
        Acts::ProtoLayerHelper::Config plhConfigLB;
        auto protoLayerHelperLB = std::make_shared<const Acts::ProtoLayerHelper>(plhConfigLB,
            Acts::getDefaultLogger(lbc.configurationName + "ProtoLayerHelper", layerLogLevel));

        //-------------------------------------------------------------------------------------
        lbc.layerCreator =
            (layerCreatorLB != nullptr) ? layerCreatorLB : layerCreator;
        lbc.protoLayerHelper =
            (protoLayerHelperLB != nullptr) ? protoLayerHelperLB : protoLayerHelper;

        auto layerBuilder = std::make_shared<const Acts::TGeoLayerBuilder>(lbc,
            Acts::getDefaultLogger(lbc.configurationName + "LayerBuilder",layerLogLevel));
        // remember the layer builder
        tgLayerBuilders.push_back(layerBuilder);

        // build the pixel volume
        Acts::CylinderVolumeBuilder::Config volumeConfig;
        volumeConfig.trackingVolumeHelper = cylinderVolumeHelper;
        volumeConfig.volumeName = lbc.configurationName;
        volumeConfig.buildToRadiusZero = (volumeBuilders.size() == 0);
        volumeConfig.layerEnvelopeR = {1. * Acts::UnitConstants::mm,
                                       5. * Acts::UnitConstants::mm};
        auto ringLayoutConfiguration =
            [&](const std::vector<Acts::TGeoLayerBuilder::LayerConfig>& lConfigs)
            -> void {
          for (const auto& lcfg : lConfigs) {
            for (const auto& scfg : lcfg.splitConfigs) {
              if (scfg.first == Acts::BinningValue::binR and scfg.second > 0.) {
                volumeConfig.ringTolerance =
                    std::max(volumeConfig.ringTolerance, scfg.second);
                volumeConfig.checkRingLayout = true;
              }
            }
          }
        };
        ringLayoutConfiguration(lbc.layerConfigurations[0]);
        ringLayoutConfiguration(lbc.layerConfigurations[2]);
        volumeConfig.layerBuilder = layerBuilder;
        auto volumeBuilder = std::make_shared<const Acts::CylinderVolumeBuilder>(
            volumeConfig,
            Acts::getDefaultLogger(lbc.configurationName + "VolumeBuilder",
                                   volumeLogLevel));
        volumeBuilders.push_back(volumeBuilder);
    }

    // create the tracking geometry
    Acts::TrackingGeometryBuilder::Config tgConfig;
    // Add the builders
    tgConfig.materialDecorator = matDeco;

    for (auto& vb : volumeBuilders) 
    {
        tgConfig.trackingVolumeBuilders.push_back(
            [=](const auto& gcontext, const auto& inner, const auto&) {
              return vb->trackingVolume(gcontext, inner);
            });
    }
    // Add the helper
    tgConfig.trackingVolumeHelper = cylinderVolumeHelper;
    auto cylinderGeometryBuilder =
        std::make_shared<const Acts::TrackingGeometryBuilder>(tgConfig,
            Acts::getDefaultLogger("TrackerGeometryBuilder", volumeLogLevel));
    // get the geometry
    _trackingGeometry = cylinderGeometryBuilder->trackingGeometry(_geometryContext);
    // collect the detector element store
    for (auto& lBuilder : tgLayerBuilders)
    {
        auto detElements = lBuilder->detectorElements();
        _detectorStore.insert(_detectorStore.begin(), detElements.begin(),
                              detElements.end());
    }

    // Restore old gGeoManager
    if (gGeoManagerOld != nullptr) gGeoManager = gGeoManagerOld;

}

/* ********************************************************************************************
 *
 * Magnetic field from DD4Hep
 *
 * ******************************************************************************************** */
void ACTSBaseTracker::buildBfield()
{
    dd4hep::Detector& lcdd = dd4hep::Detector::getInstance();
    _magneticField = std::make_shared<Acts::DD4hepFieldAdapter>(lcdd.field());
}

LCCollection* ACTSBaseTracker::getCollection(const string& collectionName, LCEvent* evt)
{
    try
    {
        return evt->getCollection(collectionName);
    }
    catch (DataNotAvailableException& e)
    {
        streamlog_out(DEBUG5) << "- cannot get collection. Collection "
                              << collectionName << " is unavailable" << std::endl;
    }
    return nullptr;
}

Acts::Vector3 ACTSBaseTracker::magneticFieldValue(const Acts::Vector3 position)
{
    auto result = _magneticField->getField(position, _magCache);
    if (!result.ok())
    {
        throw std::runtime_error("Field lookup error: " + result.error().value());
    }
    return *result;
}

void ACTSBaseTracker::end() {
    streamlog_out(MESSAGE) << " end()  " << name() << " processed "
                           << _eventNumber << " events in " << _runNumber
                           << " runs " << std::endl;
}

/* ********************************************************************************************
 *
 * Tracks conversion
 *
 * ******************************************************************************************** */
EVENT::Track* ACTSBaseTracker::convert_track(const TrackResult& fitter_res)
{
    IMPL::TrackImpl* track = new IMPL::TrackImpl;

    track->setChi2(fitter_res.chi2());
    track->setNdf(fitter_res.nDoF());
    track->setNholes(fitter_res.nHoles());

    /* ********************************************************************************************
     *  Track state at IP
     ******************************************************************************************* */
    const Acts::Vector3 zeroPos(0, 0, 0);
    const Acts::BoundVector& params = fitter_res.parameters();
    const Acts::BoundMatrix& covariance = fitter_res.covariance();
    EVENT::TrackState* trackStateAtIP = convert_state(EVENT::TrackState::AtIP,
                                                      params, covariance, zeroPos);
    track->trackStates().push_back(trackStateAtIP);

    /* ********************************************************************************************
     *  Track states in the middle
     ******************************************************************************************* */
    EVENT::TrackerHitVec hitsOnTrack;
    EVENT::TrackStateVec statesOnTrack;

    for (const auto& trk_state : fitter_res.trackStatesReversed())
    {
        if (!trk_state.hasUncalibratedSourceLink()) continue;

        auto sl = trk_state.getUncalibratedSourceLink()
                           .get<MarlinACTS::SourceLink>();
        EVENT::TrackerHit* curr_hit = sl.lciohit();
        hitsOnTrack.push_back(curr_hit);

        const Acts::Vector3 hitPos(curr_hit->getPosition()[0],
                                   curr_hit->getPosition()[1],
                                   curr_hit->getPosition()[2]);

        EVENT::TrackState* trackState = convert_state(EVENT::TrackState::AtOther,
                                                      trk_state.parameters(), trk_state.covariance(),
                                                      hitPos);
        statesOnTrack.push_back(trackState);
    }

    /* ********************************************************************************************
     *  Track state at calorimeter
     ******************************************************************************************* */

    double d0 = params[Acts::eBoundLoc0];
    double z0 = params[Acts::eBoundLoc1];
    double phi = params[Acts::eBoundPhi];
    double theta = params[Acts::eBoundTheta];
    Acts::Vector3 pos(d0 * cos(phi), d0 * sin(phi), z0);

    Acts::CurvilinearTrackParameters start(Acts::VectorHelpers::makeVector4(pos, params[Acts::eBoundTime]),
                                           phi, theta, params[Acts::eBoundQOverP],
                                           covariance, Acts::ParticleHypothesis::pion());

    Propagator::template Options<Acts::ActionList<Acts::MaterialInteractor>,
                                 Acts::AbortList<Acts::EndOfWorldReached>>
	    caloPropOptions(geometryContext(), magneticFieldContext());
    caloPropOptions.pathLimit = 20 * Acts::UnitConstants::m;               // TODO missing conf. parameter

    auto resultProp = propagator->propagate(start, *caloSurface, caloPropOptions);
    if (!resultProp.ok())
    {
        resultProp = propagator->propagate(start, 
                                           (theta > M_PI/2) ? *caloCapL : *caloCapR, caloPropOptions);
    }

    if (resultProp.ok())
    {
        auto end_parameters = resultProp.value().endParameters;
        const Acts::BoundMatrix& atCaloCovariance = *(end_parameters->covariance());

        EVENT::TrackState* trackStateAtCalo = convert_state(EVENT::TrackState::AtCalorimeter,
                                                            end_parameters->parameters(),
                                                            atCaloCovariance,
                                                            zeroPos);   //TODO same field as in IP, is it correct??
        track->trackStates().push_back(trackStateAtCalo);
    }
    else
    {
        streamlog_out(DEBUG) << "Failed propagation to calorimeter!" << std::endl;
    }

    /* ********************************************************************************************
     *  Track conversion
     ******************************************************************************************* */
    std::reverse(hitsOnTrack.begin(), hitsOnTrack.end());
    std::reverse(statesOnTrack.begin(), statesOnTrack.end());

    UTIL::CellIDDecoder<lcio::TrackerHit> decoder(
        lcio::LCTrackerCellID::encoding_string());
    EVENT::IntVec& subdetectorHitNumbers = track->subdetectorHitNumbers();
    subdetectorHitNumbers.resize(7, 0);

    for (EVENT::TrackerHit* hit : hitsOnTrack)
    {
        track->addHit(hit);

        uint32_t sysid = decoder(hit)["system"];
        if (subdetectorHitNumbers.size() <= sysid)
        {
            subdetectorHitNumbers.resize(sysid + 1, 0);
        }
        subdetectorHitNumbers[sysid]++;
    }

    if (statesOnTrack.size() > 0)
    {
        dynamic_cast<IMPL::TrackStateImpl*>(statesOnTrack.back())
        ->setLocation(EVENT::TrackState::AtLastHit);
        dynamic_cast<IMPL::TrackStateImpl*>(statesOnTrack.front())
        ->setLocation(EVENT::TrackState::AtFirstHit);
    }

    EVENT::TrackStateVec& myTrackStates = track->trackStates();
    myTrackStates.insert(myTrackStates.end(), statesOnTrack.begin(),
                   statesOnTrack.end());

    return track;
}

EVENT::TrackState* ACTSBaseTracker::convert_state(int location,
    const Acts::BoundVector& value, const Acts::BoundMatrix& cov,
    Acts::Vector3 mag_pos)
{
    double Bz = magneticFieldValue(mag_pos)[2] / Acts::UnitConstants::T;

    IMPL::TrackStateImpl* trackState = new IMPL::TrackStateImpl();
    trackState->setLocation(location);

    // Central values
    double d0 = value[Acts::eBoundLoc0];
    double z0 = value[Acts::eBoundLoc1];
    double phi = value[Acts::eBoundPhi];
    double theta = value[Acts::eBoundTheta];
    double qoverp = value[Acts::eBoundQOverP];

    double p = 1e3 / qoverp;
    double omega = (0.3 * Bz) / (p * std::sin(theta));
    double lambda = M_PI / 2 - theta;
    double tanlambda = std::tan(lambda);

    trackState->setPhi(phi);
    trackState->setTanLambda(tanlambda);
    trackState->setOmega(omega);
    trackState->setD0(d0);
    trackState->setZ0(z0);

    // Uncertainties (covariance matrix)
    Acts::ActsMatrix<6, 6> jac = Acts::ActsMatrix<6, 6>::Zero();

    jac(0, Acts::eBoundLoc0) = 1;
    jac(1, Acts::eBoundPhi) = 1;
    jac(2, Acts::eBoundTheta) = omega / std::tan(theta);
    jac(2, Acts::eBoundQOverP) = omega / qoverp;
    jac(3, Acts::eBoundLoc1) = 1;
    jac(4, Acts::eBoundTheta) = std::pow(1 / std::cos(lambda), 2);

    Acts::ActsMatrix<6, 6> trcov = (jac * cov * jac.transpose());

    EVENT::FloatVec lcioCov(15, 0);
    lcioCov[0] = trcov(0, 0);
    lcioCov[1] = trcov(0, 1);
    lcioCov[2] = trcov(1, 1);
    lcioCov[3] = trcov(0, 2);
    lcioCov[4] = trcov(1, 2);
    lcioCov[5] = trcov(2, 2);
    lcioCov[6] = trcov(0, 3);
    lcioCov[7] = trcov(1, 3);
    lcioCov[8] = trcov(2, 3);
    lcioCov[9] = trcov(3, 3);
    lcioCov[10] = trcov(0, 4);
    lcioCov[11] = trcov(1, 4);
    lcioCov[12] = trcov(2, 4);
    lcioCov[13] = trcov(3, 4);
    lcioCov[14] = trcov(4, 4);

    trackState->setCovMatrix(lcioCov);

    return trackState;
}

