#include "ACTSBaseTracker.h"

#include <DD4hep/DD4hepUnits.h>
#include <DD4hep/Detector.h>

#include <UTIL/LCTrackerConf.h>

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
#include <Acts/Utilities/BinningType.hpp>

#include "Helpers.h"

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

    registerProcessorParameter("DetectorSchema", "Detector schema name (MuColl_v1, MAIA_v0, MuSIC_v1, MuSIC_v2).",
                               _detSchema, _detSchema);

    registerInputCollections(LCIO::TRACKERHITPLANE, "TrackerHitCollectionNames",
                             "Name of the TrackerHit input collections.",_inputTrackerHitCollections, {});

    registerOutputCollection(LCIO::TRACK, "TrackCollectionName", "Name of track output collection.",
                             _outputTrackCollection, string("Tracks"));

    _magCache = _magneticField->makeCache(_magneticFieldContext);
}

const Acts::Surface* ACTSBaseTracker::findSurface(const EVENT::TrackerHit* hit) const
{
    uint64_t moduleGeoId = _geoIDMappingTool->getGeometryID(hit);
    return _trackingGeometry->findSurface(moduleGeoId);
}

void ACTSBaseTracker::init()
{
    _matFile = MarlinACTS::findFile(_matFile);
    _tgeoFile = MarlinACTS::findFile(_tgeoFile);
    _tgeodescFile = MarlinACTS::findFile(_tgeodescFile);

    _runNumber = 0;
    _eventNumber = 0;

    printParameters();

    streamlog_out(MESSAGE) << " -------------------------------------" << std::endl;

    streamlog_out(MESSAGE) << " -- Building magnetic field" << std::endl;
    buildBfield();
    streamlog_out(MESSAGE) << " -- Building tracking detector" << std::endl;
    buildDetector();

    streamlog_out(MESSAGE) << " -------------------------------------" << std::endl;

    DetSchema dSchema;
    if (_detSchema == "MuSIC_v1") dSchema = DetSchema::MuSIC_v1;
    else if (_detSchema == "MuSIC_v2") dSchema = DetSchema::MuSIC_v2;
    else if (_detSchema == "MAIA_v0") dSchema = DetSchema::MAIA_v0;
    else if (_detSchema == "MuColl_v1") dSchema = DetSchema::MuColl_v1;
    else dSchema = DetSchema::MuColl_v1;

    _geoIDMappingTool = std::make_shared<GeometryIdMappingTool>(
        lcio::LCTrackerCellID::encoding_string(), dSchema);
}

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

