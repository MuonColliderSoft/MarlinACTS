#ifndef ACTSBASETRACKER_H
#define ACTSBASETRACKER_H 1

#include <Acts/Geometry/GeometryContext.hpp>
#include <Acts/Geometry/TrackingGeometry.hpp>
#include <Acts/MagneticField/MagneticFieldContext.hpp>
#include <Acts/MagneticField/MagneticFieldProvider.hpp>
#include <Acts/Plugins/TGeo/TGeoDetectorElement.hpp>
#include <Acts/Utilities/CalibrationContext.hpp>
#include <Acts/Propagator/EigenStepper.hpp>
#include <Acts/Propagator/Navigator.hpp>
#include <Acts/Propagator/Propagator.hpp>
#include <Acts/Surfaces/PerigeeSurface.hpp>
#include <Acts/EventData/VectorTrackContainer.hpp>
#include <Acts/EventData/VectorMultiTrajectory.hpp>

#include <marlin/Processor.h>

#include "GeometryIdMappingTool.h"

#include <IMPL/LCCollectionVec.h>
#include <EVENT/Track.h>
#include <EVENT/TrackState.h>

#include <memory>

using MarlinACTS::GeometryIdMappingTool;
using std::string;
using std::vector;

class ACTSBaseTracker : public marlin::Processor
{
public:
    using DetectorElementPtr = std::shared_ptr<const Acts::TGeoDetectorElement>;
    using DetectorStore = vector<DetectorElementPtr>;

    ACTSBaseTracker(const ACTSBaseTracker&) = delete;
    ACTSBaseTracker& operator=(const ACTSBaseTracker&) = delete;
    ACTSBaseTracker(const string& procname);

    virtual void init();

    virtual void processRunHeader(LCRunHeader* run) { _runNumber++; }

    virtual void processEvent(LCEvent* evt) { _eventNumber++; }

    virtual void check(LCEvent* evt) {}

    virtual void end();

protected:
    using Stepper = Acts::EigenStepper<>;
    using Navigator = Acts::Navigator;
    using Propagator = Acts::Propagator<Stepper, Navigator>;
    using PropagatorPtr = std::shared_ptr<Propagator>;
    using TrackResult = Acts::TrackContainer<Acts::VectorTrackContainer,
                                             Acts::VectorMultiTrajectory,
                                             std::shared_ptr>::TrackProxy;


    string _matFile {};
    string _tgeoFile = "data/MuSIC_v2.root";
    string _tgeodescFile = "data/MuSIC_v2.json";
    string _detSchema = "MuSIC_v2";

    vector<string> _inputTrackerHitCollections;
    string _outputTrackCollection;

    uint32_t _eventNumber;
    uint32_t _runNumber;
    uint32_t _fitFails;

    float _caloFaceR = 1857; //mm
    float _caloFaceZ = 2307; //mm

    PropagatorPtr propagator;

    std::shared_ptr<LCCollectionVec> trackCollection;

    std::shared_ptr<Acts::PerigeeSurface> perigeeSurface;

    std::shared_ptr<GeometryIdMappingTool> geoIDMappingTool() const { return _geoIDMappingTool; }

    const Acts::MagneticFieldContext& magneticFieldContext() const { return _magneticFieldContext; }

    const Acts::GeometryContext& geometryContext() const { return _geometryContext; }

    const Acts::CalibrationContext& calibrationContext() const { return _calibrationContext; }

    std::shared_ptr<Acts::MagneticFieldProvider> magneticField() const { return _magneticField; }

    std::shared_ptr<const Acts::TrackingGeometry> trackingGeometry() const { return _trackingGeometry; }

    const Acts::Surface* findSurface(const EVENT::TrackerHit* hit) const;

    LCCollection* getCollection(const string& collectionName, LCEvent* evt);

    Acts::Vector3 magneticFieldValue(const Acts::Vector3 position);

    virtual void storeData(LCEvent* evt) { evt->addCollection(trackCollection.get(), _outputTrackCollection); }

    virtual EVENT::Track* convert_track(const TrackResult& fitter_res);
    virtual EVENT::TrackState* convert_state(int location, const Acts::BoundVector& value,
                                             const Acts::BoundMatrix& cov, Acts::Vector3 mag_pos);


private:
    void buildDetector();

    void buildBfield();

    std::shared_ptr<GeometryIdMappingTool> _geoIDMappingTool;

    Acts::MagneticFieldContext _magneticFieldContext;
    std::shared_ptr<Acts::MagneticFieldProvider> _magneticField;
    Acts::MagneticFieldProvider::Cache _magCache;

    Acts::GeometryContext _geometryContext;
    DetectorStore _detectorStore;
    std::shared_ptr<const Acts::TrackingGeometry> _trackingGeometry = nullptr;

    Acts::CalibrationContext _calibrationContext;
};

#endif //ACTSBASETRACKER_H
