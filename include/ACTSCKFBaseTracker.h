#ifndef ACTSCKFBASETRACKER_H
#define ACTSCKFBASETRACKER_H 1

#include "ACTSBaseTracker.h"
#include "SourceLink.h"
#include "Measurement.h"
#include "CollectionHandler.h"

#include "Acts/Propagator/MaterialInteractor.hpp"
#include <Acts/TrackFinding/CombinatorialKalmanFilter.hpp>
#include <Acts/TrackFitting/GainMatrixUpdater.hpp>
#include <Acts/EventData/VectorTrackContainer.hpp>
#include <Acts/EventData/VectorMultiTrajectory.hpp>

#include "tbb/tbb.h"

using std::vector;

class ACTSTrackFinderThread;

/* ********************************************************************************************
 * Track finder abstract class
 ******************************************************************************************* */

class ACTSCKFBaseTracker : public ACTSBaseTracker
{
    friend ACTSTrackFinderThread;
public:
    ACTSCKFBaseTracker(const ACTSCKFBaseTracker&) = delete;
    ACTSCKFBaseTracker& operator=(const ACTSCKFBaseTracker&) = delete;
    ACTSCKFBaseTracker(const string& procname);

    virtual void init() override;

    virtual void processEvent(LCEvent* evt) override;

protected:

    using TrackContainer = Acts::TrackContainer<Acts::VectorTrackContainer,
                                                Acts::VectorMultiTrajectory,
											    std::shared_ptr>;
    using TrackFinderOptions = Acts::CombinatorialKalmanFilterOptions<TrackContainer>;

    using CKF = Acts::CombinatorialKalmanFilter<Propagator, TrackContainer>;
    using CKFPtr = std::shared_ptr<CKF>;
    using SeedParamList = vector<Acts::BoundTrackParameters>;

    using ExtrapolationOpts = Propagator::template Options<Acts::ActorList<Acts::MaterialInteractor,
                                                                           Acts::EndOfWorldReached>>;

    virtual SeedParamList getSeeds(const MarlinACTS::MeasurementContainer& m_list,
                                   LCEvent* evt) = 0;

    double _CKF_chi2CutOff = 15;
    int32_t _CKF_numMeasurementsCutOff = 10;

    double _initialTrackError_d0;
    double _initialTrackError_phi;
    double _initialTrackError_relP;
    double _initialTrackError_lambda;
    double _initialTrackError_z0;
    double _initialTrackError_time;

    float theta_tolerance;
    bool mtckf_mode;

    CKFPtr trackFinder;

private:
    void runCKFonSeed(Acts::BoundTrackParameters& seed,
                      TrackFinderOptions& ckfOptions,
                      ExtrapolationOpts& extrapolationOptions,
                      MarlinACTS::CollectionHandler& coll_handler);

    tbb::queuing_mutex p_mutex;
};

/* ********************************************************************************************
 * Track finder friend class for MT CKF
 ******************************************************************************************* */

class ACTSTrackFinderThread
{
public:
    ACTSTrackFinderThread(ACTSCKFBaseTracker& proc,
                      vector<Acts::BoundTrackParameters>& pseeds,
                      MarlinACTS::CollectionHandler& tColl,
                      ACTSCKFBaseTracker::ExtrapolationOpts& extOpts,
                      ACTSCKFBaseTracker::TrackFinderOptions& tFinderOpts) :
        processor(proc),
        paramseeds(pseeds),
        trackCollection(tColl),
        extrapolOptions(extOpts),
        ckfFinderOptions(tFinderOpts)  {}
    virtual ~ACTSTrackFinderThread() {}

    void operator()(const tbb::blocked_range<std::size_t>& prange) const;
private:
    ACTSCKFBaseTracker& processor;
    vector<Acts::BoundTrackParameters>& paramseeds;
    MarlinACTS::CollectionHandler& trackCollection;
    ACTSCKFBaseTracker::ExtrapolationOpts& extrapolOptions;
    ACTSCKFBaseTracker::TrackFinderOptions& ckfFinderOptions;
};


#endif //ACTSCKFBASETRACKER_H
