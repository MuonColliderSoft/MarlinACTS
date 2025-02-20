#ifndef ACTSCKFBASETRACKER_H
#define ACTSCKFBASETRACKER_H 1

#include "ACTSBaseTracker.h"
#include "SourceLink.h"

#include "Acts/Propagator/MaterialInteractor.hpp"
#include <Acts/TrackFinding/CombinatorialKalmanFilter.hpp>
#include <Acts/TrackFitting/GainMatrixUpdater.hpp>
#include <Acts/EventData/VectorTrackContainer.hpp>
#include <Acts/EventData/VectorMultiTrajectory.hpp>

class ACTSCKFBaseTracker : public ACTSBaseTracker
{
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
    using TrackFinderOptions = Acts::CombinatorialKalmanFilterOptions<
        MarlinACTS::SourceLinkAccessor::Iterator, TrackContainer>;

    using CKF = Acts::CombinatorialKalmanFilter<Propagator, TrackContainer>;
    using CKFPtr = std::shared_ptr<CKF>;

    double _CKF_chi2CutOff = 15;
    int32_t _CKF_numMeasurementsCutOff = 10;

    double _initialTrackError_d0;
    double _initialTrackError_phi;
    double _initialTrackError_relP;
    double _initialTrackError_lambda;
    double _initialTrackError_z0;
    double _initialTrackError_time;

    CKFPtr trackFinder;
};

#endif //ACTSCKFBASETRACKER_H
