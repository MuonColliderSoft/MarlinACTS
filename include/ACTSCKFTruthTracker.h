#ifndef ACTSCKFTRUTHTRACKER_H
#define ACTSCKFTRUTHTRACKER_H 1

#include "ACTSCKFBaseTracker.h"

class ACTSCKFTruthTracker : public ACTSCKFBaseTracker
{
public:
    ACTSCKFTruthTracker(const ACTSCKFTruthTracker&) = delete;
    ACTSCKFTruthTracker& operator=(const ACTSCKFTruthTracker&) = delete;
    ACTSCKFTruthTracker();

    virtual marlin::Processor *newProcessor() { return new ACTSCKFTruthTracker; }

    virtual void init() override;

protected:
    virtual SeedParamList getSeeds(const MarlinACTS::MeasurementContainer& m_list,
                                   LCEvent* evt) override;

private:
    std::string _inputParticleCollection;
};

#endif //ACTSCKFTRUTHTRACKER_H
