#ifndef COLLECTIONHANDLER_H
#define COLLECTIONHANDLER_H 1

#include <EVENT/LCEvent.h>
#include <IMPL/TrackImpl.h>
#include <IMPL/TrackerHitPlaneImpl.h>
#include <IMPL/LCCollectionVec.h>

#include <string>
#include <vector>
#include <queue>

namespace MarlinACTS {

class CollectionHandler
{
public:

    using TrackPtr = EVENT::Track*;

    CollectionHandler(std::string coll_name, bool dedup, float tolerance);
    virtual ~CollectionHandler();
    void addTrack(TrackPtr track);
    void process();
    void saveCollection(EVENT::LCEvent* evt);

private:

    struct track_compare
    {
        bool operator()(const TrackPtr trk1,
                        const TrackPtr trk2) const
        {
            // reversed order for a priority queue
	    return trk1->getTrackState(EVENT::TrackState::AtIP)->getTanLambda() >
                   trk2->getTrackState(EVENT::TrackState::AtIP)->getTanLambda();
        }
    };

    bool tracks_equal(const TrackPtr trk1, const TrackPtr trk2);

    bool track_quality_compare(const TrackPtr trk1, const TrackPtr trk2);

    std::string collection_name;
    bool dedup_on;
    float theta_tolerance;

    std::priority_queue<TrackPtr, std::vector<TrackPtr>, track_compare> t_queue;

    IMPL::LCCollectionVec* trackCollection;
};


} // namespace MarlinACTS
#endif //COLLECTIONHANDLER_H
