#include "CollectionHandler.h"

#include <IMPL/LCFlagImpl.h>
#include <EVENT/LCIO.h>

#include <algorithm>
#include <cmath>

using namespace MarlinACTS;
using EVENT::TrackState;
using vector_idx = std::vector<EVENT::Track*>::size_type;

CollectionHandler::CollectionHandler(std::string coll_name, bool dedup, float tolerance) :
    collection_name(coll_name),
    dedup_on(dedup),
    theta_tolerance(tolerance),
    t_queue()
{
    trackCollection = new IMPL::LCCollectionVec(EVENT::LCIO::TRACK);
    IMPL::LCFlagImpl trkFlag(0);
    trkFlag.setBit(EVENT::LCIO::TRBIT_HITS);
    trackCollection->setFlag(trkFlag.getFlag());
}

CollectionHandler::~CollectionHandler()
{
    if (trackCollection != nullptr)
    {
        delete trackCollection;
    }
}

void CollectionHandler::addTrack(TrackPtr track)
{
    if (dedup_on)
    {
        t_queue.push(track);
    }
    else
    {
        trackCollection->addElement(track);
    }
}

void CollectionHandler::process()
{
    if (!dedup_on) return;

    std::vector<TrackPtr> track_list;
    vector_idx t_bound = 0;

    while (!t_queue.empty())
    {
        TrackPtr curr_track = t_queue.top();

        bool foundAnEqual = false;

        float l_angle = std::atan(curr_track->getTrackState(TrackState::AtIP)->getTanLambda());
        l_angle -= theta_tolerance;
        float low_theta = std::tan(std::max(l_angle, -0.5f * static_cast<float>(M_PI)));

        while (track_list.size() > t_bound and
            low_theta > track_list[t_bound]->getTrackState(TrackState::AtIP)->getTanLambda())
        {
            t_bound++;
        }

        for (vector_idx idx = t_bound; idx < track_list.size(); ++idx)
        {
            // Skip tracks that are not equal
            if (!tracks_equal(curr_track, track_list[idx])) continue;

            foundAnEqual = true;

            // Replace if my track is better
            if (track_quality_compare(curr_track, track_list[idx]))
            {
                track_list[idx] = curr_track;
                break;
            }
        }

        if (!foundAnEqual) track_list.push_back(curr_track);

        t_queue.pop();
    }

    for (TrackPtr t_ptr : track_list) trackCollection->addElement(t_ptr);
}

void CollectionHandler::saveCollection(EVENT::LCEvent* evt)
{
    evt->addCollection(trackCollection, collection_name);
    trackCollection = nullptr;
}

bool CollectionHandler::tracks_equal(const TrackPtr trk1, const TrackPtr trk2)
{
    const EVENT::TrackerHitVec& hits1 = trk1->getTrackerHits();
    const EVENT::TrackerHitVec& hits2 = trk2->getTrackerHits();

    // Number of overlapping hist
    uint32_t hitOlap = 0;
    for (const EVENT::TrackerHit* hit1 : hits1)
    {
        if (std::find(hits2.begin(), hits2.end(), hit1) != hits2.end()) hitOlap++;
    }

    // Smaller track count
    uint32_t size = std::min(hits1.size(), hits2.size());

    return 2 * hitOlap > size;  // half of smaller track belong to larger track
}

bool CollectionHandler::track_quality_compare(const TrackPtr trk1, const TrackPtr trk2)
{
    // If number of hits are different, then the one with more
    // hits should be chosen first
    if (trk1->getTrackerHits().size() != trk2->getTrackerHits().size())
    return trk1->getTrackerHits().size() > trk2->getTrackerHits().size();

    // Same number of hits means I want smaller chi2
    return trk1->getChi2() < trk2->getChi2();
}

