#ifndef SOURCELINK_H
#define SOURCELINK_H 1

#include <EVENT/TrackerHit.h>

#include "Acts/EventData/SourceLink.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "GeometryContainers.h"

namespace MarlinACTS {

class SourceLink final
{
public:
    SourceLink(Acts::GeometryIdentifier gid, std::size_t index, EVENT::TrackerHit* lciohit) :
        m_geometryId(gid), m_index(index), m_lciohit(lciohit) {}

    SourceLink() = default;
    SourceLink(const SourceLink&) = default;
    SourceLink(SourceLink&&) = default;
    SourceLink& operator=(const SourceLink&) = default;
    SourceLink& operator=(SourceLink&&) = default;

    constexpr Acts::GeometryIdentifier geometryId() const { return m_geometryId; }

    constexpr std::size_t index() const { return m_index; }

    constexpr EVENT::TrackerHit* lciohit() const { return m_lciohit; }

private:
    Acts::GeometryIdentifier m_geometryId;
    std::size_t m_index = -1;
    EVENT::TrackerHit* m_lciohit = nullptr;

    friend constexpr bool operator==(const SourceLink& lhs, const SourceLink& rhs)
    {
        return (lhs.m_geometryId == rhs.m_geometryId) and
               (lhs.m_index == rhs.m_index) and (lhs.m_lciohit == rhs.m_lciohit);
    }
    friend constexpr bool operator!=(const SourceLink& lhs, const SourceLink& rhs)
    {
        return not(lhs == rhs);
    }
};

using SourceLinkContainer = GeometryIdMultiset<SourceLink>;

/*
    Accessor for the above source link container
    It wraps up a few lookup methods to be used in the Combinatorial Kalman Filter
*/
struct SourceLinkAccessor : GeometryIdMultisetAccessor<SourceLink>
{
    using BaseIterator = GeometryIdMultisetAccessor<SourceLink>::Iterator;
    using Iterator = Acts::SourceLinkAdapterIterator<BaseIterator>;

    std::pair<Iterator, Iterator> range(const Acts::Surface& surface) const
    {
        assert(container != nullptr);
        auto [begin, end] = container->equal_range(surface.geometryId());
        return {Iterator{begin}, Iterator{end}};
    }
};

struct SurfaceAccessor
{
    const std::shared_ptr<const Acts::TrackingGeometry> trackingGeometry;

    const Acts::Surface* operator()(const Acts::SourceLink& sourceLink) const
    {
        const auto& mySourceLink = sourceLink.get<SourceLink>();
        return trackingGeometry->findSurface(mySourceLink.geometryId());
    }
};

}  // namespace MarlinACTS

#endif // SOURCELINK_H
