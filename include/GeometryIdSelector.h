#ifndef GEOMETRYIDSELECTOR_H
#define GEOMETRYIDSELECTOR_H 1

#include <Acts/Geometry/GeometryIdentifier.hpp>
#include <vector>

namespace MarlinACTS {

class GeometryIdSelector
{
public:
    using Mask = Acts::GeometryIdentifier::Value;

    GeometryIdSelector() = default;
    GeometryIdSelector(const std::vector<Acts::GeometryIdentifier>& selection);

    bool check(const Acts::GeometryIdentifier& geoID);

    static Mask makeMask(const Acts::GeometryIdentifier& id);

private:
    std::vector<std::pair<Acts::GeometryIdentifier, Mask>> m_selection;
};

}  // namespace ACTSTracking

#endif // GEOMETRYIDSELECTOR_H
