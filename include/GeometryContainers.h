#ifndef GEOMCONTAINERS_H
#define GEOMCONTAINERS_H 1

#include <Acts/Geometry/GeometryIdentifier.hpp>
#include <boost/container/flat_set.hpp>

namespace MarlinACTS {

struct GeometryIdGetter
{
    constexpr Acts::GeometryIdentifier operator()(Acts::GeometryIdentifier geometryId) const
    {
        return geometryId;
    }
    constexpr Acts::GeometryIdentifier operator()(Acts::GeometryIdentifier::Value encoded) const
    {
        return Acts::GeometryIdentifier(encoded);
    }
    template <typename T>
    constexpr Acts::GeometryIdentifier operator()(const std::pair<Acts::GeometryIdentifier,
                                                  T>& mapItem) const
    {
        return mapItem.first;
    }
    template <typename T>
    inline auto operator()(const T& thing) const
      -> decltype(thing.geometryId(), Acts::GeometryIdentifier())
    {
        return thing.geometryId();
    }
};

namespace detail {
struct CompareGeometryId
{
    using is_transparent = void;

    template <typename Left, typename Right>
    constexpr bool operator()(Left&& lhs, Right&& rhs) const
    {
        return GeometryIdGetter()(lhs) < GeometryIdGetter()(rhs);
    }
};
}  // namespace detail

template <typename T>
using GeometryIdMultiset = boost::container::flat_multiset<T, detail::CompareGeometryId>;

template <typename T>
using GeometryIdMultimap = GeometryIdMultiset<std::pair<Acts::GeometryIdentifier, T>>;

template <typename T>
struct GeometryIdMultisetAccessor
{
    using Container = GeometryIdMultiset<T>;
    using Key = Acts::GeometryIdentifier;
    using Value = typename GeometryIdMultiset<T>::value_type;
    using Iterator = typename GeometryIdMultiset<T>::const_iterator;

    const Container* container = nullptr;

    size_t count(const Acts::GeometryIdentifier& geoId) const
    {
        assert(container != nullptr);
        return container->count(geoId);
    }

    std::pair<Iterator, Iterator> range(const Acts::GeometryIdentifier& geoId) const
    {
        assert(container != nullptr);
        return container->equal_range(geoId);
    }

    std::pair<Iterator, Iterator> range(const Acts::Surface& surface) const
    {
        assert(container != nullptr);
        auto [begin, end] = container->equal_range(surface.geometryId());
        return {Iterator{begin}, Iterator{end}};
    }

    // get the element using the iterator
    const Value& at(const Iterator& it) const { return *it; }
};
}  // namespace MarlinACTS

#endif //GEOMCONTAINERS_H
