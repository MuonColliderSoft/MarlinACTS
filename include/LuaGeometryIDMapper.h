#ifndef LUAIDMAP_H
#define LUAIDMAP_H 1

#include <Acts/Geometry/GeometryIdentifier.hpp>

#include <EVENT/SimTrackerHit.h>
#include <EVENT/TrackerHit.h>
#include <UTIL/CellIDDecoder.h>

#include "lua.hpp"

#include <string>
#include <unordered_map>
#include <optional>

using std::string;
using std::optional;

using ObjectMap = std::unordered_map<uint32_t, uint64_t>;
using OptMapID = optional<uint64_t>;

namespace MarlinACTS {

class LuaGeometryIDMapper
{
public:
    LuaGeometryIDMapper(const string& script, const string& encoderString);
    virtual ~LuaGeometryIDMapper();

    Acts::GeometryIdentifier getGeometryID(const lcio::SimTrackerHit* hit);
    Acts::GeometryIdentifier getGeometryID(const lcio::TrackerHit* hit);

    Acts::GeometryIdentifier getGeometryID(uint32_t systemID, uint32_t layerID,
            int32_t sideID, uint32_t ladderID, uint32_t moduleID);
private:
    static const int32_t VertexEndCapNegative;
    static const int32_t VertexBarrel;
    static const int32_t VertexEndCapPositive;
    static const int32_t InnerTrackerEndCapNegative;
    static const int32_t InnerTrackerBarrel;
    static const int32_t InnerTrackerEndCapPositive;
    static const int32_t OuterInnerTrackerEndCapNegative;
    static const int32_t OuterInnerTrackerBarrel;
    static const int32_t OuterInnerTrackerEndCapPositive;
    static const int32_t OuterTrackerEndCapNegative;
    static const int32_t OuterTrackerBarrel;
    static const int32_t OuterTrackerEndCapPositive;

    uint32_t hash(int32_t sideID, uint32_t systemID, uint32_t layerID);
    OptMapID call_map(string fnct, int32_t systemID, uint32_t layerID);
    OptMapID call_map(string fnct, int32_t systemID, uint32_t layerID,
                      uint32_t ladderID, uint32_t moduleID);
    uint64_t map_volume(uint32_t systemID, uint32_t layerID);
    uint64_t map_layer(uint32_t systemID, uint32_t layerID);
    uint64_t map_sensor(uint32_t systemID, uint32_t layerID, uint32_t ladderID, uint32_t moduleID);

    string _encoderString;
    lua_State* l_ctx;

    // Caches for volumes, layers and sensitives
    ObjectMap volume_map;
    ObjectMap layer_map;
    ObjectMap nlad_map;
};

}  // namespace MarlinACTS

#endif // LUAIDMAP_H
