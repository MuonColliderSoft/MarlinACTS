#include "LuaGeometryIDMapper.h"

using MarlinACTS::LuaGeometryIDMapper;

const int32_t LuaGeometryIDMapper::VertexEndCapNegative = -2;
const int32_t LuaGeometryIDMapper::VertexBarrel = 1;
const int32_t LuaGeometryIDMapper::VertexEndCapPositive = 2;
const int32_t LuaGeometryIDMapper::InnerTrackerEndCapNegative = -4;
const int32_t LuaGeometryIDMapper::InnerTrackerBarrel = 3;
const int32_t LuaGeometryIDMapper::InnerTrackerEndCapPositive = 4;
const int32_t LuaGeometryIDMapper::OuterInnerTrackerEndCapNegative = -8;
const int32_t LuaGeometryIDMapper::OuterInnerTrackerBarrel = 7;
const int32_t LuaGeometryIDMapper::OuterInnerTrackerEndCapPositive = 8;
const int32_t LuaGeometryIDMapper::OuterTrackerEndCapNegative = -6;
const int32_t LuaGeometryIDMapper::OuterTrackerBarrel = 5;
const int32_t LuaGeometryIDMapper::OuterTrackerEndCapPositive = 6;

LuaGeometryIDMapper::LuaGeometryIDMapper(const string& script, const string& encoderString) :
    _encoderString(encoderString),
    l_ctx(luaL_newstate())
{
    luaL_openlibs(l_ctx);
    if (luaL_loadstring(l_ctx, script.c_str()) != LUA_OK)
    {
        // TODO missing error handling
    }
}

LuaGeometryIDMapper::~LuaGeometryIDMapper()
{
    if(l_ctx) lua_close(l_ctx);
}

Acts::GeometryIdentifier LuaGeometryIDMapper::getGeometryID(const lcio::SimTrackerHit* hit)
{
    UTIL::CellIDDecoder<lcio::SimTrackerHit> decoder(_encoderString);
    uint32_t systemID = decoder(hit)["system"];
    uint32_t layerID = decoder(hit)["layer"];
    int32_t sideID = decoder(hit)["side"];
    uint32_t ladderID = decoder(hit)["module"];
    uint32_t moduleID = decoder(hit)["sensor"];

    return getGeometryID(systemID, layerID, sideID, ladderID, moduleID);
}

Acts::GeometryIdentifier LuaGeometryIDMapper::getGeometryID(const lcio::TrackerHit* hit)
{
    UTIL::CellIDDecoder<lcio::TrackerHit> decoder(_encoderString);
    uint32_t systemID = decoder(hit)["system"];
    uint32_t layerID = decoder(hit)["layer"];
    int32_t sideID = decoder(hit)["side"];
    uint32_t ladderID = decoder(hit)["module"];
    uint32_t moduleID = decoder(hit)["sensor"];

    return getGeometryID(systemID, layerID, sideID, ladderID, moduleID);
}

Acts::GeometryIdentifier LuaGeometryIDMapper::getGeometryID(uint32_t systemID,
        uint32_t layerID, int32_t sideID, uint32_t ladderID, uint32_t moduleID)
{
    uint64_t geometry_id = 0;

    // endcap is split in +/- sides by ACTS
    int32_t signSystemID = (sideID < 0) ? -systemID : systemID;

    uint64_t volume_id = map_volume(signSystemID, layerID);
    geometry_id |= volume_id << (14 * 4);

    uint64_t layer_id = map_layer(signSystemID, layerID);
    geometry_id |= layer_id << (9 * 4);


    return Acts::GeometryIdentifier { geometry_id };
}




uint32_t LuaGeometryIDMapper::hash(uint32_t systemID, uint32_t layerID)
{
    return 0;
}

uint64_t LuaGeometryIDMapper::map_volume(uint32_t systemID, uint32_t layerID)
{
    uint32_t key = hash(systemID, layerID);
    if (volume_map.find(key) != volume_map.end()) return volume_map[key];

    lua_getglobal(l_ctx, "getVolumeID");
    if (not lua_isfunction(l_ctx, -1))
    {
        // TODO missing error handling
    }
    lua_pushinteger(l_ctx, systemID);
    lua_pushinteger(l_ctx, layerID);
    if (lua_pcall(l_ctx, 2, 1, 0) == LUA_OK and lua_isinteger(l_ctx, -1))
    {
        uint64_t res = lua_tointeger(l_ctx, -1);
        volume_map.emplace(key, res);
        lua_pop(l_ctx, 1);
        return res;
    }
    lua_pop(l_ctx, 1);
    // TODO missing error handling
    return -1;
}

uint64_t LuaGeometryIDMapper::map_layer(uint32_t systemID, uint32_t layerID)
{
    uint32_t key = hash(systemID, layerID);
    if (layer_map.find(key) != layer_map.end()) return layer_map[key];

    lua_getglobal(l_ctx, "getLayerID");
    if (not lua_isfunction(l_ctx, -1))
    {
        // TODO missing error handling
    }
    lua_pushinteger(l_ctx, systemID);
    lua_pushinteger(l_ctx, layerID);
    if (lua_pcall(l_ctx, 2, 1, 0) == LUA_OK and lua_isinteger(l_ctx, -1))
    {
        uint64_t res = lua_tointeger(l_ctx, -1);
        layer_map.emplace(key, res);
        lua_pop(l_ctx, 1);
        return res;
    }
    lua_pop(l_ctx, 1);
    // TODO missing error handling
    return -1;
}









