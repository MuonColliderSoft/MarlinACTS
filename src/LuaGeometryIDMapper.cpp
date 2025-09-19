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
        throw "Script loading failure";
    }
    if (lua_pcall(l_ctx, 0, 0, 0) == LUA_OK)
    {
        lua_pop(l_ctx, lua_gettop(l_ctx));
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
    int32_t signSystemID = (sideID < 0) ? -1 * systemID : systemID;

    uint32_t key = hash(sideID, systemID, layerID);

    // map volume id
    if (volume_map.find(key) == volume_map.end())
    {
        auto vol_id = call_map("getVolumeID", signSystemID, layerID);
        if (vol_id) volume_map.emplace(key, vol_id.value());
    }
    geometry_id |= volume_map[key] << (14 * 4);

    // map layer id
    if (layer_map.find(key) == layer_map.end())
    {
        auto lay_id = call_map("getLayerID", signSystemID, layerID);
        if (lay_id) layer_map.emplace(key, lay_id.value());
    }
    geometry_id |= layer_map[key] << (9 * 4);

    // map sensitive id
    uint64_t sens_id = 0;
    if (nlad_map.find(key) != nlad_map.end())
    {
        sens_id = nlad_map[key] * ladderID + moduleID + 1;
    }
    else
    {
        auto lad_num = call_map("getLadderNumber", signSystemID, layerID);
        if (lad_num)
        {
            nlad_map.emplace(key, lad_num.value());
            sens_id = nlad_map[key] * ladderID + moduleID + 1;
        }
        else
        {
            auto res = call_map("getSensorID", signSystemID, layerID, ladderID, moduleID);
            if (res) sens_id = res.value();
        }
    }
    geometry_id |= sens_id << (2 * 4);

    return Acts::GeometryIdentifier { geometry_id };
}

uint32_t LuaGeometryIDMapper::hash(int32_t sideID, uint32_t systemID, uint32_t layerID)
{
    uint32_t res = systemID * 10000 + layerID;
    return sideID < 0 ? 1000000 + res : res;
}

OptMapID LuaGeometryIDMapper::call_map(string fnct, int32_t systemID, uint32_t layerID)
{
    lua_getglobal(l_ctx, fnct.c_str());
    if (not lua_isfunction(l_ctx, -1))
    {
        throw "Not a function";
    }
    lua_pushinteger(l_ctx, systemID);
    lua_pushinteger(l_ctx, layerID);
    if (lua_pcall(l_ctx, 2, 1, 0) == LUA_OK)
    {
        if (lua_isinteger(l_ctx, -1))
        {
            uint64_t res = lua_tointeger(l_ctx, -1);
            lua_pop(l_ctx, 1);
            return res;
        }
        else if (lua_isnil(l_ctx, -1))
        {
            return std::nullopt;
        }
    }
    lua_pop(l_ctx, 1);
    throw "Call failure";
}

OptMapID LuaGeometryIDMapper::call_map(string fnct, int32_t systemID, uint32_t layerID,
                                       uint32_t ladderID, uint32_t moduleID)
{
    lua_getglobal(l_ctx, fnct.c_str());
    if (not lua_isfunction(l_ctx, -1))
    {
        throw "Not a function";
    }
    lua_pushinteger(l_ctx, systemID);
    lua_pushinteger(l_ctx, layerID);
    lua_pushinteger(l_ctx, ladderID);
    lua_pushinteger(l_ctx, moduleID);
    if (lua_pcall(l_ctx, 4, 1, 0) == LUA_OK)
    {
        if (lua_isinteger(l_ctx, -1))
        {
            uint64_t res = lua_tointeger(l_ctx, -1);
            lua_pop(l_ctx, 1);
            return res;
        }
        else if (lua_isnil(l_ctx, -1))
        {
            return std::nullopt;
        }
    }
    lua_pop(l_ctx, 1);
    throw "Call failure";
}



