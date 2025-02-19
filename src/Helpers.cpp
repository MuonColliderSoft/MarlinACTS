#include "Helpers.h"

#include <filesystem>

namespace MarlinACTS {

std::string findFile(const std::string& inpath)
{
    if (inpath.empty()) return inpath;

    if (inpath[0] == '/') return inpath;

    if (std::filesystem::exists(inpath)) return inpath;

    return inpath;
}

} //namespace MarlinACTS

