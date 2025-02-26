#ifndef HELPERS_H
#define HELPERS_H 1

#include <string>
#include <tuple>

namespace MarlinACTS {

std::string findFile(const std::string& inpath);

using ZBinSchema = std::tuple<float, int, int, int, int, bool>;

ZBinSchema parseZBinSchema(std::string schema_str);

} //namespace MarlinACTS

#endif //HELPERS_H
