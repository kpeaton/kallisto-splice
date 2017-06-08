#ifndef EXONMAP_H
#define EXONMAP_H

#include <string>
#include <utility>
#include <vector>
#include <map>

typedef std::map<std::string, std::tuple<std::string, std::string, std::vector<std::vector<int>>>> ExonMap;

#endif // EXONMAP_H