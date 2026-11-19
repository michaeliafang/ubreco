#ifndef UBRECO_EXACTHITMAPPING_H
#define UBRECO_EXACTHITMAPPING_H
#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <map>
#include <vector>

namespace showerreco { namespace recovery {
// Full copied-hit fingerprint; no geometry tolerance or event-wide pair scan.
using HitFingerprint = std::array<double,21>;
inline bool finiteFingerprint(const HitFingerprint& key) {
  for (double x : key) if (!std::isfinite(x)) return false;
  return true;
}
constexpr std::size_t noHit = std::numeric_limits<std::size_t>::max();
struct HitMapping {
  std::vector<std::size_t> sourceToSelected;
  std::size_t matched = 0, ambiguous = 0, missing = 0, invalid = 0;
};
inline HitMapping mapCopiedHits(const std::vector<HitFingerprint>& source,
                               const std::vector<HitFingerprint>& selected) {
  using Bucket = std::vector<std::size_t>;
  std::map<HitFingerprint,Bucket> src, dst;
  for (std::size_t i=0; i<source.size(); ++i)
    if (finiteFingerprint(source[i])) src[source[i]].push_back(i);
  for (std::size_t i=0; i<selected.size(); ++i)
    if (finiteFingerprint(selected[i])) dst[selected[i]].push_back(i);
  HitMapping result;
  result.sourceToSelected.assign(source.size(),noHit);
  for (const auto& key : selected) if (!finiteFingerprint(key)) ++result.invalid;
  for (const auto& entry : dst) {
    const auto found = src.find(entry.first);
    if (found == src.end()) { result.missing += entry.second.size(); continue; }
    if (found->second.size()!=1 || entry.second.size()!=1) {
      result.ambiguous += entry.second.size(); continue;
    }
    result.sourceToSelected[found->second.front()] = entry.second.front();
    ++result.matched;
  }
  return result;
}
} }
#endif
