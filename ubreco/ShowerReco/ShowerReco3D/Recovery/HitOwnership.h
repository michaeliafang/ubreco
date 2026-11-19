#ifndef UBRECO_RECOVERY_HITOWNERSHIP_H
#define UBRECO_RECOVERY_HITOWNERSHIP_H
#include <set>
#include <vector>

namespace showerreco {
namespace recovery {
// Identity must include product ID as well as hit key (art::Ptr provides this).
template <class Hit> std::vector<Hit> uniqueHitsInOrder(const std::vector<Hit>& hits) {
  std::set<Hit> seen;
  std::vector<Hit> result;
  for (const auto& hit : hits) if (seen.insert(hit).second) result.push_back(hit);
  return result;
}

template <class Hit> bool disjointHits(const std::vector<Hit>& a, const std::vector<Hit>& b) {
  const std::set<Hit> first(a.begin(), a.end());
  for (const auto& hit : b) if (first.count(hit)) return false;
  return true;
}
} // namespace recovery
} // namespace showerreco
#endif
