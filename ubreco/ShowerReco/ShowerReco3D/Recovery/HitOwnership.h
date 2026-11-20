#ifndef UBRECO_RECOVERY_HITOWNERSHIP_H
#define UBRECO_RECOVERY_HITOWNERSHIP_H
#include <set>
#include <array>
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
// Combine supplied shower hits without recovering unrelated singleton clusters.
// Deduplicate by complete hit identity before grouping by detector plane.
template <class Hit, class PlaneOf>
std::array<std::vector<Hit>,3> hitsByPlane(const std::vector<Hit>& hits, PlaneOf planeOf) {
  std::array<std::vector<Hit>,3> result;
  for (const auto& hit : uniqueHitsInOrder(hits)) result.at(planeOf(hit)).push_back(hit);
  return result;
}
} // namespace recovery
} // namespace showerreco
#endif
