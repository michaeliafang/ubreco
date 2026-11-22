#ifndef UBRECO_RECOVERY_HITOWNERSHIP_H
#define UBRECO_RECOVERY_HITOWNERSHIP_H
#include <set>
#include <array>
#include <cstddef>
#include <utility>
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
// Attach a detached fragment to the cmerger cluster with the most hits owned
// by a proposed parent, independently in each plane. Other clusters survive.
// Build transactionally: missing parent evidence leaves the input untouched.
template <class Hit, class PlaneOf>
bool attachFragment(std::vector<std::vector<Hit>>& clusters,
                    const std::vector<Hit>& parent, const std::vector<Hit>& fragment,
                    PlaneOf planeOf) {
  const auto parts=hitsByPlane(fragment,planeOf);
  const std::set<Hit> parentSet(parent.begin(),parent.end());
  const std::set<Hit> fragmentSet(fragment.begin(),fragment.end());
  std::array<std::size_t,3> target{}, support{};
  for (std::size_t c=0; c<clusters.size(); ++c) {
    const auto hits=uniqueHitsInOrder(clusters[c]);
    if (hits.empty()) continue;
    const auto pl=planeOf(hits.front());
    if (pl>=3) return false;
    std::size_t count=0;
    for (const auto& h : hits) {
      if (planeOf(h)!=pl) return false;
      count+=parentSet.count(h);
    }
    if (count>support[pl]) { support[pl]=count; target[pl]=c; }
  }
  for (unsigned pl=0; pl<3; ++pl) if (!parts[pl].empty() && !support[pl]) return false;
  for (auto& cluster : clusters) {
    std::vector<Hit> retained;
    for (const auto& hit : uniqueHitsInOrder(cluster))
      if (!fragmentSet.count(hit)) retained.push_back(hit);
    cluster=std::move(retained);
  }
  for (unsigned pl=0; pl<3; ++pl) if (!parts[pl].empty()) {
    auto& cluster=clusters[target[pl]];
    cluster.insert(cluster.end(),parts[pl].begin(),parts[pl].end());
  }
  return true;
}

} // namespace recovery
} // namespace showerreco
#endif
