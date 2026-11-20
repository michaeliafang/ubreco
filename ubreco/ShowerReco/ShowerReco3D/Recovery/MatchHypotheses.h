#ifndef UBRECO_RECOVERY_MATCHHYPOTHESES_H
#define UBRECO_RECOVERY_MATCHHYPOTHESES_H
#include <array>
#include <vector>
#include <set>
namespace showerreco { namespace recovery {
using ClusterMatches = std::vector<std::vector<unsigned int>>;
struct MatchHypothesis {
  ClusterMatches matches;
  unsigned planeMask = 7, swapMask = 0;
  bool twoPlaneFallback = false;
};
// Only permute an unambiguous pair with exactly one cluster per represented
// plane per shower. Fix the first plane's ordering to remove global relabeling.
inline std::vector<MatchHypothesis> pairingHypotheses(const ClusterMatches& matches,
                                                     const std::vector<unsigned>& planes) {
  if (matches.size()!=2) return {};
  std::array<std::array<unsigned,3>,2> byPlane{};
  std::array<unsigned,2> masks{};
  std::set<unsigned> seen;
  for (unsigned s=0; s<2; ++s) for (auto ci : matches[s]) {
    if (ci>=planes.size() || planes[ci]>=3 || !seen.insert(ci).second) return {};
    const unsigned pl=planes[ci], bit=1u<<pl;
    if (masks[s]&bit) return {};
    masks[s]|=bit; byPlane[s][pl]=ci;
  }
  if (masks[0]!=masks[1] || (masks[0]&(masks[0]-1))==0) return {};
  const unsigned fixed=masks[0]&(~masks[0]+1);
  std::vector<MatchHypothesis> result;
  for (unsigned swap=0; swap<8; ++swap) {
    if ((swap&~masks[0]) || (swap&fixed)) continue;
    MatchHypothesis h; h.matches.resize(2); h.planeMask=masks[0]; h.swapMask=swap;
    for (unsigned s=0; s<2; ++s) for (unsigned pl=0; pl<3; ++pl)
      if (masks[0]&(1u<<pl)) h.matches[s].push_back(byPlane[s^unsigned(bool(swap&(1u<<pl)))][pl]);
    result.push_back(h);
  }
  return result;
}
inline std::vector<MatchHypothesis> twoPlaneHypotheses(const ClusterMatches& matches,
                                                      const std::vector<unsigned>& planes,
                                                      bool allowSwaps) {
  if (matches.size()!=2) return {};
  std::vector<MatchHypothesis> result;
  for (unsigned mask : {5u,6u}) { // collection + one induction plane only
    ClusterMatches retained(2);
    bool omitted=false;
    for (unsigned s=0; s<2; ++s) for (auto ci : matches[s]) {
      if (ci>=planes.size() || planes[ci]>=3) return {};
      if (mask&(1u<<planes[ci])) retained[s].push_back(ci);
      else omitted=true;
    }
    if (!omitted) continue; // already a two-plane match, not a fallback
    for (auto h : pairingHypotheses(retained,planes)) {
      if (h.planeMask!=mask || (!allowSwaps && h.swapMask)) continue;
      h.twoPlaneFallback=true;
      result.push_back(h);
    }
  }
  return result;
}

} }
#endif
