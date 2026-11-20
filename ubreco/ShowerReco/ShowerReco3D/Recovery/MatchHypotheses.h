#ifndef UBRECO_RECOVERY_MATCHHYPOTHESES_H
#define UBRECO_RECOVERY_MATCHHYPOTHESES_H
#include <array>
#include <algorithm>
#include <cmath>
#include <vector>
#include <set>
namespace showerreco { namespace recovery {
using ClusterMatches = std::vector<std::vector<unsigned int>>;
struct MatchHypothesis {
  ClusterMatches matches;
  unsigned planeMask = 7, swapMask = 0;
  bool twoPlaneFallback = false;
  bool clusterSearch = false;
  std::size_t selectedHits = 0;
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

// Independent of the matcher's grouping: retain the largest few clusters per
// plane, choose two per view, and test both relative assignments. Charge
// coverage is measured against ALL eligible clusters, including those beyond
// the search cap, so truncation cannot hide substantial extra activity.
inline std::vector<MatchHypothesis> substantialPairHypotheses(
    const std::vector<unsigned>& planes, const std::vector<std::size_t>& hits,
    const std::vector<double>& charge, std::size_t minHits,
    std::size_t maxPerPlane, double minChargeCoverage) {
  if (planes.size()!=hits.size() || planes.size()!=charge.size() || maxPerPlane<2) return {};
  std::array<std::vector<unsigned>,3> pool;
  std::array<double,3> total{};
  for (unsigned c=0; c<planes.size(); ++c) {
    if (planes[c]>=3 || !std::isfinite(charge[c]) || charge[c]<0.) return {};
    if (hits[c]<minHits) continue;
    total[planes[c]]+=charge[c];
    pool[planes[c]].push_back(c);
  }
  for (auto& view : pool) {
    std::stable_sort(view.begin(),view.end(),[&](unsigned a,unsigned b) {
      if (hits[a]!=hits[b]) return hits[a]>hits[b];
      return charge[a]>charge[b];
    });
    if (view.size()>maxPerPlane) view.resize(maxPerPlane);
  }
  std::vector<MatchHypothesis> result;
  for (unsigned induction : {0u,1u}) {
    const auto& a=pool[induction]; const auto& b=pool[2];
    for (std::size_t i=0; i<a.size(); ++i) for (std::size_t j=i+1; j<a.size(); ++j) {
      if (!(total[induction]>0.) || charge[a[i]]+charge[a[j]]<minChargeCoverage*total[induction]) continue;
      for (std::size_t k=0; k<b.size(); ++k) for (std::size_t l=k+1; l<b.size(); ++l) {
        if (!(total[2]>0.) || charge[b[k]]+charge[b[l]]<minChargeCoverage*total[2]) continue;
        for (unsigned swap=0; swap<2; ++swap) {
          MatchHypothesis h;
          h.matches={{a[i],b[swap ? l : k]},{a[j],b[swap ? k : l]}};
          h.planeMask=(1u<<induction)|4u; h.swapMask=swap ? 4u : 0u;
          h.twoPlaneFallback=true; h.clusterSearch=true;
          h.selectedHits=hits[a[i]]+hits[a[j]]+hits[b[k]]+hits[b[l]];
          result.push_back(std::move(h));
        }
      }
    }
  }
  std::stable_sort(result.begin(),result.end(),[](const auto& a,const auto& b) {
    return a.selectedHits>b.selectedHits;
  });
  return result;
}

} }
#endif
