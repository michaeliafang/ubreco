#ifndef UBRECO_SPACEPOINTCONSISTENCY_H
#define UBRECO_SPACEPOINTCONSISTENCY_H
#include <array>
#include <cstddef>
#include <map>
#include <vector>

namespace showerreco { namespace recovery {
struct AssignedPointHit { std::size_t hit; unsigned int plane, shower; };
struct PointConsistency {
  std::size_t informative = 0;
  double weight = 0., conflictingWeight = 0.;
  double fraction() const { return weight > 0. ? conflictingWeight/weight : 0.; }
};
// A point is informative only when its assigned hits span at least two planes.
// Shared-hit ambiguities downweight points: one hit appearing in many candidate
// points must not create many independent votes against a split.
inline PointConsistency pointConsistency(const std::vector<std::vector<AssignedPointHit>>& points) {
  std::map<std::size_t,std::size_t> multiplicity;
  for (const auto& p : points) for (const auto& h : p) ++multiplicity[h.hit];
  PointConsistency result;
  for (const auto& p : points) {
    std::array<bool,3> planes{{false,false,false}};
    double weight = 1.;
    unsigned int first = 0;
    bool have = false, conflict = false;
    for (const auto& h : p) {
      if (h.plane >= 3) continue;
      planes[h.plane] = true;
      if (!have) { first=h.shower; have=true; }
      conflict = conflict || h.shower!=first;
      const double w = 1./multiplicity.at(h.hit);
      if (w<weight) weight=w;
    }
    if (unsigned(planes[0])+unsigned(planes[1])+unsigned(planes[2]) < 2) continue;
    ++result.informative;
    result.weight += weight;
    if (conflict) result.conflictingWeight += weight;
  }
  return result;
}
} }
#endif
