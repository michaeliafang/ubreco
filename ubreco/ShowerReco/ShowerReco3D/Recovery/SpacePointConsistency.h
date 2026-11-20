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
struct PointDiagnostics {
  std::array<PointConsistency,3> pairs; // 01, 02, 12; independently normalized
  PointConsistency unique;
  std::array<double,4> flipConflict{}; // fix plane 0, flip plane 1/2 labels
};
inline PointDiagnostics diagnosePoints(const std::vector<std::vector<AssignedPointHit>>& points,
                                       bool twoShowers) {
  PointDiagnostics result;
  const unsigned pairs[3][2] = {{0,1},{0,2},{1,2}};
  for (unsigned k=0; k<3; ++k) {
    std::vector<std::vector<AssignedPointHit>> selected;
    for (const auto& point : points) {
      std::vector<AssignedPointHit> hits;
      for (const auto& h : point)
        if (h.plane==pairs[k][0] || h.plane==pairs[k][1]) hits.push_back(h);
      selected.push_back(std::move(hits));
    }
    result.pairs[k] = pointConsistency(selected);
  }
  // A conservative uniqueness proxy, independent of proposed shower labels.
  // This is not a guarantee of physical correctness or solver fit quality.
  std::map<std::size_t,std::size_t> uses;
  for (const auto& point : points) for (const auto& h : point) ++uses[h.hit];
  std::vector<std::vector<AssignedPointHit>> unique;
  for (const auto& point : points) {
    std::array<unsigned,3> counts{};
    bool good = true;
    for (const auto& h : point) {
      if (h.plane>=3 || uses[h.hit]!=1) { good=false; break; }
      if (++counts[h.plane]>1) { good=false; break; }
    }
    if (good) unique.push_back(point);
  }
  result.unique = pointConsistency(unique);
  if (twoShowers) for (unsigned mask=0; mask<4; ++mask) {
    auto flipped = points;
    for (auto& point : flipped) for (auto& h : point)
      if ((h.plane==1 && (mask&1)) || (h.plane==2 && (mask&2))) h.shower=1-h.shower;
    result.flipConflict[mask] = pointConsistency(flipped).fraction();
  }
  return result;
}

} }
#endif
