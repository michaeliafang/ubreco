#ifndef UBRECO_SHOWERRECO_VERTEXSPLIT_H
#define UBRECO_SHOWERRECO_VERTEXSPLIT_H

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <map>
#include <stdexcept>
#include <utility>
#include <vector>

namespace showerreco {
namespace recovery {

// Coordinates are in cm relative to the projected interaction vertex.
struct Hit2D { double w, t, charge; };
struct SplitConfig {
  std::size_t minHits = 10;
  std::size_t minUpstreamHits = 6;
  std::size_t maxHits = 6000;
  double neighborRadius = 0.8;
  double innerRadius = 0.;
  double upstreamRadius = 10.;
  double maxStartRadius = 5.;
  std::size_t minStartHits = 3;
  double minTrackSpan = 3.;
  double maxTrackGap = 1.5;
  double trackWidth = 0.6;
  double minTrackFraction = 0.6;
  std::size_t maxEdges = 200000;
  double angleStepDegrees = 2.;
  double minOpeningDegrees = 10.;
  double minChargeFraction = 0.20;
  double maxCutFraction = 0.15;

  void validate() const {
    if (minStartHits == 0 || maxEdges == 0 ||
        !std::isfinite(maxStartRadius) || maxStartRadius <= innerRadius || maxStartRadius > upstreamRadius ||
        !std::isfinite(minTrackSpan) || minTrackSpan <= 0. || minTrackSpan > upstreamRadius ||
        !std::isfinite(maxTrackGap) || maxTrackGap <= 0. ||
        !std::isfinite(trackWidth) || trackWidth <= 0. ||
        !std::isfinite(minTrackFraction) || minTrackFraction <= 0. || minTrackFraction > 1. ||
        minHits < 2 || minUpstreamHits == 0 || maxHits < 2 * minHits ||
        !std::isfinite(neighborRadius) || neighborRadius <= 0. ||
        !std::isfinite(innerRadius) || innerRadius < 0. ||
        !std::isfinite(upstreamRadius) || upstreamRadius <= innerRadius ||
        !std::isfinite(angleStepDegrees) || angleStepDegrees < 0.1 || angleStepDegrees > 90. ||
        !std::isfinite(minOpeningDegrees) || minOpeningDegrees <= 0. || minOpeningDegrees >= 180. ||
        !std::isfinite(minChargeFraction) || minChargeFraction < 0. || minChargeFraction > 0.5 ||
        !std::isfinite(maxCutFraction) || maxCutFraction < 0. || maxCutFraction > 1.)
      throw std::invalid_argument("Invalid single-shower recovery split configuration");
  }
};
struct Split {
  std::array<std::vector<std::size_t>, 2> parts;
  double cutFraction = std::numeric_limits<double>::infinity();
  double openingDegrees = 0.;
  bool valid() const { return !parts[0].empty() && !parts[1].empty(); }
};

struct SplitDiagnostics {
  std::size_t directions = 0, trackPairs = 0, balancedPairs = 0, connectedPairs = 0;
  bool edgeLimit = false;
};

struct EarlyTrack {
  double w = 0., t = 0., originW = 0., originT = 0.;
  bool valid = false;
};

// Weighted 2D line fit with a free offset. Weight is r^2/(r^2+1 cm^2),
// so the closest hits count as start evidence without dominating direction.
// The ray origin is the point on the fitted line closest to the vertex.
inline EarlyTrack fitEarlyTrack(const std::vector<Hit2D>& hits,
                               const std::vector<std::size_t>& indices) {
  EarlyTrack result;
  double sum = 0., x = 0., y = 0.;
  for (auto i : indices) {
    const auto& h = hits[i];
    const double r2 = h.w*h.w+h.t*h.t, weight = r2/(r2+1.);
    sum += weight; x += weight*h.w; y += weight*h.t;
  }
  if (sum <= 0.) return result;
  x /= sum; y /= sum;
  double xx = 0., xy = 0., yy = 0.;
  for (auto i : indices) {
    const auto& h = hits[i];
    const double r2 = h.w*h.w+h.t*h.t, weight = r2/(r2+1.);
    xx += weight*(h.w-x)*(h.w-x);
    xy += weight*(h.w-x)*(h.t-y);
    yy += weight*(h.t-y)*(h.t-y);
  }
  if (xx+yy <= 0.) return result;
  const double angle = 0.5*std::atan2(2.*xy,xx-yy);
  result.w = std::cos(angle); result.t = std::sin(angle);
  if (result.w*x+result.t*y < 0.) { result.w = -result.w; result.t = -result.t; }
  const double along = result.w*x+result.t*y;
  result.originW = x-along*result.w;
  result.originT = y-along*result.t;
  result.valid = true;
  return result;
}

inline EarlyTrack earlyTrack(const std::vector<Hit2D>& hits,
                            const std::vector<std::size_t>& upstream,
                            const std::vector<unsigned int>& side,
                            unsigned int label, const SplitConfig& cfg) {
  std::vector<std::size_t> indices;
  for (auto i : upstream) if (side[i] == label) indices.push_back(i);
  const auto total = indices.size();
  if (total < cfg.minUpstreamHits) return {};
  auto result = fitEarlyTrack(hits,indices);
  if (!result.valid) return result;
  // One bounded refinement removes off-axis branch/noise hits.
  indices.erase(std::remove_if(indices.begin(),indices.end(),[&](auto i) {
    const double w = hits[i].w-result.originW, t = hits[i].t-result.originT;
    return w*result.w+t*result.t <= 0. || std::abs(w*result.t-t*result.w) > cfg.trackWidth;
  }),indices.end());
  result = fitEarlyTrack(hits,indices);
  if (!result.valid) return result;
  result.valid = false;
  std::vector<std::pair<double,double>> support; // distance along ray, radius
  for (auto i : upstream) if (side[i] == label) {
    const auto& h = hits[i];
    const double w = h.w-result.originW, t = h.t-result.originT;
    const double along = w*result.w+t*result.t;
    if (along > 0. && std::abs(w*result.t-t*result.w) <= cfg.trackWidth)
      support.emplace_back(along,std::hypot(h.w,h.t));
  }
  std::sort(support.begin(),support.end());
  std::size_t begin = 0, near = 0;
  for (std::size_t j = 0; j < support.size(); ++j) {
    if (j && support[j].first-support[j-1].first > cfg.maxTrackGap) { begin = j; near = 0; }
    near += support[j].second <= cfg.maxStartRadius;
    const auto count = j-begin+1;
    if (count >= cfg.minUpstreamHits && near >= cfg.minStartHits &&
        double(count) >= cfg.minTrackFraction*total &&
        support[j].first-support[begin].first >= cfg.minTrackSpan) {
      result.valid = true;
      break;
    }
  }
  return result;
}

// A deliberately separate scoring graph: Euclidean neighbors, not the original
// RMS/time-overlap compatibility graph. No changes to proximity clustering.
inline Split splitAtVertex(const std::vector<Hit2D>& hits, const SplitConfig& cfg,
                           SplitDiagnostics* diagnostics = nullptr) {
  cfg.validate();
  SplitDiagnostics local;
  auto& diag = diagnostics ? *diagnostics : local;
  diag = {};
  Split best;
  if (hits.size() < 2 * cfg.minHits || hits.size() > cfg.maxHits) return best;
  using Cell = std::pair<int, int>;
  std::map<Cell, std::vector<std::size_t>> cells;
  std::vector<std::pair<std::size_t, std::size_t>> edges;
  std::vector<double> radius(hits.size());
  std::vector<std::size_t> upstream;
  double totalCharge = 0.;
  for (std::size_t i = 0; i < hits.size(); ++i) {
    const auto& h = hits[i];
    if (!std::isfinite(h.w) || !std::isfinite(h.t) || !std::isfinite(h.charge)) return best;
    // Guard conversion to integer cell coordinates for malformed input.
    if (std::abs(h.w / cfg.neighborRadius) > 1.e8 ||
        std::abs(h.t / cfg.neighborRadius) > 1.e8) return best;
    radius[i] = std::hypot(h.w,h.t);
    if (radius[i] > 0. && radius[i] >= cfg.innerRadius && radius[i] <= cfg.upstreamRadius)
      upstream.push_back(i);
    totalCharge += std::max(0., h.charge);
    Cell cell{int(std::floor(h.w / cfg.neighborRadius)), int(std::floor(h.t / cfg.neighborRadius))};
    for (int dw = -1; dw <= 1; ++dw) for (int dt = -1; dt <= 1; ++dt) {
      auto found = cells.find({cell.first + dw, cell.second + dt});
      if (found == cells.end()) continue;
      for (auto j : found->second) {
        const double x = h.w - hits[j].w, y = h.t - hits[j].t;
        if (x*x + y*y <= cfg.neighborRadius * cfg.neighborRadius) {
          // Two genuine starts can touch near the vertex. Do not penalize that
          // contact, but keep those hits in trajectory fitting and ownership.
          if (radius[i] <= cfg.maxStartRadius && radius[j] <= cfg.maxStartRadius) continue;
          if (edges.size() >= cfg.maxEdges) { diag.edgeLimit = true; return best; }
          edges.emplace_back(j, i);
        }
      }
    }
    cells[cell].push_back(i);
  }
  if (edges.empty() || totalCharge <= 0.) return best;
  constexpr double pi = 3.14159265358979323846;
  std::vector<unsigned int> side(hits.size());
  // Scan unoriented lines through the vertex. Neither side is required to have
  // half the hits or charge. Every hit belongs to exactly one side.
  for (double deg = 0.; deg < 180.; deg += cfg.angleStepDegrees) {
    const double nx = -std::sin(deg*pi/180.), ny = std::cos(deg*pi/180.);
    ++diag.directions;
    for (std::size_t i = 0; i < hits.size(); ++i)
      side[i] = (nx*hits[i].w + ny*hits[i].t >= 0.) ? 1 : 0;
    auto left = earlyTrack(hits,upstream,side,0,cfg);
    auto right = earlyTrack(hits,upstream,side,1,cfg);
    if (!left.valid || !right.valid) continue;
    // The fitted rays guide ALL hits, including downstream ownership.
    const auto distance = [](const Hit2D& h, const EarlyTrack& track) {
      const double w = h.w-track.originW, t = h.t-track.originT;
      const double along = w*track.w+t*track.t;
      const double across = w*track.t-t*track.w;
      return across*across + (along < 0. ? along*along : 0.);
    };
    for (std::size_t i = 0; i < hits.size(); ++i)
      side[i] = distance(hits[i],left) <= distance(hits[i],right) ? 0 : 1;
    // Validate sustained starts again after ownership changes. No iteration.
    left = earlyTrack(hits,upstream,side,0,cfg);
    right = earlyTrack(hits,upstream,side,1,cfg);
    if (!left.valid || !right.valid) continue;
    const double cosine = std::max(-1.,std::min(1.,left.w*right.w+left.t*right.t));
    const double opening = std::acos(cosine)*180./pi;
    if (opening < cfg.minOpeningDegrees) continue;
    ++diag.trackPairs;
    std::array<std::size_t, 2> count{{0,0}}, internal{{0,0}};
    std::array<double, 2> charge{{0.,0.}};
    for (std::size_t i = 0; i < hits.size(); ++i) {
      ++count[side[i]];
      charge[side[i]] += std::max(0.,hits[i].charge);
    }
    if (std::min(count[0], count[1]) < cfg.minHits ||
        std::min(charge[0], charge[1]) < cfg.minChargeFraction * totalCharge) continue;
    ++diag.balancedPairs;
    std::size_t cut = 0;
    for (const auto& edge : edges) {
      if (side[edge.first] == side[edge.second]) ++internal[side[edge.first]];
      else ++cut;
    }
    if (std::min(internal[0], internal[1]) == 0) continue;
    const double fraction = double(cut) / (2.*std::min(internal[0], internal[1]) + cut);
    if (fraction > cfg.maxCutFraction) continue;
    ++diag.connectedPairs;
    if (fraction > best.cutFraction ||
        (fraction == best.cutFraction && opening <= best.openingDegrees)) continue;
    best.parts = {};
    best.cutFraction = fraction;
    best.openingDegrees = opening;
    for (std::size_t i = 0; i < hits.size(); ++i) best.parts[side[i]].push_back(i);
  }
  return best;
}

} // namespace recovery
} // namespace showerreco
#endif
