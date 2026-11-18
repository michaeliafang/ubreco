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
  std::size_t minUpstreamHits = 3;
  std::size_t maxHits = 6000;
  double neighborRadius = 0.8;
  double innerRadius = 5.;
  double upstreamRadius = 20.;
  double angleStepDegrees = 2.;
  double minOpeningDegrees = 10.;
  double minChargeFraction = 0.05;
  double maxCutFraction = 0.15;

  void validate() const {
    if (minHits < 2 || minUpstreamHits == 0 || maxHits < 2 * minHits ||
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

// A deliberately separate scoring graph: Euclidean neighbors, not the original
// RMS/time-overlap compatibility graph. No changes to proximity clustering.
inline Split splitAtVertex(const std::vector<Hit2D>& hits, const SplitConfig& cfg) {
  cfg.validate();
  Split best;
  if (hits.size() < 2 * cfg.minHits || hits.size() > cfg.maxHits) return best;
  using Cell = std::pair<int, int>;
  std::map<Cell, std::vector<std::size_t>> cells;
  std::vector<std::pair<std::size_t, std::size_t>> edges;
  double totalCharge = 0.;
  for (std::size_t i = 0; i < hits.size(); ++i) {
    const auto& h = hits[i];
    if (!std::isfinite(h.w) || !std::isfinite(h.t) || !std::isfinite(h.charge)) return best;
    // Guard conversion to integer cell coordinates for malformed input.
    if (std::abs(h.w / cfg.neighborRadius) > 1.e8 ||
        std::abs(h.t / cfg.neighborRadius) > 1.e8) return best;
    totalCharge += std::max(0., h.charge);
    Cell cell{int(std::floor(h.w / cfg.neighborRadius)), int(std::floor(h.t / cfg.neighborRadius))};
    for (int dw = -1; dw <= 1; ++dw) for (int dt = -1; dt <= 1; ++dt) {
      auto found = cells.find({cell.first + dw, cell.second + dt});
      if (found == cells.end()) continue;
      for (auto j : found->second) {
        const double x = h.w - hits[j].w, y = h.t - hits[j].t;
        if (x*x + y*y <= cfg.neighborRadius * cfg.neighborRadius) edges.emplace_back(j, i);
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
    std::array<std::size_t, 2> count{{0,0}}, nearCount{{0,0}}, internal{{0,0}};
    std::array<double, 2> charge{{0.,0.}}, ux{{0.,0.}}, uy{{0.,0.}};
    for (std::size_t i = 0; i < hits.size(); ++i) {
      const auto& h = hits[i];
      const unsigned int s = (nx*h.w + ny*h.t >= 0.) ? 1 : 0;
      side[i] = s;
      ++count[s];
      charge[s] += std::max(0., h.charge);
      const double r = std::hypot(h.w, h.t);
      if (r >= cfg.innerRadius && r <= cfg.upstreamRadius && r > 0.) {
        ++nearCount[s];
        ux[s] += h.w/r;
        uy[s] += h.t/r;
      }
    }
    if (std::min(count[0], count[1]) < cfg.minHits ||
        std::min(nearCount[0], nearCount[1]) < cfg.minUpstreamHits ||
        std::min(charge[0], charge[1]) < cfg.minChargeFraction * totalCharge) continue;
    const double norm = std::hypot(ux[0],uy[0]) * std::hypot(ux[1],uy[1]);
    if (norm <= 0.) continue;
    const double cosine = std::max(-1., std::min(1., (ux[0]*ux[1]+uy[0]*uy[1])/norm));
    const double opening = std::acos(cosine)*180./pi;
    if (opening < cfg.minOpeningDegrees) continue;
    std::size_t cut = 0;
    for (const auto& edge : edges) {
      if (side[edge.first] == side[edge.second]) ++internal[side[edge.first]];
      else ++cut;
    }
    if (std::min(internal[0], internal[1]) == 0) continue;
    const double fraction = double(cut) / (2.*std::min(internal[0], internal[1]) + cut);
    if (fraction > cfg.maxCutFraction) continue;
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
