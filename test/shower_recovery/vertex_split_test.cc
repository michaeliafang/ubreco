#include "ubreco/ShowerReco/ShowerReco3D/Recovery/VertexSplit.h"
#include <algorithm>
#include <cstdlib>
#include <iostream>
#include <random>
#include <set>

// Unlike assert(), test checks must run when prof/release builds define NDEBUG.
#define CHECK(condition) do { \
  if (!(condition)) { \
    std::cerr << __FILE__ << ':' << __LINE__ << ": check failed: " \
              << #condition << '\n'; \
    std::exit(EXIT_FAILURE); \
  } \
} while (false)

using namespace showerreco::recovery;
constexpr double pi = 3.14159265358979323846;

void ray(std::vector<Hit2D>& hits, double degrees, double start, double end) {
  const double angle = degrees*pi/180.;
  for (double r = start; r <= end; r += 0.4)
    hits.push_back({r*std::cos(angle),r*std::sin(angle),1.});
}
void partition(const Split& split, std::size_t n) {
  CHECK(split.valid());
  std::set<std::size_t> seen;
  for (const auto& part : split.parts)
    for (auto index : part) { CHECK(index < n); CHECK(seen.insert(index).second); }
  CHECK(seen.size() == n);
}
int main() {
  SplitConfig cfg;
  std::vector<Hit2D> hits;
  ray(hits,-20.,6.,35.);
  const std::size_t first = hits.size();
  ray(hits,20.,6.,22.); // Unequal shower sizes are intentional.
  const std::size_t secondEnd = hits.size();
  // A narrow bridge between the rays, away from the vertex.
  const double x = 15.*std::cos(20.*pi/180.);
  for (double y = -15.*std::sin(20.*pi/180.); y <= 15.*std::sin(20.*pi/180.); y += 0.4)
    hits.push_back({x,y,1.});
  auto split = splitAtVertex(hits,cfg);
  partition(split,hits.size());
  std::vector<int> label(hits.size());
  for (int s = 0; s < 2; ++s) for (auto i : split.parts[s]) label[i]=s;
  for (std::size_t i=0; i<first; ++i) CHECK(label[i] == label[0]);
  for (std::size_t i=first; i<secondEnd; ++i) CHECK(label[i] != label[0]);

  std::mt19937 rng(42);
  auto shuffled = hits;
  std::shuffle(shuffled.begin(),shuffled.end(),rng);
  auto shuffledSplit = splitAtVertex(shuffled,cfg);
  partition(shuffledSplit,shuffled.size());
  CHECK(std::abs(shuffledSplit.cutFraction-split.cutFraction) < 1.e-12);

  // Rotation across the angular wrap-around must still allow a valid split.
  for (auto& h : hits) { h.w=-h.w; h.t=-h.t; }
  partition(splitAtVertex(hits,cfg),hits.size());

  std::vector<Hit2D> straight;
  ray(straight,0.,6.,40.);
  CHECK(!splitAtVertex(straight,cfg).valid());
  // A late detached branch has no independent support in the upstream annulus.
  ray(straight,35.,25.,40.);
  CHECK(!splitAtVertex(straight,cfg).valid());
  CHECK(!splitAtVertex({},cfg).valid());
  CHECK(!splitAtVertex({{1.,1.,1.}},cfg).valid());
  auto invalid = hits;
  invalid[0].w = std::numeric_limits<double>::quiet_NaN();
  CHECK(!splitAtVertex(invalid,cfg).valid());
  auto limited = cfg;
  limited.maxHits = 2*cfg.minHits;
  CHECK(!splitAtVertex(hits,limited).valid());
  auto strict = cfg;
  strict.maxCutFraction = 0.;
  CHECK(!splitAtVertex(hits,strict).valid());
  bool threw=false;
  try { auto bad=cfg; bad.neighborRadius=0.; splitAtVertex(hits,bad); }
  catch (const std::invalid_argument&) { threw=true; }
  CHECK(threw);
  // Equal hit counts, unequal charge: the 20% cut is a charge requirement.
  std::vector<Hit2D> unequalCharge;
  ray(unequalCharge,-25.,6.,18.);
  const auto boundary = unequalCharge.size();
  ray(unequalCharge,25.,6.,18.);
  for (std::size_t i=boundary; i<unequalCharge.size(); ++i)
    unequalCharge[i].charge = 0.19/0.81;
  CHECK(!splitAtVertex(unequalCharge,cfg).valid());
  for (std::size_t i=boundary; i<unequalCharge.size(); ++i)
    unequalCharge[i].charge = 0.21/0.79;
  CHECK(splitAtVertex(unequalCharge,cfg).valid());

  // A balanced pair starting beyond 10 cm must fail the tightened support cut.
  std::vector<Hit2D> distant;
  ray(distant,-25.,12.,25.);
  ray(distant,25.,12.,25.);
  CHECK(!splitAtVertex(distant,cfg).valid());
  auto oldRadius = cfg;
  oldRadius.upstreamRadius = 20.;
  CHECK(splitAtVertex(distant,oldRadius).valid());
  std::cout << "Vertex splitter tests passed\n";
}
