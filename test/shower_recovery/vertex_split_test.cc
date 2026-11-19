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
  ray(hits,-20.,1.,35.);
  const std::size_t first = hits.size();
  ray(hits,20.,1.,22.); // Unequal shower sizes are intentional.
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
  ray(unequalCharge,-25.,1.,18.);
  const auto boundary = unequalCharge.size();
  ray(unequalCharge,25.,1.,18.);
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
  oldRadius.maxStartRadius = 15.;
  CHECK(splitAtVertex(distant,oldRadius).valid());
  // A valid pair touching at the vertex must not fail because of local edges.
  std::vector<Hit2D> touching;
  ray(touching,-12.,0.,30.);
  ray(touching,12.,0.,30.);
  SplitDiagnostics diagnostic;
  partition(splitAtVertex(touching,cfg,&diagnostic),touching.size());
  CHECK(diagnostic.trackPairs > 0);

  // Starts just beyond 5 cm are rejected even though they extend through 10 cm.
  std::vector<Hit2D> late;
  ray(late,-20.,5.2,25.);
  ray(late,20.,5.2,25.);
  CHECK(!splitAtVertex(late,cfg).valid());
  auto relaxedStart = cfg;
  relaxedStart.maxStartRadius = 7.;
  CHECK(splitAtVertex(late,relaxedStart).valid());

  // Three near-vertex hits plus distant rays are not two sustained starts.
  std::vector<Hit2D> fragments;
  for (double angle : {-20.,20.}) {
    ray(fragments,angle,1.,1.8);
    ray(fragments,angle,7.,30.);
  }
  CHECK(!splitAtVertex(fragments,cfg).valid());

  // Modest transverse jitter and short local gaps retain track evidence.
  auto irregular = touching;
  for (auto& h : irregular) h.t += 0.06*std::sin(7.*h.w);
  irregular.erase(std::remove_if(irregular.begin(),irregular.end(),[](const auto& h) {
    const double r = std::hypot(h.w,h.t);
    return r > 6. && r < 6.7;
  }),irregular.end());
  partition(splitAtVertex(irregular,cfg),irregular.size());

  // A modest vertex displacement need not force both fits through the origin.
  auto shifted = touching;
  for (auto& h : shifted) h.t -= 2.;
  partition(splitAtVertex(shifted,cfg),shifted.size());

  // Missing-hit band wider than the old gap cut retains both starts.
  auto deadBand = touching;
  deadBand.erase(std::remove_if(deadBand.begin(),deadBand.end(),[](const auto& h) {
    const double r = std::hypot(h.w,h.t);
    return r > 3. && r < 5.2;
  }),deadBand.end());
  partition(splitAtVertex(deadBand,cfg),deadBand.size());

  // Off-axis activity no longer dilutes a clear line's supporting-hit fraction.
  std::vector<Hit2D> branchy;
  for (double x=0.5; x<9.; x+=0.4) {
    branchy.push_back({x,0.,1.});
    branchy.push_back({x,2.,1.});
    branchy.push_back({x,-2.,1.});
  }
  std::vector<std::size_t> early;
  for (std::size_t i=0; i<branchy.size(); ++i) early.push_back(i);
  std::vector<unsigned int> sides(branchy.size(),0);
  CHECK(earlyTrack(branchy,early,sides,0,cfg).valid);
  auto oldFraction = cfg; oldFraction.minTrackFraction = 0.6;
  CHECK(!earlyTrack(branchy,early,sides,0,oldFraction).valid);

  // Direct density tests isolate ranking from the start-fitting and graph cuts.
  EarlyTrack a{std::cos(20.*pi/180.),-std::sin(20.*pi/180.),0.,0.,true};
  EarlyTrack b{a.w,-a.t,0.,0.,true};
  std::vector<Hit2D> valley;
  ray(valley,-20.,1.,12.);
  ray(valley,20.,1.,12.);
  const auto clear = valleyEvidence(valley,a,b,cfg);
  CHECK(clear.bands == 4 && clear.score > 0.99);
  auto crossing = valley;
  for (double y=-2.; y<=2.; y+=0.1) crossing.push_back({5.,y,1.});
  const auto bridged = valleyEvidence(crossing,a,b,cfg);
  CHECK(bridged.score > 0. && bridged.score < clear.score);
  auto filled = valley;
  for (double x=2.; x<10.; x+=0.1) filled.push_back({x,0.,1.});
  CHECK(valleyEvidence(filled,a,b,cfg).score < bridged.score);
  std::vector<Hit2D> oneSide;
  ray(oneSide,-20.,1.,12.);
  CHECK(valleyEvidence(oneSide,a,b,cfg).score == 0.);
  CHECK(valleyEvidence({},a,b,cfg).score == 0.);
  auto rotated = valley;
  for (auto& h : rotated) { const auto w=h.w; h.w=-h.t; h.t=w; }
  EarlyTrack ar{-a.t,a.w,0.,0.,true}, br{-b.t,b.w,0.,0.,true};
  CHECK(std::abs(valleyEvidence(rotated,ar,br,cfg).score-clear.score) < 1.e-12);
  Split preferred, fallback;
  preferred.parts = {{{0},{1}}}; fallback.parts = preferred.parts;
  preferred.valleyScore = clear.score; preferred.cutFraction = 0.1;
  fallback.valleyScore = bridged.score; fallback.cutFraction = 0.01;
  CHECK(betterSplit(preferred,fallback));
  preferred.valleyScore = fallback.valleyScore = 0.;
  CHECK(betterSplit(fallback,preferred)); // absent valley evidence falls back to cut

  // Dense inputs have a bounded graph allocation and an explicit rejection.
  auto bounded = cfg;
  bounded.maxEdges = 1;
  CHECK(!splitAtVertex(touching,bounded,&diagnostic).valid());
  CHECK(diagnostic.edgeLimit);
  threw = false;
  try { auto bad=cfg; bad.maxStartRadius=11.; splitAtVertex(touching,bad); }
  catch (const std::invalid_argument&) { threw=true; }
  CHECK(threw);
  std::cout << "Vertex splitter tests passed\n";
}
