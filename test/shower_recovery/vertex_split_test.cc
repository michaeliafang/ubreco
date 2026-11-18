#include "ubreco/ShowerReco/ShowerReco3D/Recovery/VertexSplit.h"
#include <algorithm>
#include <cassert>
#include <iostream>
#include <random>
#include <set>

using namespace showerreco::recovery;
constexpr double pi = 3.14159265358979323846;

void ray(std::vector<Hit2D>& hits, double degrees, double start, double end) {
  const double angle = degrees*pi/180.;
  for (double r = start; r <= end; r += 0.4)
    hits.push_back({r*std::cos(angle),r*std::sin(angle),1.});
}
void partition(const Split& split, std::size_t n) {
  assert(split.valid());
  std::set<std::size_t> seen;
  for (const auto& part : split.parts)
    for (auto index : part) { assert(index < n); assert(seen.insert(index).second); }
  assert(seen.size() == n);
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
  for (std::size_t i=0; i<first; ++i) assert(label[i] == label[0]);
  for (std::size_t i=first; i<secondEnd; ++i) assert(label[i] != label[0]);

  std::mt19937 rng(42);
  auto shuffled = hits;
  std::shuffle(shuffled.begin(),shuffled.end(),rng);
  auto shuffledSplit = splitAtVertex(shuffled,cfg);
  partition(shuffledSplit,shuffled.size());
  assert(std::abs(shuffledSplit.cutFraction-split.cutFraction) < 1.e-12);

  // Rotation across the angular wrap-around must still allow a valid split.
  for (auto& h : hits) { h.w=-h.w; h.t=-h.t; }
  partition(splitAtVertex(hits,cfg),hits.size());

  std::vector<Hit2D> straight;
  ray(straight,0.,6.,40.);
  assert(!splitAtVertex(straight,cfg).valid());
  // A late detached branch has no independent support in the upstream annulus.
  ray(straight,35.,25.,40.);
  assert(!splitAtVertex(straight,cfg).valid());
  assert(!splitAtVertex({},cfg).valid());
  assert(!splitAtVertex({{1.,1.,1.}},cfg).valid());
  auto invalid = hits;
  invalid[0].w = std::numeric_limits<double>::quiet_NaN();
  assert(!splitAtVertex(invalid,cfg).valid());
  auto limited = cfg;
  limited.maxHits = 2*cfg.minHits;
  assert(!splitAtVertex(hits,limited).valid());
  auto strict = cfg;
  strict.maxCutFraction = 0.;
  assert(!splitAtVertex(hits,strict).valid());
  bool threw=false;
  try { auto bad=cfg; bad.neighborRadius=0.; splitAtVertex(hits,bad); }
  catch (const std::invalid_argument&) { threw=true; }
  assert(threw);
  std::cout << "Vertex splitter tests passed\n";
}
