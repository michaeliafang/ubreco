#include "ubreco/ShowerReco/ShowerReco3D/Recovery/HitOwnership.h"
#include <cstdlib>
#include <iostream>
#include <utility>
#define CHECK(c) do { if (!(c)) { std::cerr << "Failed line " << __LINE__ << ": " << #c << '\n'; std::exit(EXIT_FAILURE); } } while(false)
int main() {
  using namespace showerreco::recovery;
  using Hit = std::pair<unsigned,unsigned>; // product ID, key
  const Hit a{1,0}, b{1,1}, c{1,2}, otherProduct{2,0};
  const std::vector<Hit> first{a,b,a,b}, second{b,c}, independent{otherProduct};
  const auto clean = uniqueHitsInOrder(first);
  CHECK((clean == std::vector<Hit>{a,b}));
  CHECK((uniqueHitsInOrder(second) == second)); // shared b remains in both candidates
  CHECK(!disjointHits(clean,second)); // rejected if both become output showers
  CHECK(disjointHits(clean,independent)); // unused overlapping candidate is irrelevant
  CHECK(disjointHits(std::vector<Hit>{a},independent)); // same key, different product
  CHECK(uniqueHitsInOrder(std::vector<Hit>{}).empty());
  CHECK(disjointHits(clean,std::vector<Hit>{}));
  CHECK(disjointHits(std::vector<Hit>{a,b},std::vector<Hit>{c}));
  // Same input grouping used by combine/split and fragment/parent attachment.
  const auto grouped=hitsByPlane(std::vector<Hit>{a,b,c,otherProduct,a},
    [](const Hit& h) { return h.second%3; });
  CHECK((grouped[0]==std::vector<Hit>{a,otherProduct}));
  CHECK((grouped[1]==std::vector<Hit>{b}));
  CHECK((grouped[2]==std::vector<Hit>{c}));
  std::vector<Hit> combined;
  for (const auto& plane : grouped) combined.insert(combined.end(),plane.begin(),plane.end());
  CHECK(combined.size()==4); // detached-fragment hits and product identity retained
  CHECK(disjointHits(grouped[0],grouped[1]));
  // Two existing clusters in plane 1 must survive an attachment in plane 2.
  const auto planeOf=[](int h) { return unsigned(h/100); };
  std::vector<std::vector<int>> clusters{{100,101},{110,111},{200,201},{210},{220}};
  const auto before=clusters;
  CHECK(attachFragment(clusters,std::vector<int>{100,101,200,201},
    std::vector<int>{220,220},planeOf));
  CHECK(clusters[0]==before[0] && clusters[1]==before[1]);
  CHECK((clusters[2]==std::vector<int>{200,201,220}));
  CHECK(clusters[3]==before[3] && clusters[4].empty());
  std::set<int> allBefore, allAfter;
  for (const auto& c0 : before) allBefore.insert(c0.begin(),c0.end());
  for (const auto& c0 : clusters) allAfter.insert(c0.begin(),c0.end());
  CHECK(allBefore==allAfter);
  CHECK(disjointHits(clusters[2],clusters[3]));
  // An alternative parent gets the fragment without changing either other view.
  auto alternate=before;
  CHECK(attachFragment(alternate,std::vector<int>{110,111,210},std::vector<int>{220},planeOf));
  CHECK((alternate[3]==std::vector<int>{210,220}));
  CHECK(alternate[2]==before[2]);
  // Missing parent evidence in a fragment's plane must not partly mutate inputs.
  auto unsupported=before;
  CHECK(!attachFragment(unsupported,std::vector<int>{100,101},std::vector<int>{220},planeOf));
  CHECK(unsupported==before);
  // Union input includes unmatched cmerger hits, not just baseline shower hits.
  std::vector<int> unionHits;
  for (const auto& c0 : before) unionHits.insert(unionHits.end(),c0.begin(),c0.end());
  const auto unionPlanes=hitsByPlane(unionHits,planeOf);
  CHECK((unionPlanes[1]==std::vector<int>{100,101,110,111}));
  CHECK((unionPlanes[2]==std::vector<int>{200,201,210,220}));
  std::vector<int> reference;
  std::set<int> selected;
  for (int h=100; h<200; ++h) { reference.push_back(h); if (h<198) selected.insert(h); }
  reference.push_back(100); // duplicate association must not change denominator
  reference.push_back(200); // omitted plane must not enter the denominator
  selected.insert(999); // non-reference hit must not inflate numerator
  auto coverage=retainedPlaneCoverage(reference,selected,2u,planeOf);
  CHECK(coverage.retained==98 && coverage.total==100 && coverage.passes(.98));
  selected.erase(197);
  coverage=retainedPlaneCoverage(reference,selected,2u,planeOf);
  CHECK(!coverage.passes(.98));
  selected.insert(197); // adding the fragment makes attachment worth trying
  CHECK(retainedPlaneCoverage(reference,selected,2u,planeOf).passes(.98));
  CHECK(!retainedPlaneCoverage(reference,selected,6u,planeOf).passes(.98));
  CHECK(!retainedPlaneCoverage(reference,selected,1u,planeOf).passes(.98));
  CHECK(!retainedPlaneCoverage(std::vector<int>{},selected,7u,planeOf).passes(.98));
  // A valid below-target pair is retained; only a better valid result replaces it.
  HitCoverage saved{97,100};
  CHECK(!saved.passes(.98));
  const HitCoverage worse{95,100}, better{99,100}, empty{};
  CHECK(!(worse.fraction()>saved.fraction()));
  CHECK(!(empty.fraction()>saved.fraction()));
  CHECK(better.fraction()>saved.fraction() && better.passes(.98));
  std::cout << "Hit ownership tests passed\n";
}
