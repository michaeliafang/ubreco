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
  std::cout << "Hit ownership tests passed\n";
}
