#include "ubreco/ShowerReco/ShowerReco3D/Recovery/RecoveryEligibility.h"
#include "ubreco/ShowerReco/ShowerReco3D/Recovery/MatchHypotheses.h"
#include <cstdlib>
#include <iostream>
#include <limits>
#define CHECK(c) do { if (!(c)) { std::cerr << "Failed line " << __LINE__ << ": " << #c << '\n'; std::exit(EXIT_FAILURE); } } while(false)
int main() {
  using namespace showerreco::recovery;
  CHECK(recoveryTarget({2.,50.},{80.,20.},20.,0.3)==0);
  CHECK(recoveryTarget({50.,2.},{20.,80.},20.,0.3)==1);
  CHECK(recoveryTarget({2.,20.},{80.,20.},20.,0.3)==-1);
  CHECK(recoveryTarget({2.,50.},{70.,30.},20.,0.3)==-1);
  CHECK(recoveryTarget({40.,50.},{80.,20.},20.,0.3)==-1);
  CHECK(recoveryTarget({2.,3.},{80.,20.},20.,0.3)==-1);
  CHECK(recoveryTarget({2.,50.},{0.,20.},20.,0.3)==-1);
  CHECK(recoveryTarget({2.,50.},{80.,std::numeric_limits<double>::quiet_NaN()},20.,0.3)==-1);
  const ClusterMatches base{{0,1,2},{3,4,5}};
  const std::vector<unsigned> planes{0,1,2,0,1,2};
  const auto variants=pairingHypotheses(base,planes);
  CHECK(variants.size()==4 && variants[0].matches==base);
  CHECK((variants[2].matches==ClusterMatches{{0,1,5},{3,4,2}}));
  CHECK(pairingHypotheses({{1,2},{4,5}},planes).size()==2);
  CHECK(pairingHypotheses({{0},{3}},planes).empty());
  CHECK(pairingHypotheses({{0,1},{3,5}},planes).empty());
  CHECK(pairingHypotheses({{0,1},{0,4}},planes).empty());
  CHECK(pairingHypotheses({{0,3},{1,4}},planes).empty());
  CHECK(pairingHypotheses({{0,1},{3,4},{2,5}},planes).empty());
  std::cout << "Recovery option tests passed\n";
}
