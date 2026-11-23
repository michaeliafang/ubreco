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
  const auto fallback=twoPlaneHypotheses(base,planes,true);
  CHECK(fallback.size()==4);
  for (const auto& h : fallback) {
    CHECK(h.twoPlaneFallback && (h.planeMask==5 || h.planeMask==6));
    CHECK(h.matches[0].size()==2 && h.matches[1].size()==2);
    for (const auto& match : h.matches) {
      bool collection=false;
      for (auto ci : match) collection=collection || planes[ci]==2;
      CHECK(collection);
    }
  }
  CHECK(twoPlaneHypotheses(base,planes,false).size()==2);
  CHECK(twoPlaneHypotheses({{0,1},{3,4}},planes,true).empty());
  CHECK(twoPlaneHypotheses({{0,2},{3,5}},planes,true).empty());
  CHECK(twoPlaneHypotheses({{0,1,2},{3,4,5},{6}},planes,true).empty());
  CHECK(detachedFragment({2.,50.},{75.,25.},20.,.30)==1);
  CHECK(detachedFragment({2.,50.},{70.,30.},20.,.30)==-1);
  CHECK(detachedFragment({2.,3.,50.},{45.,36.,19.},20.,.20)==2);
  CHECK(detachedFragment({2.,3.,50.},{45.,35.,20.},20.,.20)==-1);
  CHECK(detachedFragment({2.,30.,50.},{45.,36.,19.},20.,.20)==-1);
  CHECK(detachedFragment({2.,3.,20.},{45.,36.,19.},20.,.20)==-1);
  CHECK(detachedFragment({2.,3.,50.},{45.,0.,19.},20.,.20)==-1);
  CHECK(detachedFragment({2.,3.,50.},{45.,36.,std::numeric_limits<double>::quiet_NaN()},20.,.20)==-1);
  CHECK(detachedFragment({2.,3.},{45.},20.,.20)==-1);
  CHECK(attachmentAngle({1.,0.,0.},{50.,0.,0.})==0.);
  CHECK(std::abs(attachmentAngle({1.,0.,0.},{0.,50.,0.})-90.)<1.e-8);
  CHECK(attachmentAngle({1.,0.,0.},{-50.,0.,0.})==180.);
  CHECK(attachmentAngle({0.,0.,0.},{50.,0.,0.})==180.);
  RecoveryBudget budget;
  CHECK(budget.canTry(4,4));
  budget.trials=1; budget.showers=2;
  CHECK(budget.canTry(4,4) && budget.canReconstruct(2,4));
  CHECK(!budget.canReconstruct(3,4));
  budget.trials=2; budget.showers=4;
  CHECK(!budget.canTry(4,4) && !budget.canReconstruct(2,4));
  budget.showers=5;
  CHECK(!budget.canTry(4,4) && !budget.canReconstruct(2,4));
  // 4827-like case: two large clusters in each retained plane plus a small
  // extra group. Search is independent of the matcher's three-group result.
  const std::vector<unsigned> searchPlanes{1,1,1,2,2};
  auto searched=substantialPairHypotheses(searchPlanes,{201,740,12,450,471},
    {201.,740.,12.,450.,471.},10,3,.8);
  CHECK(searched.size()==2);
  CHECK((searched[0].matches==ClusterMatches{{1,4},{0,3}}));
  CHECK((searched[1].matches==ClusterMatches{{1,3},{0,4}}));
  CHECK(searched[0].selectedHits==1862);
  // A genuinely substantial third cluster cannot be silently dropped.
  CHECK(substantialPairHypotheses(searchPlanes,{300,300,300,450,471},
    {300.,300.,300.,450.,471.},10,3,.8).empty());
  // The charge denominator must include clusters beyond the top-N cap.
  CHECK(substantialPairHypotheses({1,1,1,1,2,2},{400,390,380,370,450,471},
    {400.,390.,380.,370.,450.,471.},10,2,.8).empty());
  CHECK(substantialPairHypotheses(searchPlanes,{201,740,12,450,471},
    {201.,740.,12.,0.,0.},10,3,.8).empty());
  // Search bound with three eligible clusters in each of three views.
  const auto bounded=substantialPairHypotheses({0,0,0,1,1,1,2,2,2},
    {100,100,100,100,100,100,100,100,100},
    {100.,100.,100.,100.,100.,100.,100.,100.,100.},10,3,.5);
  CHECK(bounded.size()==36);
  for (const auto& h : bounded) {
    CHECK(h.clusterSearch && h.twoPlaneFallback && h.matches.size()==2);
    std::set<unsigned> unique;
    for (const auto& match : h.matches) for (auto ci : match) CHECK(unique.insert(ci).second);
  }
  std::vector<std::array<std::size_t,3>> choices{{{0,0,1}},{{0,1,0}},{{0,1,1}}};
  promoteMultiPlaneRetry(choices);
  CHECK((choices[0]==std::array<std::size_t,3>{{0,0,1}}));
  CHECK((choices[1]==std::array<std::size_t,3>{{0,1,1}}));
  CHECK((choices[2]==std::array<std::size_t,3>{{0,1,0}}));
  choices={{{0,0,1}},{{0,1,0}}};
  promoteMultiPlaneRetry(choices);
  CHECK((choices[1]==std::array<std::size_t,3>{{0,1,0}}));
  choices.clear(); promoteMultiPlaneRetry(choices); CHECK(choices.empty());
  std::cout << "Recovery option tests passed\n";
}
