#include "ubreco/ShowerReco/ShowerReco3D/Recovery/ExactHitMapping.h"
#include "ubreco/ShowerReco/ShowerReco3D/Recovery/SpacePointConsistency.h"
#include <cstdlib>
#include <cmath>
#include <iostream>
#define CHECK(c) do { if (!(c)) { std::cerr << "Failed line " << __LINE__ << ": " << #c << '\n'; std::exit(EXIT_FAILURE); } } while(false)
int main() {
  using namespace showerreco::recovery;
  HitFingerprint a{}, b{}, c{}; a[0]=1; b[0]=2; c[0]=3;
  auto m=mapCopiedHits({a,b,c},{c,a});
  CHECK(m.matched==2 && m.missing==0);
  CHECK(m.sourceToSelected[0]==1 && m.sourceToSelected[1]==noHit && m.sourceToSelected[2]==0);
  m=mapCopiedHits({a,a,b},{a,b});
  CHECK(m.ambiguous==1 && m.matched==1);
  m=mapCopiedHits({a,b},{a,a,c});
  CHECK(m.ambiguous==2 && m.missing==1 && m.matched==0);
  c[1]=std::numeric_limits<double>::quiet_NaN();
  m=mapCopiedHits({a},{c}); CHECK(m.invalid==1);
  b=a; b[10]=1; // Same position/charge, different fit uncertainty: not identical.
  m=mapCopiedHits({a},{b}); CHECK(m.missing==1);
  auto score=pointConsistency({{{0,0,0},{1,1,0}},{{2,0,1},{3,2,1}}});
  CHECK(score.informative==2 && score.fraction()==0.);
  score=pointConsistency({{{0,0,0},{1,1,1}},{{2,0,1},{3,2,1}}});
  CHECK(score.fraction()==0.5);
  score=pointConsistency({{{0,0,0},{1,0,1}}}); CHECK(score.informative==0);
  score=pointConsistency({{{0,0,0},{1,1,1}},{{0,0,0},{1,1,1}}});
  CHECK(score.weight==1. && score.fraction()==1.); // repeated shared hits do not multiply weight
  // Global shower renaming cannot change conflict. Re-pairing one plane can.
  std::vector<std::vector<AssignedPointHit>> paired{
    {{0,0,0},{1,1,1}}, {{2,0,0},{3,1,1}}, {{4,0,0},{5,1,0}}};
  const auto originalConflict = pointConsistency(paired).fraction();
  for (auto& point : paired) for (auto& hit : point) hit.shower = 1-hit.shower;
  CHECK(pointConsistency(paired).fraction() == originalConflict);
  for (auto& point : paired) for (auto& hit : point)
    if (hit.plane == 1) hit.shower = 1-hit.shower;
  CHECK(std::abs(pointConsistency(paired).fraction()-(1.-originalConflict)) < 1.e-12);
  const std::vector<std::vector<AssignedPointHit>> evidence{
    {{0,0,0},{1,1,1},{2,2,0}}, {{3,0,1},{4,1,0},{5,2,1}}};
  const auto detail = diagnosePoints(evidence,true);
  CHECK(detail.pairs[0].fraction()==1. && detail.pairs[1].fraction()==0. && detail.pairs[2].fraction()==1.);
  CHECK(detail.unique.informative==2 && detail.unique.fraction()==1.);
  CHECK(detail.flipConflict[0]==1. && detail.flipConflict[1]==0.);
  auto repeated = evidence; repeated.push_back(evidence[0]);
  CHECK(diagnosePoints(repeated,true).unique.informative==1);
  CHECK(diagnosePoints({},true).unique.weight==0.);
  CHECK(pointConsistency(evidence).fraction()==1.); // diagnostics leave input untouched
  std::cout << "Spacepoint mapping and consistency tests passed\n";
}
