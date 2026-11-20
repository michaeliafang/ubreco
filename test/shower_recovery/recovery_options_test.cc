#include "ubreco/ShowerReco/ShowerReco3D/Recovery/RecoveryEligibility.h"
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
  std::cout << "Recovery option tests passed\n";
}
