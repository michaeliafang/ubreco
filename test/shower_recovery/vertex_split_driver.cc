// Text adapter for synthetic-pattern studies: invokes the production splitter.
#include "ubreco/ShowerReco/ShowerReco3D/Recovery/VertexSplit.h"
#include <iomanip>
#include <iostream>
#include <vector>
int main() {
  std::size_t n;
  if (!(std::cin >> n) || n > 100000) return 1;
  std::vector<showerreco::recovery::Hit2D> hits(n);
  for (auto& h : hits) if (!(std::cin >> h.w >> h.t >> h.charge)) return 2;
  auto split = showerreco::recovery::splitAtVertex(hits, {});
  std::vector<int> labels(n, -1);
  if (split.valid()) {
    for (int s = 0; s < 2; ++s) for (auto i : split.parts[s]) {
      if (i >= n || labels[i] != -1) return 3;
      labels[i] = s;
    }
    for (auto s : labels) if (s == -1) return 4;
  }
  std::cout << split.valid() << ' ' << std::setprecision(17)
            << split.cutFraction << ' ' << split.openingDegrees << '\n';
  for (auto s : labels) std::cout << s << '\n';
}
