// Compile with SPLITTER_HEADER set to an older header to compare implementations.
#ifndef SPLITTER_HEADER
#define SPLITTER_HEADER "ubreco/ShowerReco/ShowerReco3D/Recovery/VertexSplit.h"
#endif
#include SPLITTER_HEADER
#include <chrono>
#include <iostream>

int main() {
  using namespace showerreco::recovery;
  for (const std::size_t n : {300u,1000u,3000u,6000u}) {
    std::vector<Hit2D> hits;
    for (std::size_t i=0; i<n; ++i) {
      const double r = 0.1+80.*(i/2)/double(n/2);
      const double angle = (i%2 ? 0.22 : -0.22);
      // Forward branching beyond the two early starts, plus small variation.
      const double branch = r>12. ? 0.04*std::sin(r*0.3) : 0.;
      hits.push_back({r*std::cos(angle+branch),r*std::sin(angle+branch)+0.03*std::sin(i),1.});
    }
    std::vector<double> times;
    std::size_t accepted=0;
    for (int repeat=0; repeat<7; ++repeat) {
      const auto start=std::chrono::steady_clock::now();
      const auto result=splitAtVertex(hits,SplitConfig{});
      const double ms=std::chrono::duration<double,std::milli>(std::chrono::steady_clock::now()-start).count();
      if (repeat) times.push_back(ms); // warmup excluded
      accepted += result.valid();
    }
    std::sort(times.begin(),times.end());
    std::cout << n << " hits: median_ms=" << (times[2]+times[3])/2.
              << " accepted=" << accepted << "/7\n";
  }
}
