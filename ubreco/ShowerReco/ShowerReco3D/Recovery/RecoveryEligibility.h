#ifndef UBRECO_RECOVERY_ELIGIBILITY_H
#define UBRECO_RECOVERY_ELIGIBILITY_H
#include <array>
#include <algorithm>
#include <cmath>
#include <vector>
#include <cstddef>
namespace showerreco { namespace recovery {
// Return the near-vertex shower to recover; the other is the detached fragment.
// Charge is positive collection-plane integral, deduplicated before calling.
inline int recoveryTarget(const std::array<double,2>& distances,
                          const std::array<double,2>& charge,
                          double detachedDistance, double maxFraction) {
  for (unsigned i=0; i<2; ++i)
    if (!std::isfinite(distances[i]) || distances[i]<0. ||
        !std::isfinite(charge[i]) || charge[i]<=0.) return -1;
  const double total = charge[0]+charge[1];
  if (!std::isfinite(total)) return -1;
  for (unsigned i=0; i<2; ++i)
    if (distances[i]>detachedDistance && distances[1-i]<=detachedDistance &&
        charge[i]/total<maxFraction) return 1-i;
  return -1;
}
// Exactly one distant, small fragment; every other shower must be near the vertex.
inline int detachedFragment(const std::vector<double>& distances,
                            const std::vector<double>& charge,
                            double detachedDistance, double maxFraction) {
  if (distances.size()!=charge.size() || (charge.size()!=2 && charge.size()!=3)) return -1;
  double total=0.;
  int fragment=-1;
  for (std::size_t i=0; i<charge.size(); ++i) {
    if (!std::isfinite(distances[i]) || distances[i]<0. ||
        !std::isfinite(charge[i]) || charge[i]<=0.) return -1;
    total+=charge[i];
    if (distances[i]>detachedDistance) {
      if (fragment>=0) return -1;
      fragment=static_cast<int>(i);
    }
  }
  if (!std::isfinite(total) || fragment<0 || charge[fragment]/total>=maxFraction) return -1;
  return fragment;
}
// Angle between the main shower direction and the vector to the fragment start.
// Invalid/zero vectors cannot establish an attachment.
inline double attachmentAngle(const std::array<double,3>& direction,
                              const std::array<double,3>& displacement) {
  double dot=0., a=0., b=0.;
  for (unsigned i=0; i<3; ++i) {
    if (!std::isfinite(direction[i]) || !std::isfinite(displacement[i])) return 180.;
    dot+=direction[i]*displacement[i]; a+=direction[i]*direction[i]; b+=displacement[i]*displacement[i];
  }
  if (!std::isfinite(a*b) || a<=0. || b<=0.) return 180.;
  const double cosine=std::max(-1.,std::min(1.,dot/std::sqrt(a*b)));
  return std::acos(cosine)*180./std::acos(-1.);
}
// Shared across attempts in an event, including both possible fragment parents.
struct RecoveryBudget {
  std::size_t trials=0, showers=0;
  bool canTry(std::size_t maxTrials, std::size_t maxShowers) const {
    return trials<maxTrials && showers<=maxShowers && maxShowers-showers>=2;
  }
  bool canReconstruct(std::size_t count, std::size_t maxShowers) const {
    return showers<=maxShowers && count<=maxShowers-showers;
  }
};
} }
#endif
