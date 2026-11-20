#ifndef UBRECO_RECOVERY_ELIGIBILITY_H
#define UBRECO_RECOVERY_ELIGIBILITY_H
#include <array>
#include <cmath>
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
} }
#endif
