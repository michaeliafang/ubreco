# Standalone splitter checks

These run the actual `Recovery/VertexSplit.h` implementation, without art,
ROOT, or LArSoft. They do not validate event associations or full reconstruction.
Use the separate `dt-shr-gpt` Conda environment.

From the repository root:

```bash
conda run -n dt-shr-gpt x86_64-conda-linux-gnu-c++ -std=c++17 -O2 -Wall -Wextra -Werror -pedantic -I . test/shower_recovery/vertex_split_test.cc -o /tmp/vertex_split_test
conda run -n dt-shr-gpt /tmp/vertex_split_test
conda run -n dt-shr-gpt x86_64-conda-linux-gnu-c++ -std=c++17 -O2 -Wall -Wextra -Werror -pedantic -I . test/shower_recovery/vertex_split_driver.cc -o /tmp/vertex_split_driver
MPLCONFIGDIR=/tmp/dt-shr-mpl XDG_CACHE_HOME=/tmp/dt-shr-cache conda run -n dt-shr-gpt python test/shower_recovery/synthetic_patterns.py --driver /tmp/vertex_split_driver --output /tmp/dt-shr-synthetic --seeds 20
```

The assertion tests check basic geometric cases and partition integrity.
The synthetic-pattern study is exploratory: undesirable splits are **reported**,
not hidden by an assertion or tuned away. Nine geometries each have 20 fixed-seed
variations. Shapes include widening clouds around a continuous backbone,
lognormal charge variations, 0.3 cm wire quantization, missing-hit bands, bridges,
crossing branches, and an offset vertex. Noise hits, detector response, induction
signal shapes, hit RMS, electronics, and 3D reconstruction are not simulated.
All generated hits are supplied as one candidate cluster; the study separately
reports connected components of the splitter's Euclidean scoring graph. These
are not necessarily components of the original proximity compatibility graph.

Known synthetic shower identity is used only for evaluation, never as splitter
input. Branch hits retain their parent's identity. For two-shower cases, the two
output labels are permuted to maximize correctly assigned charge. Hit accuracy,
charge accuracy, branch accuracy, and minimum per-shower hit recall use that same
permutation. Accuracy summaries are conditional on an accepted split; acceptance
counts are reported separately. A split of a one-shower control is a false split.
Examples always use seed 1000, rather than selecting favorable examples.

The study writes per-trial CSV, JSON summaries, labeled example CSVs, and PNG
truth/assignment comparisons. Results are toy stress measurements, not physical
signal efficiencies or neutrino/EXT background rates. Splitter defaults remain
unchanged by the study.

Input hit ownership: recovery removes repeated hit pointers within each cluster
before computing charge or split geometry, preserving first-occurrence order.
Overlapping input clusters remain independent candidates. A final usable shower
pair is rejected if its members share any hit pointer; unused clusters do not
trigger this guard. At least one selected split must supply its actual daughter
clusters to different output showers. Hit identity includes product ID and key.
Diagnostics report `duplicates_removed`, `overlapping_input_hits` (extra cluster
references after within-cluster deduplication), and `ownership_rejected`.
`hit_ownership_test.cc` checks the standalone deduplication and disjointness
helpers. Full matching, reconstruction and event associations require GPVM tests.

With `SingleShowerRecovery.SaveSplitCandidates: true` (enabled in the recovery
FCL), the reconstruction module also writes clusters and cluster-to-hit
associations with instance name `showersplitCandidates`. In the split-label job,
the full input tag is `showersplittest:showersplitCandidates`. This is diagnostic
output, not another shower producer. Consecutive clusters (0/1, 2/3, etc.) are
the two daughters of each retained candidate, ordered by plane then candidate
rank. The log maps each pair to its original cmerger cluster. Candidates are
saved before cross-plane matching, including when the baseline shower is retained.
These are the best split per source cluster, subject to MaxCandidatesPerPlane,
not every angular trial. Empty output means no retained candidates or recovery
was not attempted. Alternatives may share hits; they are not simultaneous showers.
Disable SaveSplitCandidates to omit these extra products. Existing shower and
accepted recovery products are unaffected. Validate the art products/associations
and gallery display on GPVM; standalone tests cannot exercise their publication.

Per reconstructed recovery candidate, logs now show its planes, reconstruction
pass flag, squared start/direction magnitudes and final usability. Detailed tool
exceptions remain in the full log (they may be excluded by the diagnostic grep).

Start/valley update: defaults now allow 3 cm support gaps, disable the early-hit
fraction requirement (MinTrackFraction: 0), and use a 5 degree minimum opening.
The count/span/near-vertex support requirements and 20% daughter charge cut remain.
The gap relaxation is not dead-wire-aware; detached fragments can still fail.

Among candidates passing all existing cuts, ranking now prefers a larger valley
score, then a smaller graph cut fraction. It no longer prefers larger opening
angles on ties. Four equal longitudinal bands cover 2--10 cm by default along
the fitted direction bisector. Each compares equal-area strips centered on the
two rays and their midpoint. The common half-width is min(TrackWidth, separation/4);
regions narrower than 0.15 cm contribute no evidence. A band needs at least two
hits in each side strip; its score is max(0, 1 - middle_count/min(side_counts)).
The overall score averages all four bands, including zero for unsupported bands.
Thus an empty side cannot create positive evidence, and one crossing branch
need not erase evidence elsewhere. The score uses hit counts, not charge.
`valley_score` and `valley_bands` are logged for each evaluated large cluster's
best split. No valley threshold rejects candidates. This is not a dead-wire or
induction-artifact discriminator; those effects can still bias the score.
Tests cover missing-hit gaps, off-axis dilution, empty/single-sided bands,
filled valleys, a localized bridge, rotation and ranking fallback.

The split-label and existing-SPS jobs now schedule a display-only copier under
module label `showersplitCandidates`. Its clusters and cluster-to-hit associations
both use the default instance, avoiding viewers' loss of instance labels. The
source `showersplittest:showersplitCandidates` remains intact. The copier logs
cluster and association counts under `[ShowerSplitCandidateDisplay]`; it does not
rerun splitting or alter reconstruction. Full art publication needs GPVM validation.

Investigation diagnostics: DiagnoseTwoShowers (enabled in recovery FCL) generates
split candidates for exactly two usable baseline showers, but returns before
cross-plane matching or reconstruction and never replaces their output. Thus 614
can be inspected without enabling two-shower recovery. No distance/charge-based
two-shower eligibility cuts are applied in this diagnostic mode.

Spacepoint diagnostics report separately normalized plane-pair conflict/weight
for 01, 02, 12; a conservative subset with one assigned hit per represented plane
and no assigned-hit reuse across points; and four hypothetical label permutations
when there are exactly two matched candidates. Flip masks 0/1/2/3 mean unchanged,
flip plane 1, flip plane 2, flip both (plane 0 fixed). These hypotheses do not
rematch clusters or establish valid showers. The uniqueness subset is based on
assigned mapped hits only, not all original solver associations, and is not proof
of physical correspondence. Zero conflict with zero weight means no evidence.
None of these diagnostics changes the existing aggregate conflict veto. Their
additional work is linear in assigned associations per fixed diagnostic pass.

## Initial independently switchable multi-shower recovery (historical)

`RecoverTwoShowers` enables replacement only for exactly two usable baseline
showers with disjoint hit lists: one starts >20 cm from the reconstructed vertex
and has <30% of their summed positive collection-plane hit integral; the other
starts <=20 cm away. Hits are deduplicated; missing/invalid collection charge
fails eligibility. Using only collection charge avoids summing unequal view
coverage. These distances and fractions are configurable. No inference is made
from cone length or display angle.

In the initial implementation (superseded by the combine/merge update below), the detached shower's hits were excluded from trial inputs. Coverage is measured
against the remaining baseline shower. A successful replacement must have two
usable, disjoint showers, both with 3D starts <=5 cm from the vertex, and must use
both daughters of an actual split. Failed recovery preserves every baseline
product. Source cmerger products are not modified. This is an experimental
selection: excluding detached hits may discard genuine branch charge.

`DiagnoseThreeShowers` adds candidate-only output for exactly three baseline
showers (e.g. 4867), without matching or changing shower output. Ineligible
baseline pairs still use `DiagnoseTwoShowers`. Near-vertex singleton/fragment
handling remains unchanged and deferred. New options default off in C++ but
are enabled in the experimental recovery FCL; set each false for independent
comparisons. Local tests cover eligibility, not the art integration.

## Alternative cross-plane pairings

`TryAlternativePairings` enumerates the distinct swaps for exactly two matched
showers with one cluster in each common plane (four possibilities for three
planes; two for two planes). Swaps are actual cluster reassignment before
reconstruction, not merely relabeled diagnostic points. They must pass the
configured matching tool, hit exclusivity, and the existing spacepoint conflict
cut with sufficient weight. Surviving assignments are ranked by conflict and
reconstructed within the shared MaxRecoShowers budget. No third match is dropped
to force a two-shower result. With UseSpacePoints false, alternatives are disabled.
`pairing_swap_mask` uses plane bits 1/2/4 (unlike the old diagnostic flip mask,
which used 1/2 to mean planes 1/2); `reconstruct_swap_mask` records the choice.
The original pairing remains the sparse-evidence fallback.

## Supported two-plane fallback

This describes the original exactly-two-group fallback. The independent
substantial-cluster search below extends this to other matcher group counts.

`AllowTwoPlaneFallback` tries collection+plane0 or collection+plane1 subsets of
exactly two matched candidates; it never drops an entire third matched shower.
Each retained plane must have exactly one cluster per shower, >=MinHits hits per
cluster, and each shower must carry >=MinChargeFraction of that plane's pair
charge. The configured matching tool must accept each pair. Spacepoint conflict
must pass with >=MinSpacePointWeight overall, and EACH shower must have
>=MinTwoPlaneShowerWeight consistent point weight (default 3). Multiplicity
weights are computed before separating showers. Missing evidence cannot count
as agreement. Induction-only reconstruction is not enabled.

Within each splitting trial, full-view hypotheses precede two-plane fallbacks;
within each group, lower conflict ranks first. All trials/options share the
existing reconstruction budget, so a fallback is not guaranteed to be reached
if earlier reconstructions exhaust it. Coverage uses only original hits in the
retained planes for fallback; logs report considered, retained and omitted hit
counts. The omitted plane is not assigned to the recovered showers. Source
clusters and diagnostic products remain available. Pairing alternatives and
fallback can each be disabled independently; with UseSpacePoints false both
are inactive. A low conflict score still is not proof of correct reconstruction.

## November 20 change ledger and rollback

- `fc78169`: gated actual two-shower recovery; candidate-only three-shower tests.
- `e475a65`: spacepoint-ranked alternative cross-plane cluster assignments.
- `4cd5c98`: supported collection-plus-induction fallback and retained-view coverage.
- Reconstruction failure logging: `[ShowerRecoFailure]` identifies candidate,
  failing tool and exception reason. Correlate with surrounding event/trial logs.

For independent comparisons, override these under
`physics.producers.showersplittest.SingleShowerRecovery`:
`RecoverTwoShowers`, `TryAlternativePairings`, `AllowTwoPlaneFallback`, and
`DiagnoseThreeShowers`. For the combine/merge update below, also disable
`RecoverThreeShowers` and `SearchSubstantialClusterPairs`; disabling
`RecoverTwoShowers` also disables combine-then-split. Reverting
commits is an alternative, in reverse dependency order; do not rewrite history.

Suggested GPVM checks: 616 preserves its successful recovery; 650 exercises
pairing changes; 602 exercises retained planes 1+2; 614/4881/4890 exercise eligible
two-shower recovery if they satisfy the explicit charge/distance cuts; 4867 gets
candidates (and now tests fragment attachment when `RecoverThreeShowers` is enabled); 4901 reports the failed
reconstruction tool. Event displays alone do not establish eligibility fractions.
785 still cannot borrow singleton start hits: fragment handling, induction-streak
suppression and induction-only reconstruction remain deferred. No cmerger edits.

Local validation used all four standalone tests with C++17, -O2, -DNDEBUG,
-Wall -Wextra -Werror -pedantic. These validate helper logic; they do not compile
art/LArSoft integration or demonstrate physics performance. No claim about full
job runtime or successful reconstruction is possible until GPVM testing.


## Combine/split, three-shower attachment, and substantial-cluster search

The current experimental FCL enables three independently switchable additions:

- `CombineTwoShowerHits`: for eligible two-shower events, combine **both** baseline
  showers' associated hits into one deduplicated cluster per plane, then apply the
  existing splitter. Eligibility remains a detached start strictly beyond 20 cm
  with strictly less than 30% of total positive collection charge, and the other
  shower within 20 cm. No detached hits are excluded. `MinOriginalHitFraction`
  (0.8) is checked for each baseline shower separately as well as their union.
  For an accepted two-plane fallback, coverage uses only retained planes and
  omitted hits are logged. Switching this off uses ordinary cmerger inputs with
  aggregate coverage of both showers; switch `RecoverTwoShowers` off to disable
  two-shower replacement entirely.
- `RecoverThreeShowers`: exactly one detached shower must be beyond 20 cm and
  below `ThreeShowerMaxFragmentChargeFraction` (0.20), with both others within
  20 cm. Try attaching all fragment hits to either main shower, grouped by plane.
  Rank parents by the angle from the main shower direction to the fragment start;
  require at most `MaxFragmentAttachmentAngle` (35 degrees). This is a tunable
  geometric criterion, not a validated physics threshold. Keep main assignments
  fixed: no splitting, swapping main identities, or dropping a plane in this path.
  Require the configured matching tool, mapped-spacepoint support for both
  showers, the conflict cut, per-plane charge balance, per-original-shower hit
  coverage, and two usable reconstructed showers starting within 5 cm. If neither
  attachment passes, preserve the original three showers. `DiagnoseThreeShowers`
  still writes split diagnostics from the original cmerger inputs independently;
  those diagnostic clusters are not the attachment hypotheses.
- `SearchSubstantialClusterPairs`: permit collection-plus-induction matching
  independent of how many groups the original matcher returned. Consider the top
  `MaxPairClustersPerPlane` (3) clusters by hit count, with charge as a tie-breaker.
  Choose two per plane and test both assignments: at most 36 lightweight options
  per split trial. Each pair must retain `MinPairPlaneChargeCoverage` (0.80) of
  positive charge in **each retained plane**, measured against all clusters with
  at least `MinHits`, even those outside the top-three cap. Overlapping input
  clusters can make this denominator conservative. Require both split daughters
  to contribute, disjoint hit ownership, existing matcher/spacepoint/support/
  charge-balance checks, and original-hit coverage. Existing matching hypotheses
  are tried first; extra options prioritize selected hit count, then conflict.
  Logs show original group count and unused retained-plane cluster/hit counts.

The common event budget stays at four matching trials and four additional
reconstructed proto-showers, shared by both possible fragment-parent attempts.
This bounds reconstruction work, not wall time; measure timing on GPVM. No change
was made to singleton near-vertex evidence, dead-wire handling, induction-only
reconstruction, or `FilterShowers` angle cuts. Product labels remain unchanged.

All three additions default off in C++ and are enabled in the experimental FCL.
For isolated comparisons, toggle `CombineTwoShowerHits`, `RecoverThreeShowers`,
and `SearchSubstantialClusterPairs` separately. Standalone tests cover strict
eligibility boundaries, invalid inputs, attachment angles, shared budgets,
hit conservation/deduplication, both pairing assignments with extra groups,
charge coverage beyond the search cap, and the bounded option count. They do not
compile the art module or validate real event physics.

GPVM checks: run 5512/97/4867 for three-shower attachment and 5512/96/4827 for
independent plane-pair selection; compare 12/602, 12/614, 12/616, 97/4881, and
97/4890 against the visually successful previous output. Inspect `recovery_mode`,
`fragment_parent`, `attachment_angle`, `initial_match_groups`,
`substantial_pair_search`, `coverage_group`, `event_reco_showers`, and `status`.
The source builds and event outcomes still require GPVM validation.
