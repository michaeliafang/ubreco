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
