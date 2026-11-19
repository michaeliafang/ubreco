# Experimental single-shower recovery

## Status (2026-11-17)

The C++ splitter has passed standalone tests, including branching two-shower
patterns. The art integration has been reviewed against this repository's
ClusterMaker, matcher, ProtoShowerCMTool, and shower publication code.
It has **not been compiled or run inside LArSoft**: the local environment has no
available art/LArSoft installation or `lar` executable. Framework compatibility,
FCL resolution, event associations, and full reconstruction still require a GPVM
build and event run. Standalone success is not framework validation.

## Integration

This is an optional retry inside the existing `ShrReco3D` producer, not a separate
art module. `VertexSplit.h` is the detector-independent algorithm;
`SingleShowerRecovery.h` adapts hits and clusters to the existing matcher and
shower reconstruction manager. The CMake files install the headers and FCL and
link the producer against the cluster-matching library.

Recovery is disabled unless `SingleShowerRecovery.Enable` is true. It runs only
after the normal pass produces exactly one usable shower and requires one input
vertex. It considers clusters in all three planes, including unmatched clusters,
and retries the existing cross-plane matching and shower reconstruction. It
accepts exactly two usable showers, requires the configured original-hit retention,
and requires a selected split to contribute to both showers. Otherwise the normal
result is retained. Events with zero or at least two usable showers are untouched.
The existing collection-plane requirement remains in force.

Accepted recovery clusters, PFParticles, and their hit/cluster associations use
the `showerreco3d:recovery` instance. Reconstructed showers retain the usual
`showerreco3d` label. The accompanying `NuSliceReBuilderProducer` change handles
the recovered cluster-hit product by ProductID; it must be built and deployed too.
Recovery products are empty when recovery is enabled but not accepted.

## GPVM check

Copy all modified and new files, including the untracked Recovery directory,
test directory, and FCL files. Build and set up this checkout with the normal
MicroBooNE development workflow for this release; a standalone Conda compiler
cannot build the art integration. Both `ShrReco3D` and
`NuSliceReBuilderProducer` must compile successfully before running events.

The separate job FCL is
`run_fullshowerreco_overlay_numi_reprocess_recovery.fcl`. It includes the existing
NuMI overlay reprocessing job, enables the recovery table, and copies the
baseline match and priority tool settings. Its process name is
`ShowerReco3DRecovery`. `single_shower_recovery.fcl` holds the algorithm settings;
it is not a complete job on its own.

After building and setting up the modified release, run on the same original
input used for the baseline reprocessing job (replace `input.root`):

```bash
lar -c run_fullshowerreco_overlay_numi_reprocess_recovery.fcl -s input.root -n 10 -o recovery.root
```

Use a separate working directory for the baseline comparison to keep histogram
outputs separate:

```bash
lar -c run_fullshowerreco_overlay_numi_reprocess.fcl -s input.root -n 10 -o baseline.root
```

Do not feed the newly recovered output back into this same recovery process.
The input must contain the products required by the original reprocessing job.
Check that both jobs finish without exceptions, then compare identical event IDs:
events with two baseline showers must remain unchanged, unsuccessful recovery
must preserve the single shower, and successful recovery must publish two showers
with valid cluster-hit associations through the rebuilt Pandora output. The
`SingleShowerRecovery` log reports accepted recovery or baseline retention. A
small sample without a recovery acceptance does not validate the accepted path;
include at least one known merged event.

## Synthetic-test interpretation

The previous small study (before the tightened cuts below) uses common-vertex showers with multiple forward-going
branches. Twelve two-shower toys accepted splits; three single-shower controls
did not. Crossing-branch hit ownership and vertex errors remain limitations.
These counts do not establish signal efficiency or background purity.

For future synthetic work, make shower separation and spread increase farther
from the vertex, while retaining branching in both showers and occasional branches
pointing into the other shower. The user requested recording this correction
without another local synthetic run.

## Tightened recovery defaults (2026-11-18)

Each proposed daughter now needs at least 20% of its parent 2D cluster's
nonnegative integrated hit charge (previously 5%). This is not a hit-count
fraction or a fraction of all event charge. Vertex support now requires at least
three hits in the 5–10 cm annulus (previously 5–20 cm). The absolute minimum
remains ten hits per daughter. The 10 cm upper radius is an initial tuning choice;
these thresholds need event validation. No final 3D shower-size cut is added.

Release-mode standalone checks cover rejection below 20% charge, acceptance
above it, and rejection of balanced daughters beginning beyond 10 cm. The earlier
synthetic-study acceptance rates do not validate these new defaults.

## Gallery product-label diagnostic

Use `run_fullshowerreco_overlay_numi_reprocess_splitlabel.fcl` to run the same
recovery configuration under module label `showersplittest`, with process
`ShowerSplitLabelTest`. This is a separate job option; the regular recovery FCL
continues to produce `showerreco3d`. The diagnostic replaces the producer in the
scheduled path, and configures the downstream rebuilder's `ShowerProducer`
parameter to read both showers and their associations from `showersplittest`.
That parameter defaults to `showerreco3d` for all existing jobs.

After rebuilding/installing, run in a separate output directory:

```bash
lar -c run_fullshowerreco_overlay_numi_reprocess_splitlabel.fcl -s input.root -n 10 -o splitlabel.root
```

Select `showersplittest` as the shower product in gallery. Its appearance tests
product selection, not whether a split was accepted: it contains the final
showers even when recovery is skipped or rejected. Recovered clusters and
PFParticles use `showersplittest:recovery`. The rebuilt interaction also remains
available as `pandora`. Existing shower products retained from input may still
appear in the display. Return to the regular recovery FCL to use the standard
`showerreco3d` output label again. Local source checks do not replace a GPVM build
and FCL/event validation of this diagnostic.

## Early-trajectory revision (pending GPVM validation)

The new 2D candidate search includes hits inside 5 cm. Each proposed daughter
needs a sustained early track: at least six hits, at least three within 5 cm,
at least 3 cm of longitudinal span, and no longitudinal gap above 1.5 cm in
the supporting run. Inliers must lie within 0.6 cm of the fitted line and the
run must contain at least 60% of that daughter's hits in the 0–10 cm region.
These are projected 2D distances; they are not a measured 3D start radius.
Hits at exactly zero radius remain assigned but cannot determine a direction.

A bounded scan initializes two sides. Weighted line fits and one inlier-refit
per side identify directions with free offsets. Weights r²/(r²+1 cm²) reduce
leverage of the closest hits. Fitted rays assign all hits by nearest-ray distance;
track support is rechecked after assignment. This is a simple trajectory proxy,
not a full track reconstruction or robust multi-branch pattern recognizer.
The 10-degree opening, 20% charge, ten-hit daughter, and final two-shower
requirements remain. Contact edges with both endpoints within 5 cm no longer
penalize the cut, but these hits remain in fitting and ownership.

MaxEdges=200000 rejects dense scoring graphs. MaxMatchTrials=4 and
MaxRecoShowers=4 cap per-event retries and total extra proto-showers passed to
3D reconstruction. Entire match sets that exceed the remaining budget are
skipped, never truncated; this avoids hiding a third successful shower.
Multi-plane split combinations are tried first. These caps can miss a later
successful configuration, and do not guarantee a wall-clock slowdown factor.
The optional existing-spacepoint consistency consumer is described below.

Logs show baseline_usable and baseline_reco_ms, per-cluster track_pairs,
balanced_pairs, connection_pass and edge_limit, followed by trial/reconstruction
counts, coverage/contribution failures, budget skips, and recovery_ms. The
baseline time measures the reconstruction manager only, whereas recovery_ms
includes splitting and rematching. Use full-job TimeTracker comparisons on the
same input to establish total overhead.

Release-mode standalone checks pass for touching starts, uneven charge,
late-start and gapped-fragment rejection, jitter, vertex displacement,
partition completeness, invalid configurations, and graph limits. On one local
six-repeat median benchmark, old/new times in milliseconds were 0.368/0.338
(300 hits), 1.356/1.026 (1000), 7.129/4.403 (3000), and 12.428/10.852 (6000).
These test only the splitter, not art, matching, or full reconstruction.
The small branching study accepted 12/12 two-shower toys and 0/3 single controls;
mean charge ownership by case was 90.6%, 87.4%, 92.5%, 89.6%. This does not show
uniform ownership improvement over the earlier splitter; crossing branches
remain a limitation and real-event success has not been established.

Existing NuGraph spacepoint reuse is documented in
`test/shower_recovery/SPACEPOINT_PLAN.md`. The separate SpacePointSolver
preprocessing job has been removed: use the persisted `sps` collection instead.

## Existing spacepoint consumer

The dedicated `run_fullshowerreco_overlay_numi_reprocess_sps.fcl` enables the
mapper and recovery's `UseSpacePoints` option. Other recovery FCLs keep it off.
The consumer uses persisted spacepoint hit associations to check each proposed
cross-plane match before spending its extra 3D reconstruction budget. Points
with assigned hits in at least two planes provide evidence; a point is conflicting
when those hits belong to different proposed showers. Contributions are weighted
by the inverse of the largest participating hit's point multiplicity to reduce
repeated evidence from ambiguous points. Less than six effective points causes
2D fallback. Otherwise a conflicting fraction above 25% rejects that match set.
These are exploratory parameters, not physics-validated thresholds.

This first consumer checks consistency of hit assignments across views; it does
not fit 3D start trajectories or replace the 2D partition boundary. Positions and
charges remain available in the original sps products. Ghost points may still
bias the consistency test, so compare with UseSpacePoints=false on identical
events. Logs record point count, effective weight, conflicting fraction and
rejection decision, in addition to the geometric rejection and timing fields.
