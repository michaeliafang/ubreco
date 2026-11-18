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
