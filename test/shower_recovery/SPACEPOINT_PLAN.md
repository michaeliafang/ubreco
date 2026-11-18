# SpacePointSolver feasibility and preprocessing

Before this change, ubreco had no SpacePointSolver job or consumer.
The new preprocessing job is a wrapper; it does not supply the solver itself.
The upstream larreco interface accepts HitLabel, builds spacepoints, and with
Fit enabled writes SpacePoint–Hit associations and charged spacepoint products.
Its standard FCL is SpacePointSolver.fcl / standard_spacepointsolver. These
names were checked against published upstream source, not the installed GPVM
release. The experiment-specific HitReaderTool and services need validation in
that release before production use.

Sources:
- https://sbnsoftware.github.io/doxygen/df/d9e/SpacePointSolver__module_8cc_source.html
- https://sbnsoftware.github.io/doxygen/db/df2/SpacePointSolver_8fcl_source.html

The experimental job run_spacepointsolver_shower_inputs_overlay_numi.fcl reuses
this repo's services and hit selection and schedules only nugraphshowerhits,
neutrinovertex, and spacepointsolver. It does not run the splitter or rebuild
Pandora. Use the same original input as the ordinary reprocessing job.

After installing the FCL on GPVM, run in a separate output directory:

```bash
lar -c run_spacepointsolver_shower_inputs_overlay_numi.fcl -s input.root -n 10 -o shower_spacepoints.root
```

Inspect spacepointsolver (default instance, fitted result), including hit
associations and density near the vertex in the known merged events. The
pre/noreg instances are different diagnostic outputs, not interchangeable with
the fitted result's hit associations. Record TimeTracker timing and compare to
the normal reconstruction job on the same events. Solver runtime is not yet
known; computing it once lets later splitter tests reuse the result.

NuGraphShowerHitsByPfp copies hits into a new product without a hit-to-original-hit
mapping. Therefore solver hits from gaushit cannot simply be compared by art::Ptr
to cmerger hits from nugraphshowerhits. The preprocessing job solves directly on
nugraphshowerhits and retains all products. A future second-stage job must REMOVE
nugraphshowerhits and neutrinovertex from its schedule, reuse those stored
products, and retain spacepointsolver and its associations. Recreating selected
hits in a later process would break direct product identity even if values match.

No existing recovery code reads spacepointsolver yet. Running the normal recovery
job on the preprocessing file would not add 3D information to the splitter.
Implement an explicit optional consumer only after validating the point cloud.
Use its SpacePoint–Hit links to restrict seeds to the baseline shower and connect
consistent early trajectories across planes. Ambiguous hits can support multiple
candidate points; avoid treating every point as independent evidence or forcing
all hits to have a unique 3D point. Preserve 2D fallback and unmatched hits.

Potential benefit is resolving branches that cross only in one projection and
coordinating the two starts across planes. Improvement is a hypothesis, not an
established result. Sparse/ambiguous points could erase precisely the early
structure sought. The selected-hit input also cannot recover hits excluded by
the upstream hit selector; a full-hit solver plus a reliable mapping is a later
option if needed. Existing ShrReco3D SpacePointAlg runs after reconstruction on
already assigned shower hits; it is not this independent precomputed evidence.

## Configuration provenance

`SpacePointSolver.fcl` and its `standard_spacepointsolver` table come from the
external LArSoft `larreco` package, under `larreco/SpacePointSolver/` in upstream
source. The standard FCL also includes `spacepointsolverhitreadingtools.fcl`.
Neither external FCL was copied into ubreco. The new ubreco job includes the
installed `SpacePointSolver.fcl` via FHICL_FILE_PATH and overrides only HitLabel
and Fit for this study. The module output label is `spacepointsolver`; these
points are a separate reusable product, not split-shower output.

The interface was inspected through SBN's generated source documentation linked
above. This is evidence for the upstream interface, not proof that the exact
GPVM installation has that file or identical tool defaults. Verify the installed
larreco configuration and experiment-specific hit-reader settings before running.
