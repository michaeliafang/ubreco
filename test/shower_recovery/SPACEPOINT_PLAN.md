# Reuse NuGraph spacepoints for shower recovery

The original reco2 dump contains fitted `sps::OverlayRecoStage2n` spacepoints,
PointCharge values, and SpacePoint–Hit associations, plus pre/noreg instances.
The v10_04_07_20 standard configuration schedules nuslhits -> sps -> NuGraph ->
gaushit. It uses SpacePointSolver on nuslhits (0.3 cm intersection tolerances,
MaxNTriplets=100000). This was traced through:
https://github.com/uboone/uboonecode/blob/v10_04_07_20/fcl/reco/MCC10/reco_uboone_data_mcc9_10.fcl

The redundant standalone solver job was removed in a normal follow-up commit;
shared branch history is not rewritten. No new solver is run by this workflow.

## Hit mapping

ExistingSpacePointHitMapper reads the original sps and nuslhits with explicit
process tags and the current nugraphshowerhits. Both intervening selectors copy
recob::Hit unchanged. It builds indexed fingerprints including channel, complete
WireID, view, signal type, start/end ticks, peak time and uncertainty, RMS,
amplitude and uncertainty, integral and uncertainty, summed ADC, multiplicity,
local index, goodness of fit, and degrees of freedom. Exact equality is required.

Complexity is O((N+M) log(N+M)), plus a linear traversal of spacepoint associations;
there is no all-pairs hit scan. Only unique source-to-selected matches are used.
Duplicate fingerprints on either side and nonfinite data are left unresolved.
Source hits not selected as shower hits are normal and do not constitute errors.
Missing selected hits are reported separately because they may indicate changed
hit values, a different production configuration, or missing input coverage.
A mismatched source ProductID in the input associations is an explicit error.

Products under `spsshowerhits`:
- default instance: Assns<SpacePoint,Hit> linking ORIGINAL sps points to CURRENT
  nugraphshowerhits. These retain all many-to-many relationships for mapped hits.
- `source`: vector<art::Ptr<Hit>>, parallel to the selected hit collection. Entry
  j points to the original nuslhits hit for selected hit j; null means unresolved.
  This uses an existing RecoBase dictionary, avoiding a custom Hit–Hit dictionary.

No new hits, spacepoints, point charges, or ntuple are produced. A point may have
only one selected-hit association; the mapper preserves it but recovery requires
at least two planes for that point to contribute to its consistency check.

## GPVM run

Build/install the changed ubreco modules and FCL, then run on the original input
containing nuslhits and sps (do not first run the removed preprocessing job):

```bash
lar -c run_fullshowerreco_overlay_numi_reprocess_sps.fcl -s input.root -n 10 -o shower_existing_sps.root
```

The shower label remains `showersplittest` for gallery testing. The new FCL
schedules nugraphshowerhits, then spsshowerhits, then clustering and recovery.
Inspect logs for ExistingSpacePointHitMapper's matched/missing/ambiguous counts
and mapping_ms. Original unselected nuslhits need not be present in the selected
collection. All mapped association products are retained by the existing output.

Recovery uses these links as a cross-plane consistency veto on proposed matched
showers, not yet as a 3D trajectory generator. See Recovery/README.md for weights,
thresholds, fallback, and computation budgets. With insufficient mapping/evidence
it falls back to the 2D search; missing required input products are job errors.
To compare without the consistency veto, append to a local wrapper:

```fcl
physics.producers.showersplittest.SingleShowerRecovery.UseSpacePoints: false
```

The mapping still runs in that comparison so differences isolate its use in
recovery. A separate timing run with the ordinary splitlabel FCL measures the
mapping overhead itself. Verify persisted links with eventdump and, when reading
associations, use sps as the SpacePoint collection and spsshowerhits as the
association producer. Gallery need not expose an association-only producer.

Standalone tests cover reordered subsets, normal source omissions, duplicate
ambiguities, changed fingerprints, nonfinite values, cross-plane agreement,
conflicting assignments, one-plane fallback, and shared-hit weighting. Full art
compilation and event behavior require GPVM validation; they cannot be verified
in this local environment. Existing selected-hit copies from a different chain
must not be assumed identical without checking the mapper's logs.

Local benchmark: mapping 30,000 source hits to 10,000 selected hits took about
21 ms (one standalone optimized run, excluding art I/O and association traversal).
