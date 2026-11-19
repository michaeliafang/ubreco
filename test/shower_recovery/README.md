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
