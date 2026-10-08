# Frozen integration validation evidence

Read [the Chinese validation report](../../../docs/INTEGRATION_VALIDATION.zh-CN.md)
for conclusions and limits. This directory contains compact evidence, not the
large genomes or all inference tables. The complete runs remain under
`/workspace/tevox-v05-study/final/{human,mouse,rat_corrected,arabidopsis,maize}`.
Old runs and source genomes are preserved separately.

- `cohort-summary.json`: five final runs, fragment/locus counts, directional
  empty-site fragment counts and measured child resource usage.
- `cohort/*/`: immutable input/signature, completion and lossless archive
  records; preparation/gene QC; verified EQX cache reuse records where used.
- `final-cohort-integrity.json`: current bytes, original selected inputs and
  archive sources checked against recorded SHA-256. Archive creation checked
  decompressed round trips. The post-check hashes the bound compressed source;
  it does not claim that every large source was decompressed a second time.
- `*-evidence-contract.json`: every evidence row checked using bounded memory.
  `arabidopsis-final-schema.json` records the complete relational schema check
  on the final Arabidopsis run, with archived inputs restored and rehashed.
- `matrix-summary.json`, `platforms/`, `matrix-source-sha256.json`,
  `core-source-match.json`: nine real Linux userspace-container checks and
  source bindings. The host kernel is shared.
- `verified-ci.json`, `ci-verified-*`, `macos-verified-*`: successful native
  GitHub CI runs, job metadata and uploaded Intel/ARM package artifact metadata.
  The Mac checks ran on macOS 15; minimum deployment target is not an OS test.
- `sources/index.json`: maps original source/selection manifest paths to
  content-addressed copies in `sources/`; URLs and exact assembly versions
  are retained. Rewrite local paths to downloaded, checksum-matching files
  when reproducing elsewhere.
- `source-preparation-snapshot.tar.gz` and `.json`: byte-bound, source-only
  snapshot of the earlier download/selection/conversion utilities and source
  catalog (no genomes or completed inference outputs). This preserves the
  exact HG002 bigRmsk extraction, rat deduplication and plant preparation code
  referenced by the provenance. Its historical analysis wrappers target the
  earlier engine; use `studies/published-genomes-v05/run_pair.py` for final
  v0.5 inference. Source-preparation tests are recorded separately.
- `truth/giab/`, `giab-evaluation.json`: 42 externally selected GIAB absence
  labels and evaluation against final human decisions. These are positive-only.
- `truth/mouse-reference/`, `mouse-pcr-truth-final.json`,
  `mouse-pcr-evaluation.json`: predeclared Table S2/S3 PCR selection, frozen
  reference allele responses, unique transfer and final mouse state evaluation.
  B6 PCR is C57BL/6J, not C57BL/6NJ. CAST evidence is matched by strain;
  identity of the sequenced and PCR animals is not established.
- `calibration-readiness.json`: explicit incomplete calibration acceptance
  criteria; no fitted probability model is published.
- `polyploid/summary.json`, `scaling/summary.json`,
  `synthetic-membership-audit.json`: all 36 synthetic copy-context runs,
  four sparse scale runs and independent membership-universe rechecking.
  Missing-annotation B-cubed scores evaluate observed members, not recovery
  of unannotated TE nodes. Seeds preserve geometry; these are not biological
  replicates. See `validation-summary.png` / `.pdf` for the corresponding plot.
- `evidence-sha256.json`: final byte hashes of the compact evidence files.

`run-matrix.*`, `prepare-container*`, `reuse-verified-alignments.py` and
`install-cloud.sh` preserve the actual cloud commands and workspace layout.
They are evidence of this run, not general-purpose distribution installers.
The reusable public build instructions are in `docs/PLATFORMS.md`.

The five pair analyses use DNA alignment evidence. Downloaded gene annotations
are used for descriptive candidate context; they do not supply a native
MCScanX prior in these five runs. Maize gene models are AUGUSTUS predictions.
Final counts refer to annotation fragments and reconstructed loci, not proven
independent insertion events or population frequencies.
