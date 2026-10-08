# Published-genome pairs with the v0.5 engine

This workflow accepts a matched FASTA, repeat annotation and gene GFF3 for
each of two published assemblies, plus an existing minimap2 PAF in the same
query/target order. It runs the v0.5 `pair` command; legacy unique-TE output
tables are never converted into v0.5 biological states.

Requirements: the built TEvoX executable, Python 3.8+ with NumPy, and gzip.
The C core and its five public Python helpers do not require NumPy. Large
analyses need substantial RAM and disk space; see the measured figures in
[the integration report](../../docs/INTEGRATION_VALIDATION.zh-CN.md).

## Reproduce the source inputs

The [source preparation snapshot](../validation/results/source-preparation-snapshot.tar.gz)
preserves the earlier download catalogs and exact selection/conversion code,
including HG002 bigRmsk fragment extraction and SHRSP duplicate-row removal.
Unpack it outside this checkout and run its `initialize_study.py` to create a
new study workspace. Its README documents species-specific downloads and
the pinned minimap2 build; initialization itself downloads nothing.

Use that workspace's `scripts/run_pair.py` only for `--phase prepare` and
`--phase align` when reproducing the original alignments. Its historical
`finish`/`all` commands target the earlier engine. Pass the resulting directory
to the **current** v0.5 workflow:

```sh
python3 studies/published-genomes-v05/run_pair.py \
  --legacy-result /data/public-genomes/results/PAIR \
  --tevox ./tevox --output /results/v05/PAIR
```

For the Arabidopsis pair, add `--family-policy classification`. The snapshot
is source/provenance only and contains no downloaded genome or completed
inference outputs. Local file paths must be initialized before download.

## Run a pair

```sh
python3 studies/published-genomes-v05/run_pair.py \
  --manifest /data/pair.json --paf /data/A_B.paf.gz \
  --tevox ./tevox --output /results/pair
```

Minimal manifest (paths resolve relative to the JSON file):

```json
{
  "genomes": [
    {"id": "A", "fasta": "A.fa.gz", "repeats": "A.repeatMasker.out.gz",
     "genes": "A.genes.gff3.gz", "repeat_format": "repeatmasker-out"},
    {"id": "B", "fasta": "B.fa.gz", "repeats": "B.repeatMasker.out.gz",
     "genes": "B.genes.gff3.gz", "repeat_format": "repeatmasker-out"}
  ]
}
```

`normalize_args` may supply explicit provider coordinate/column definitions.
`normalize_annotations.py --help` documents RepeatMasker, GFF and BED inputs.
Never infer coordinate convention from a filename. Default `repeat-name`
families retain repeat-library names. For sample-specific EDTA family IDs,
use `--family-policy classification` only when their class/superfamily
labels have a common meaning; this provides coarser compatibility, not
evidence that sample-specific family identifiers are homologous.

For the preserved initial study directory the equivalent command is:

```sh
python3 studies/published-genomes-v05/run_pair.py \
  --legacy-result /data/frozen-study/results/human \
  --tevox ./tevox --output /results/human
```

The directory must contain `manifest.json` and `alignment.paf.gz`. It is a
source manifest/PAF container, not a requirement to run an old executable.
The manifest must refer to the exact selected source files. Raw downloads,
chromosome selection, repeat-fragment extraction and duplicate removal are
separate provenance steps; changing an assembly invalidates a matching
annotation/alignment bundle.

## Computation and interpretation

1. Freeze source input, PAF, executable and script SHA256 values. Reject a
   partial run if a source or configuration changes. Verify prepared FASTA
   cache bytes before reusing them. Completed runs require a new output path.
2. Normalize repeat records to BED8 while preserving fragment IDs, source
   metadata and exclusion reasons. SVA is retained. Gene GFF3 is used for
   post-analysis gene-body overlap and distance, not as an MCScanX prior.
3. Reconstruct exact `=`/`X` operations from the **unchanged alignment path**
   and both FASTAs. Ambiguous bases count as mismatches. This supplies local
   identity; it is not a new independent alignment or reciprocal validation.
4. Execute `tevox pair --gzip-output`. Keep all 17 output tables, native run
   metadata, streaming row-count/state checks, resource measurements and
   source-fragment gene neighborhoods. Candidate limits retain the engine's
   defaults; these runs are not an untruncated model-training export.
5. Commit `complete.json` only after the pipeline checks succeed. Native
   `.run.json` alone signifies core success, not completed gene postprocessing.

The 17 native TEvoX tables and normalized BED8 use **0-based half-open**
coordinates. The supplementary `empty_site_source_fragments.tsv` and
`empty_site_genes.tsv` explicitly use `start1`/`end1`, which are **1-based
inclusive**; the gene-context QC JSON records this convention. Do not pass
those start1/end1 columns directly to a BED consumer without conversion.

An `EMPTY_SITE_CONFIRMED` decision is an evidence-gated absence at the
counterpart locus. Counts refer to **source annotation fragments**, which
may share one structural event. They do not by themselves establish TE
insertion direction, evolutionary mechanism, population frequency, gene
regulation, precision or calibrated posterior probability. Compare samples
only after considering repeat-library, gene-model and alignment differences.

## Archive and verify

```sh
python3 studies/published-genomes-v05/archive_run.py /results/pair
python3 studies/validation/verify_completed.py /results/pair \
  --output /results/integrity.json
```

Archiving verifies FASTA restoration against the selected compressed source,
round-trips newly compressed plain tables, and writes a recovery journal
before deleting redundant plain copies. Native gzip tables are already
archived. `archive.json` records how to recover original paths/bytes; the
original `run.json` remains immutable. Verification hashes current artifacts,
source inputs and archive sources; `--roundtrip` additionally recomputes
uncompressed hashes. Copy an archived run together with its source/manifest
bundle; moving files without rebasing or preserving paths is not supported
by that historical run's absolute provenance.

The five-species executed manifests and source/input hashes are captured in
`../validation/results/`. Their recorded paths describe the validation
workspace. They are evidence snapshots, not an instruction to download a
different “latest” assembly or silently substitute missing annotations.

## Workflow tests

```sh
python3 studies/published-genomes-v05/test_inputs.py
python3 studies/published-genomes-v05/test_run_integrity.py
python3 studies/validation/test_truth_mapping.py
```

The integrity tests run a tiny complete pair from relative input paths and
confirm that changed repeats, changed PAF and tampered FASTA caches prevent
resumption. External truth and controlled simulation instructions are in
`../validation/README.md`.
