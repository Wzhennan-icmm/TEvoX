# External truth and controlled validation

All truth has a stated origin. A synthetic ancestry label, an independently
produced assembly callset and an experimental validation are different
evidence sources; combining them without these distinctions is unsupported.

## GIAB HG002 positive absence subset

`giab_empty_truth.py prepare` takes the official HG002 SV Tier1 v0.6 VCF,
v0.6.2 confident-region BED, UCSC `hg19ToHs1.over.chain.gz` and its reverse,
the exact CHM13 FASTA and the frozen normalized CHM13 TE BED8. Example:

```sh
python3 studies/validation/giab_empty_truth.py prepare \
  --vcf HG002_SVs_Tier1_v0.6.vcf.gz \
  --regions HG002_SVs_Tier1_noVDJorXorY_v0.6.2.bed \
  --forward hg19ToHs1.over.chain.gz --reverse hs1ToHg19.over.chain.gz \
  --fasta chm13.fa --te chm13.repeats.bed --output truth/giab
python3 studies/validation/giab_empty_truth.py evaluate \
  --truth truth/giab/truth.tsv --decisions results/tevox.decisions.tsv.gz \
  --output truth/giab/evaluation.json
```

Selection requires PASS homozygous alternate deletions >=50 bp, a simple
sequence-resolved allele, complete Tier1 containment with 50 bp flanks,
unique ungapped reciprocal coordinate transfer, exact CHM13 REF sequence
agreement, >=95% TE coverage and >=80% deletion coverage. Heterozygous calls
cannot establish absence in an unphased maternal assembly. No record absent
from a VCF becomes a negative label. The exported subset can measure
positive recall/callability, **not precision, FDR or membership calibration**.
Read overlap between benchmark and assembly construction is not excluded.

## Mouse CAST PCR state subset

`extract_mouse_pcr.py` reads Table S2 (`mmc3.xlsx`) and Table S3 (`mmc4.xlsx`)
from [PMID 37228752](https://doi.org/10.1016/j.xgen.2023.100291). The public
[Europe PMC supplementary ZIP](https://www.ebi.ac.uk/europepmc/webservices/rest/PMC10203049/supplementaryFiles)
contains both. Extraction requires `openpyxl`; the subsequent transfer and
evaluation use only the standard library.

```sh
python3 studies/validation/extract_mouse_pcr.py \
  --table-s2 mmc3.xlsx --table-s3 mmc4.xlsx --output truth/mouse-reference
python3 studies/validation/mouse_pcr_truth.py prepare \
  --candidates truth/mouse-reference/candidates.json \
  --fasta c57.primary.fna.gz --te c57.repeats.bed --output truth/mouse.json
python3 studies/validation/mouse_pcr_truth.py evaluate \
  --truth truth/mouse.json --decisions results/tevox.decisions.tsv.gz \
  --output truth/mouse-evaluation.json
```

The table contains 88 rows and 85 unique variant IDs; repeated PCR rows are
preserved in the audit, and conflicting CAST labels fail closed. Eligibility
uses the external TE annotation and observed CAST PCR state before reading
any predictions. GRCm39 reference sequence plus 50 bp flanks must occur
exactly once, unchanged, in the selected C57BL/6NJ assembly. TE and deletion
coverage thresholds are 95% and 80%, respectively. This binds the NJ
reference allele by sequence, while CAST SV status is supported by PCR.
**B6 PCR is C57BL/6J, not C57BL/6NJ.** A missing VCF record or computational
`MERGE_SAMPLES` entry never substitutes for a PCR label. These labels address
biological state, not candidate membership or genome-wide accuracy.

## Archived real-run verification

`verify_completed.py RUN... --output integrity.json` checks recorded output
and input hashes, the 17-table contract and lossless archive bindings. Add
`--roundtrip` to recompute decompressed archive hashes as well. The ordinary
check relies on previously verified uncompressed digests bound to unchanged
compressed sources. `check_evidence_contract.py PREFIX... --output rows.json`
streams identity, mapping and provenance invariants with bounded memory.
It does not replace the full relational schema check used on small fixtures
and the complete Arabidopsis run; restore archived input paths before running
`tests/validate_schema.py` because that validator hashes the original inputs.

## Controlled copy-context and scale evaluation

```sh
python3 studies/validation/synthetic_suite.py --mode polyploid \
  --tevox ./tevox --output simulations/polyploid
python3 studies/validation/synthetic_suite.py --mode scale \
  --tevox ./tevox --output simulations/scale
```

The generator writes ancestry truth before invoking TEvoX. The copy-context
suite uses 128 loci and 2/4/6 retained homeologous copies, with/without missing
annotation, seeds 17/31/53 and with/without a supplied MCScanX prior (36 runs).
These copy counts are modelled homeologous copies, not a measurement of
cytogenetic ploidy. DNA paths are exact and supplied. Seeds change sequence
letters while retaining geometry; they are not independent biological
replicates. One locus overlaps the central anchor, intentionally remaining
in the evaluated set even when its context is not strong enough to merge.

The scale suite uses 100/1,000/10,000/100,000 one-to-one ancestral loci,
totalling 200–200,000 TE annotations. Runtime includes core analysis and
gzip output, but not genome alignment, annotation preparation or later
Python evaluation. It measures a sparse controlled graph; dense all-to-all
or population-scale graphs need additional evaluation. Actual multi-million
TE pair runs provide a separate real-data resource measurement.

## Probability calibration acceptance criteria

The built-in score remains `BUILTIN_UNCALIBRATED_V1`. A calibrated model
requires independently adjudicated **same-locus and different-locus** labels
for the modeled candidate task, a known sampling design, leakage-safe
development/calibration/test partitions, clade/species-pair holdout and an
untruncated candidate-generation contract. Presence/absence labels cannot
automatically label correspondence between two TE annotation fragments.

The existing public helpers prepare and audit these prerequisites:

```sh
./tevox graph ... --max-candidates 0 --max-graph-candidates 0 --output untruncated
python3 scripts/tevox_export_training.py --prefix untruncated \
  --truth truth.candidates.tsv --dataset-id D1 --run-id R1 --output model_data
python3 scripts/tevox_split_audit.py --truth truth.candidates.tsv \
  --locus-truth truth.members.tsv --splits splits.tsv --strict --output split_audit.json
```

Neither those commands nor the controlled simulations establish a validated
calibration model. `results/calibration-readiness.json` records the actual
acceptance status. The reserved model format is described in
`../../docs/MODEL_FORMAT.md`; no fitted model is silently substituted for
the built-in score.
