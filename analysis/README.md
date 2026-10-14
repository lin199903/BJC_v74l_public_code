# Versioned analysis archive - v72

> **Historical archive description — read with this in mind.**
>
> This file describes the *full* v72 analysis archive, including its numeric
> result tables, draw-level outputs and integrity manifests. This Zenodo deposit
> is the **code** deposit: it carries the scripts, the frozen rule documents,
> the run and provenance logs, and the software-environment records, but not the
> result tables. Where the description below refers to a per-sample table, a
> draw-level file, a JSON analysis record or an archive manifest, that file is
> published with the paper as part of Additional file 2 rather than here.
>
> See `../README.md` for the artifact map that pairs each analysis with its
> code location and its result location, and `CURRENT_READING_GUIDE.md` for the
> reporting corrections that apply when reading historical status headings.

## v72 tumour-level burden and sparsity extension

The v72 submission retains the complete v39 payload below and adds a
source-locked extension under `v72_burden/`. The extension separates the
within-sample cellular contrast from a tumour-sample myeloid ACP5/MMP9
aggregate burden and tests whether SIGLEC15 RNA preserves the ordering of the
latter.

The extension contains:

- the frozen burden-estimand and sparsity contract, including its SHA-256 lock;
- the analysis script and input-file hash manifest;
- one row per biological sample with the burden and four target readouts;
- primary, sparse-cell, secondary-readout and leave-one-dataset-out estimates;
- all 10,000 cell-exchangeability null draws and their summary;
- a JSON analysis and software record; and
- `V72_PAYLOAD_SHA256.csv`, which inventories every added file except itself.

The sample-aggregate point estimate had already appeared in the v39 archive as
hypothesis-informing evidence. The v72 contract was frozen before the new
sparsity, readout and cell-exchangeability analyses were run. The extension
reproduces all 20 archived sample keys, the SIGLEC15 aggregate to a maximum
absolute difference of 9.8e-17, the ACP5/MMP9 burden to 1.3e-15 and the pooled
point estimate exactly. It is not a self-contained redistribution of the three
public H5AD inputs; their paths, byte sizes and SHA-256 hashes are recorded in
the input manifest.

## Preserved v39 archive

**Archive v39 corresponds to the BMC Medicine manuscript v39.**

The v39 manuscript adds pre-submission review responses: an effect-size
detection boundary, a phased-framework position, a paired-validation design
sketch and provenance chronology. Analysis payload is identical to v38.

`PROVENANCE_TIMELINE.md` records the working file-system modification times of
the protocol, script and result files, which archive packaging would otherwise
reset.

The v38 manuscript applies reviewer-simulation corrections (a supplementary
cross-reference repair, a sensitivity-boundary disclosure added to the main
text, negative-sign normalisation, and a cross-referencing note between Tables
S24 and S25); it adds, removes or alters no analysis. The payload below is
identical to archive v37.

The v37 manuscript carries claim-ceiling corrections and a short biological
framing paragraph; it adds, removes or alters no analysis. The payload below is
identical to archive v36.

The v36 manuscript restructured the presentation of these analyses; it does not
add, remove or alter any analysis. The payload below is identical to archive v35.

The v37 archive preserves the locked primary estimand, analysis rules and random seeds
recorded in v34. It carries forward the PACA-CA detectability provenance, the complete
200-replicate CD274 downsampling file and the S24 fine-myeloid-subtype sensitivity
analysis, and adds the post hoc state-restricted sensitivity analysis reported in
Supplementary Table S25.

## What this archive adds

Under `V04_V05_state_restricted/`:

- `V04_ANALYSIS_SPECIFICATION.md` - frozen design, estimands, eligibility rules and verdict rules
- `V04_STATE_ANNOTATION_PROVENANCE.csv` - per-dataset annotation provenance and barcode-mapping audit
- `V04_GATE_EVALUATION.json` - annotation-recoverability gate outcome, marker-module definition and expression layer used
- `V04_P1A_LOCALISATION.csv` - per-sample SIGLEC15-detectable fraction within TREM2_TAM versus other myeloid cells
- `V04_P1A_WITHIN_STATE.csv` - per-sample within-TREM2_TAM ACP5/MMP9 contrast, depth-adjusted and depth-matched
- `V04_P1B_AOK.csv` - per-sample abundance, occupancy and coupling values (not reported in the manuscript)
- `V04_MECHANISM_ADJUDICATION.md` - adjudication of competing explanations, including all limitations
- `V05_CLAIM_CEILING_AUDIT.md` - wording audit and provenance checks performed before submission
- `V05_ESTIMABILITY_SELECTION.csv` - estimability status by sample with the reason for non-estimability
- `V05_WITHIN_STATE_GENE_DECOMPOSITION.csv` - per-sample ACP5 and MMP9 decomposition inside TREM2_TAM
- `V04_FIG_MECHANISM_DISCRIMINATION.pdf` - summary figure for the state-restricted analyses
- analysis scripts and machine-readable intermediate results

## Standing notes

- The post hoc state-restricted analysis is labelled post hoc throughout. It is not part
  of the locked protocol and must not be described as prespecified relative to the primary
  cross-scale results.
- The nine-class myeloid stratification is a deterministic marker-module construct
  generated before the cross-scale analysis. Its marker definitions exclude SIGLEC15,
  ACP5 and MMP9. It is not an externally validated cell identity.
- The PACA-CA extraction yielded downloadable expression files for 43 of 317 listed donors
  on 6 June 2026; 33 of those 43 had SIGLEC15 FPKM greater than zero. This subset is used
  descriptively only.
- The abundance-matched CD274 analysis is represented by all 200 values in
  `results/CD274_downsampling/downsampling_replicates.csv`; the median Hallmark PI3K
  mean-z correlation is 0.5408759598 (reported as 0.541).

Earlier versions: v24, v27, v28.1, v33 and v34 (frozen). This archive supersedes them for
submission purposes.

## Scope of this archive

This archive records the analyses as they were run. It contains the scripts and the numerical outputs, not a self-contained reproduction package: some scripts read intermediate files from the directory layout in use at the time, and those intermediates are not redistributed here. The scripts compile as distributed; scripts whose inputs are the deposited source data run directly, whereas those that depend on intermediate products require the reader to supply an equivalent layout. The primary inputs are public (GEO accessions listed in the manuscript), and the numerical outputs each script produced are included.
