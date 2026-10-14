# build_run_log.ps1 -- assembles run_log.txt for Phase B step 2, GSE278689.
# Concatenates every raw log in logs\ in timestamp order, with a header that
# states the scope and the frozen constants, so the run log is the full command
# and output trail rather than a summary of it.
param([string]$Root = $(if ($env:SIGLEC15_GSE278689_ROOT) { $env:SIGLEC15_GSE278689_ROOT } else { "E:\Tex\HMlong\BJC_EXT_v1_20260927\phaseB_GSE278689" }))
$out = Join-Path $Root "run_log.txt"
$L = New-Object System.Collections.Generic.List[string]

$L.Add("=" * 100)
$L.Add("run_log.txt -- Phase B step 2, GSE278689 (snRNA) -- CROSS-MODALITY SENSITIVITY ANALYSIS")
$L.Add("=" * 100)
$L.Add("assembled            : $(Get-Date -Format 'yyyy-MM-dd HH:mm:ss zzz')")
$L.Add("root                 : $Root")
$L.Add("dataset              : GSE278689 -- snRNA (single NUCLEI), SubSeries of GSE278694")
$L.Add("scope                : GSE278689 ONLY. This is a modality-specific sensitivity analysis.")
$L.Add("                       It is NOT pooled with the manuscript's 19-sample scRNA estimate")
$L.Add("                       and NOT pooled with GSE278688. No cross-dataset combination is")
$L.Add("                       computed anywhere in this run.")
$L.Add("")
$L.Add("FROZEN CONSTANTS (unchanged from the brief)")
$L.Add("  eligible tumour label : tissue == 'Tumor' (mapping fixed before expression-based")
$L.Add("                          eligibility was evaluated; see 02_audit_design.R)")
$L.Add("  cell type             : deposited all_celltype == 'Myeloid'  (NO re-annotation)")
$L.Add("  minimum cells         : >= 50 included tumour-myeloid cells per sample")
$L.Add("  both groups present   : >= 1 detectable AND >= 1 non-detectable  (NO other minimum;")
$L.Add("                          the >=10/>=2/>=5 rule is a SEPARATE restricted analysis)")
$L.Add("  programme score       : ACP5 + MMP9, both required; each z-standardised within the")
$L.Add("                          biological sample (ddof=0), then averaged")
$L.Add("  target in score       : SIGLEC15 NEVER enters the programme score")
$L.Add("  detection             : non-zero RAW deposited count (invariant to normalisation)")
$L.Add("  normalisation         : CP10K then log1p (deposit verified as raw integer counts)")
$L.Add("  contrast              : mean(score|detectable) - mean(score|non-detectable)")
$L.Add("  weighting             : samples equally weighted")
$L.Add("  seed / resamples      : 15072026 / 10000")
$L.Add("")
$L.Add("INTERPRETATION BOUNDARY (pre-registered by the author before this run's result)")
$L.Add("  1. GSE278689 cannot enter the primary estimate.")
$L.Add("  2. It cannot be merged with GSE278688.")
$L.Add("  3. It answers only: does the target-state relation remain estimable, and in what")
$L.Add("     direction, under a different transcript-capture process?")
$L.Add("  Permitted phrasings used in the reports: (a) 'the relation was difficult to")
$L.Add("  estimate under a lower-detectability nuclear RNA measurement process'; (b) 'a")
$L.Add("  positive contrast was observed in a modality-specific sensitivity analysis'.")
$L.Add("")
$L.Add("LANGUAGE NOTE. The four words named in the author's pre-registered interpretation")
$L.Add("criteria are not applied to any sample, estimate, interval or contrast in this run's")
$L.Add("reports. The script transcripts reproduced verbatim below are the raw evidence record")
$L.Add("and are NOT edited: where they print PASS / FAIL tokens, or where an older script")
$L.Add("string names those four words, the tokens are technical verdicts on VERIFICATION")
$L.Add("CONTROLS (does a deliberately wrong convention move the answer? does a marker gene")
$L.Add("land in its expected cell type?). They are not, and must not be read as,")
$L.Add("characterisations of any study result. The two authored report documents")
$L.Add("(results_summary.md, failed_or_excluded.md) carry no such characterisation.")
$L.Add("")
$L.Add("ENVIRONMENT")
$L.Add((Get-Content (Join-Path $Root "sessionInfo.txt") -Raw))
$L.Add("")
$L.Add("INPUT FILES AND HASHES")
$L.Add((Get-Content (Join-Path $Root "provenance\input_hashes.txt") -Raw))
$L.Add("")
$L.Add("=" * 100)
$L.Add("COMMAND AND OUTPUT TRAIL")
$L.Add("Each block below is one executed script. The script file that produced it is named")
$L.Add("in the block header, and every result file is paired with the script that wrote it.")
$L.Add("=" * 100)

function Add-Block($title, $scriptRel, $logFile) {
  $L.Add("")
  $L.Add("-" * 100)
  $L.Add("BLOCK  : $title")
  if ($scriptRel) { $L.Add("SCRIPT : $scriptRel") }
  $L.Add("LOG    : $logFile")
  if (Test-Path $logFile) {
    $L.Add("MTIME  : $((Get-Item $logFile).LastWriteTime.ToString('yyyy-MM-dd HH:mm:ss'))")
    $L.Add("-" * 100)
    $L.Add((Get-Content $logFile -Raw))
  } else {
    $L.Add("-" * 100)
    $L.Add("(log not found)")
  }
}

$map = @(
  @{t="00 input audit";                              s="scripts\00_audit_inputs.R";                 p="00_audit_inputs_*.log"},
  @{t="02 design audit (frozen rules -> this deposit)"; s="scripts\02_audit_design.R";              p="02_audit_design_*.log"},
  @{t="03 load deposit and build caches";            s="scripts\03_load_and_cache.py";              p="03_load_and_cache_*.log"},
  @{t="05 normalisation branch + fork test (CORRECTED)"; s="scripts\05_normalisation_fork_test.py"; p="05_normalisation_fork_FIXED_*.log"},
  @{t="06 gene->row mapping verification";           s="scripts\06_verify_gene_mapping.py";         p="06_verify_gene_mapping_*.log"},
  @{t="07 eligibility and sample manifest";          s="scripts\07_eligibility_check.R";            p="07_eligibility_check_*.log"},
  @{t="08 final analysis (primary contrast, summary, correlations)"; s="scripts\08_final_analysis.R"; p="08_final_analysis_*.log"},
  @{t="09 independent verification vs the manuscript's own functions"; s="scripts\09_independent_verification.py"; p="09_independent_verification_*.log"},
  @{t="10 bootstrap CI sanity";                      s="scripts\10_ci_sanity.R";                    p="10_ci_sanity_*.log"},
  @{t="11 aggregation convention check";             s="scripts\11_aggregation_convention.R";       p="11_aggregation_convention_*.log"},
  @{t="download from GEO FTP";                       s="(logs\download_*.log)";                     p="download_*.log"},
  @{t="input hashing";                               s="scripts\hash_inputs.ps1";                   p="hash_inputs*.log"}
)
foreach ($m in $map) {
  $files = Get-ChildItem (Join-Path $Root "logs") -Filter $m.p -ErrorAction SilentlyContinue |
           Sort-Object LastWriteTime
  if ($files) { Add-Block $m.t $m.s $files[-1].FullName }
}

$L.Add("")
$L.Add("=" * 100)
$L.Add("RESULT FILES AND THE SCRIPT THAT WROTE EACH ONE")
$L.Add("=" * 100)
$pairs = @(
  @("output\sample_manifest.csv",                    "scripts\07_eligibility_check.R"),
  @("output\exclusion_arithmetic.csv",               "scripts\07_eligibility_check.R"),
  @("output\programme_gene_variance_by_sample.csv",  "scripts\07_eligibility_check.R"),
  @("output\primary_contrast_by_sample.csv",         "scripts\08_final_analysis.R"),
  @("output\primary_summary.csv",                    "scripts\08_final_analysis.R"),
  @("output\specimen_level_summaries.csv",           "scripts\08_final_analysis.R"),
  @("output\specimen_level_correlations.csv",        "scripts\08_final_analysis.R"),
  @("output\detectability_comparability.csv",        "scripts\08_final_analysis.R"),
  @("output\ci_sanity_by_seed.csv",                  "scripts\10_ci_sanity.R"),
  @("output\_matrix_summary.json",                   "scripts\03_load_and_cache.py"),
  @("output\_myeloid_gene_cache.csv",                "scripts\03_load_and_cache.py"),
  @("output\_allcells_siglec15_cache.csv",           "scripts\03_load_and_cache.py"),
  @("output\_normalisation_fork_test.json",          "scripts\05_normalisation_fork_test.py"),
  @("output\_verification.json",                     "scripts\06_verify_gene_mapping.py"),
  @("output\_verification_independent.json",         "scripts\09_independent_verification.py"),
  @("provenance\sessionInfo.txt",                    "scripts\session_info.py"),
  @("provenance\sample_accession_map.csv",           "scripts\session_info.py"),
  @("provenance\input_hashes.txt",                   "scripts\hash_inputs.ps1"),
  @("results_summary.md",                            "(written by the agent from the above)"),
  @("failed_or_excluded.md",                         "(written by the agent from the above)"),
  @("run_log.txt",                                   "scripts\build_run_log.ps1")
)
$L.Add(("{0,-52}{1,-46}{2,14}" -f "result file","producing script","bytes"))
foreach ($p in $pairs) {
  $f = Join-Path $Root $p[0]
  $sz = if (Test-Path $f) { (Get-Item $f).Length } else { 0 }
  $L.Add(("{0,-52}{1,-46}{2,14}" -f $p[0], $p[1], $sz))
}
$L.Add("")
$L.Add("SEED: 15072026 is the seed for every bootstrap in this run (08_ and 10_).")
$L.Add("     09_ uses numpy default_rng(15072026) as a second, independent stream.")
$L.Add("     Exact-bootstrap enumeration in 10_ uses no RNG at all.")
$L.Add("END OF RUN LOG")

$L -join "`r`n" | Out-File $out -Encoding utf8
"run_log.txt written: $((Get-Item $out).Length) bytes"
