# build_run_log.ps1 -- assembles run_log.txt from the per-step logs actually written.
param([string]$Root = $(if ($env:SIGLEC15_GSE278688_ROOT) { $env:SIGLEC15_GSE278688_ROOT } else { "E:\Tex\HMlong\BJC_EXT_v1_20260927\phaseB_GSE278688" }))
$out = Join-Path $Root "run_log.txt"

$hdr = @()
$hdr += "# ============================================================================="
$hdr += "# run_log.txt -- Phase B step 1, GSE278688 (FROZEN PRODUCTION RUN)"
$hdr += "# dataset scope : GSE278688 ONLY. No other dataset was accessed."
$hdr += "# generated     : $(Get-Date -Format 'yyyy-MM-dd HH:mm:ss K')"
$hdr += "# host          : Windows 11 x64 build 26200"
$hdr += "# shell         : Windows PowerShell 5.1.26100.9444 (pwsh 7 NOT installed)"
$hdr += "# R             : 4.5.1 (2025-06-13 ucrt)"
$hdr += "# Python        : 3.13.0"
$hdr += "# output root   : $Root"
$hdr += "#"
$hdr += "# Every command below was executed from the session working directory using"
$hdr += "# ABSOLUTE paths. Exit status is recorded per step; 'exit=0' means the runner"
$hdr += "# reported success. Where a step first failed and was then corrected, BOTH the"
$hdr += "# failing and the corrected run are listed -- failures are not hidden."
$hdr += "# ============================================================================="
$hdr += ""
$hdr += "## STEP INDEX (script -> result file -> seed)"
$hdr += ""
$hdr += "  00_audit_inputs.R                -> stdout log only            (no seed)"
$hdr += "  01_audit_design.R                -> stdout log only            (no seed)"
$hdr += "  01b_audit_design_fixed.R         -> stdout log only            (no seed)"
$hdr += "  hash_inputs.ps1                  -> provenance/input_hashes.txt (no seed)"
$hdr += "  parallel_dl32.ps1                -> data/GSE278688_sc_matrix.mtx.gz (no seed)"
$hdr += "  03_load_and_cache.py             -> output/_myeloid_gene_cache.csv,"
$hdr += "                                      output/_allcells_siglec15_cache.csv,"
$hdr += "                                      output/_matrix_summary.json   (no RNG)"
$hdr += "  05_normalisation_fork_test.py    -> output/_normalisation_fork_test.json (no RNG)"
$hdr += "  06_verify_gene_mapping.py        -> stdout log only            (no RNG)"
$hdr += "  07_eligibility_check.R           -> stdout log only            (no seed)"
$hdr += "  08_final_analysis.R              -> output/sample_manifest.csv,"
$hdr += "                                      output/primary_contrast_by_sample.csv,"
$hdr += "                                      output/primary_summary.csv,"
$hdr += "                                      output/specimen_level_summaries.csv,"
$hdr += "                                      output/specimen_level_correlations.csv"
$hdr += "                                      SEED 15072026, 10000 bootstrap resamples"
$hdr += "  09_independent_verification.py   -> output/_verification.json      (no RNG)"
$hdr += "  session_info.R / session_info.py -> sessionInfo.txt"
$hdr += ""
$hdr += "# ============================================================================="
$hdr += "# VERBATIM STEP LOGS"
$hdr += "# ============================================================================="
$hdr += ""

$hdr | Out-File $out -Encoding utf8

$logs = Get-ChildItem (Join-Path $Root "logs") -File | Sort-Object LastWriteTime
$CAP = 262144   # 256 KB: beyond this a log is excerpted, the full file stays in logs/
foreach ($l in $logs) {
  $sz = $l.Length
  "" | Out-File $out -Append -Encoding utf8
  "################################################################################" | Out-File $out -Append -Encoding utf8
"### LOG FILE : $($l.Name)" | Out-File $out -Append -Encoding utf8
"### written  : $($l.LastWriteTime.ToString('yyyy-MM-dd HH:mm:ss'))" | Out-File $out -Append -Encoding utf8
"### size     : $sz bytes" | Out-File $out -Append -Encoding utf8
  if ($sz -gt $CAP) {
    "### NOTE     : log exceeds $CAP bytes; excerpted below (first 400 / last 400 lines)." | Out-File $out -Append -Encoding utf8
    "###            The FULL file is preserved at logs\$($l.Name)" | Out-File $out -Append -Encoding utf8
    "###            The bulk is a per-cell table() dump from the SUPERSEDED first" | Out-File $out -Append -Encoding utf8
    "###            audit; the corrected audit is 01b_audit_design_fixed." | Out-File $out -Append -Encoding utf8
  }
  "################################################################################" | Out-File $out -Append -Encoding utf8
  if ($sz -gt $CAP) {
    $all = Get-Content $l.FullName -Encoding UTF8
    "--- first 400 lines ---" | Out-File $out -Append -Encoding utf8
    $all | Select-Object -First 400 | Out-File $out -Append -Encoding utf8
    "--- ... $($all.Count) lines total; omitted ... ---" | Out-File $out -Append -Encoding utf8
    "--- last 400 lines ---" | Out-File $out -Append -Encoding utf8
    $all | Select-Object -Last 400 | Out-File $out -Append -Encoding utf8
  } else {
    Get-Content $l.FullName -Encoding UTF8 | Out-File $out -Append -Encoding utf8
  }
}
"" | Out-File $out -Append -Encoding utf8
"# =============================================================================" | Out-File $out -Append -Encoding utf8
"# END OF RUN LOG" | Out-File $out -Append -Encoding utf8
"# =============================================================================" | Out-File $out -Append -Encoding utf8
Write-Host "run_log.txt written: $((Get-Item $out).Length) bytes from $($logs.Count) log files"
