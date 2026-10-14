# hash_inputs.ps1 -- SHA256 of every deposited input file. Provenance.
param([string]$Root = $(if ($env:SIGLEC15_GSE278688_ROOT) { $env:SIGLEC15_GSE278688_ROOT } else { "E:\Tex\HMlong\BJC_EXT_v1_20260927\phaseB_GSE278688" }))
$out = Join-Path $Root "provenance\input_hashes.txt"
"# SHA256 of GSE278688 deposit files -- generated $(Get-Date -Format 'yyyy-MM-dd HH:mm:ss K')" | Out-File $out -Encoding utf8
"# source: https://ftp.ncbi.nlm.nih.gov/geo/series/GSE278nnn/GSE278688/suppl/" | Out-File $out -Append -Encoding utf8
Get-ChildItem (Join-Path $Root "data") -File | Sort-Object Name | ForEach-Object {
  $h = (Get-FileHash $_.FullName -Algorithm SHA256).Hash
  $line = "{0}  {1,15}  {2}" -f $h, $_.Length, $_.Name
  $line | Out-File $out -Append -Encoding utf8
  Write-Host $line
}
# also the frozen annotation copy, to prove it is byte-identical to the deposit
$frozen = Join-Path $(if ($env:SIGLEC15_PHASEB_ROOT) { $env:SIGLEC15_PHASEB_ROOT } else { "E:\Tex\HMlong\BJC_EXT_v1_20260927" }) "annotations\GSE278688_sc_metadata.csv.gz"
if (Test-Path $frozen) {
  $h = (Get-FileHash $frozen -Algorithm SHA256).Hash
  $line = "{0}  {1,15}  {2}" -f $h, (Get-Item $frozen).Length, "annotations/GSE278688_sc_metadata.csv.gz"
  $line | Out-File $out -Append -Encoding utf8
  Write-Host $line
}
"AUDIT-ONLY: hashing does not modify the files." | Out-File $out -Append -Encoding utf8
