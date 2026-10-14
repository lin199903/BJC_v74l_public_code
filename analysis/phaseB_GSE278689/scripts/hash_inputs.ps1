# hash_inputs.ps1 -- SHA-256 of every input file used by this run.
# Writes provenance/input_hashes.txt (append-only, timestamped).
param([string]$Root = $(if ($env:SIGLEC15_GSE278689_ROOT) { $env:SIGLEC15_GSE278689_ROOT } else { "E:\Tex\HMlong\BJC_EXT_v1_20260927\phaseB_GSE278689" }))
$out = Join-Path $Root "provenance\input_hashes.txt"
"# input hashes  $(Get-Date -Format 'yyyy-MM-dd HH:mm:ss zzz')" | Out-File $out -Append -Encoding utf8
$files = @()
$files += Get-ChildItem (Join-Path $Root "data") -File -ErrorAction SilentlyContinue
$files += Get-Item (Join-Path $(if ($env:SIGLEC15_PHASEB_ROOT) { $env:SIGLEC15_PHASEB_ROOT } else { "E:\Tex\HMlong\BJC_EXT_v1_20260927" }) "annotations\GSE278689_sn_metadata.csv.gz") -ErrorAction SilentlyContinue
foreach ($f in $files) {
  $h = Get-FileHash -Algorithm SHA256 -Path $f.FullName
  $md5 = (Get-FileHash -Algorithm MD5 -Path $f.FullName).Hash.ToLower()
  ("{0}  sha256={1}  md5={2}  bytes={3}  {4}" -f $f.Name, $h.Hash.ToLower(), $md5, $f.Length, $f.FullName) |
    Out-File $out -Append -Encoding utf8
}
Get-Content $out
