# parallel_range_dl.ps1 -- download one URL in N byte ranges concurrently.
# Windows PowerShell 5.1 compatible (Start-Job based).
param(
  [Parameter(Mandatory=$true)][string]$Url,
  [Parameter(Mandatory=$true)][string]$Out,
  [int]$Parts = 12,
  [Parameter(Mandatory=$true)][string]$PartDir
)
$ErrorActionPreference = 'Stop'
New-Item -ItemType Directory -Force -Path $PartDir | Out-Null

$req = [System.Net.HttpWebRequest]::Create($Url)
$req.Method = 'HEAD'; $req.UserAgent = 'Mozilla/5.0'
$resp = $req.GetResponse()
$total = [int64]$resp.ContentLength
$resp.Close()
if ($total -le 0) { throw "no Content-Length for $Url" }

$chunk = [int64][math]::Ceiling($total / $Parts)
Write-Host "total=$total parts=$Parts chunk=$chunk"

$script = {
  param($u, $s, $e, $f)
  $tmp = "$f.tmp"
  for ($try = 1; $try -le 8; $try++) {
    try {
      $hreq = [System.Net.HttpWebRequest]::Create($u)
      $hreq.AddRange($s, $e)
      $hreq.UserAgent = 'Mozilla/5.0'
      $hreq.Timeout = 900000; $hreq.ReadWriteTimeout = 900000
      $hresp = $hreq.GetResponse()
      $in = $hresp.GetResponseStream()
      $out = [System.IO.File]::Create($tmp)
      $buf = New-Object byte[] 1048576
      $got = 0
      while (($n = $in.Read($buf, 0, $buf.Length)) -gt 0) { $out.Write($buf, 0, $n); $got += $n }
      $out.Close(); $in.Close(); $hresp.Close()
      $want = $e - $s + 1
      if ($got -eq $want) { Move-Item -Force $tmp $f; return $true }
    } catch { }
    Start-Sleep -Seconds (2 * $try)
  }
  return $false
}

$jobs = @()
for ($i = 0; $i -lt $Parts; $i++) {
  $start = [int64]($i * $chunk)
  if ($start -ge $total) { break }
  $end = [int64][math]::Min($start + $chunk - 1, $total - 1)
  $f = Join-Path $PartDir ("part_{0:D3}.bin" -f $i)
  $jobs += [pscustomobject]@{
    Index = $i; Start = $start; End = $end; File = $f
    Job = Start-Job -ScriptBlock $script -ArgumentList $Url, $start, $end, $f
  }
}

$null = Wait-Job -Job ($jobs | ForEach-Object { $_.Job })
$results = $jobs | ForEach-Object {
  $ok = Receive-Job -Job $_.Job
  Remove-Job -Job $_.Job -Force
  $have = 0
  if (Test-Path $_.File) { $have = (Get-Item $_.File).Length }
  [pscustomobject]@{ Index = $_.Index; Ok = [bool]$ok; File = $_.File
                     Want = ($_.End - $_.Start + 1); Have = $have }
}
$results | Format-Table -AutoSize
$bad = @($results | Where-Object { -not $_.Ok -or $_.Have -ne $_.Want })
if ($bad.Count -gt 0) { throw "$($bad.Count) parts failed" }

$fs = [System.IO.File]::Create($Out)
foreach ($r in ($results | Sort-Object Index)) {
  $bytes = [System.IO.File]::ReadAllBytes($r.File)
  $fs.Write($bytes, 0, $bytes.Length)
  Remove-Item -Force $r.File
}
$fs.Close()
$final = (Get-Item $Out).Length
Write-Host "wrote $Out bytes=$final expected=$total match=$($final -eq $total)"
if ($final -ne $total) { throw "size mismatch" }
