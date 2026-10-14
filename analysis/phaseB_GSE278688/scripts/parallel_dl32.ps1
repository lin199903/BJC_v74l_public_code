# parallel_dl32.ps1 -- 32-way ranged download with retries and progress log.
param(
  [Parameter(Mandatory=$true)][string]$Url,
  [Parameter(Mandatory=$true)][string]$Out,
  [Parameter(Mandatory=$true)][string]$PartDir,
  [Parameter(Mandatory=$true)][string]$Log,
  [int]$Parts = 32
)
$ErrorActionPreference = 'Continue'
New-Item -ItemType Directory -Force -Path $PartDir | Out-Null
function Log($m) { "$(Get-Date -Format 'yyyy-MM-dd HH:mm:ss') $m" | Out-File -FilePath $Log -Append -Encoding utf8 }

$req = [System.Net.HttpWebRequest]::Create($Url); $req.Method='HEAD'; $req.UserAgent='Mozilla/5.0'
$resp = $req.GetResponse(); $total = [int64]$resp.ContentLength; $resp.Close()
$chunk = [int64][math]::Ceiling($total / $Parts)
Log "START total=$total parts=$Parts chunk=$chunk partdir=$PartDir"

$worker = {
  param($u, $s, $e, $f)
  $want = $e - $s + 1
  if ((Test-Path $f) -and ((Get-Item $f).Length -eq $want)) { return 0 }   # already complete
  $tmp = "$f.tmp"
  for ($try = 1; $try -le 40; $try++) {
    try {
      $have = 0
      if (Test-Path $tmp) { $have = [int64](Get-Item $tmp).Length }
      if ($have -ge $want) { Move-Item -Force $tmp $f; return 0 }
      $hreq = [System.Net.HttpWebRequest]::Create($u)
      $hreq.AddRange([int64]($s + $have), $e)
      $hreq.UserAgent = 'Mozilla/5.0'; $hreq.Timeout = 600000; $hreq.ReadWriteTimeout = 600000
      $hresp = $hreq.GetResponse()
      $in = $hresp.GetResponseStream()
      $mode = [System.IO.FileMode]::Append
      if ($have -eq 0) { $mode = [System.IO.FileMode]::Create }
      $out = [System.IO.File]::Open($tmp, $mode, [System.IO.FileAccess]::Write)
      $buf = New-Object byte[] 524288
      while (($n = $in.Read($buf, 0, $buf.Length)) -gt 0) { $out.Write($buf, 0, $n); $have += $n }
      $out.Close(); $in.Close(); $hresp.Close()
      if ($have -ge $want) { Move-Item -Force $tmp $f; return 0 }
    } catch { }
    Start-Sleep -Seconds ([math]::Min(30, 2 * $try))
  }
  return 1
}

$jobs = @()
for ($i = 0; $i -lt $Parts; $i++) {
  $start = [int64]($i * $chunk); if ($start -ge $total) { break }
  $end = [int64][math]::Min($start + $chunk - 1, $total - 1)
  $f = Join-Path $PartDir ("part_{0:D3}.bin" -f $i)
  $jobs += [pscustomobject]@{ Index=$i; Start=$start; End=$end; File=$f
    Job = Start-Job -ScriptBlock $worker -ArgumentList $Url, $start, $end, $f }
}
Log "launched $($jobs.Count) workers"

# progress heartbeat while waiting
$sw = [System.Diagnostics.Stopwatch]::StartNew()
while ($true) {
  $done = @($jobs | Where-Object { $_.Job.State -in @('Completed','Failed','Stopped') }).Count
  if ($done -eq $jobs.Count) { break }
  $bytes = 0
  foreach ($j in $jobs) { if (Test-Path $j.File) { $bytes += (Get-Item $j.File).Length }
                          elseif (Test-Path "$($j.File).tmp") { $bytes += (Get-Item "$($j.File).tmp").Length } }
  Log ("progress {0:N0}/{1:N0} ({2:N1}%) parts_done={3}/{4} rate={5:N1} KB/s elapsed={6:N1} min" -f `
        $bytes,$total,(100*$bytes/$total),$done,$jobs.Count,($bytes/1KB/$sw.Elapsed.TotalSeconds),$sw.Elapsed.TotalMinutes)
  Start-Sleep -Seconds 60
}
$null = Wait-Job -Job ($jobs | ForEach-Object { $_.Job })
$rc = $jobs | ForEach-Object { Receive-Job -Job $_.Job }
$jobs | Remove-Job -Force
$failed = @($rc | Where-Object { $_ -ne 0 }).Count
Log "workers finished, failed=$failed"

# verify every part
$bad = 0
foreach ($j in $jobs) {
  $want = $j.End - $j.Start + 1
  $have = 0; if (Test-Path $j.File) { $have = (Get-Item $j.File).Length }
  if ($have -ne $want) { Log "PART BAD index=$($j.Index) have=$have want=$want"; $bad++ }
}
if ($bad -gt 0) { Log "ABORT: $bad bad parts"; exit 1 }

# concatenate
$fs = [System.IO.File]::Create($Out)
foreach ($j in ($jobs | Sort-Object Index)) {
  $in = [System.IO.File]::OpenRead($j.File)
  $in.CopyTo($fs, 4194304); $in.Close()
  Remove-Item -Force $j.File
}
$fs.Close()
$fin = (Get-Item $Out).Length
Log ("DONE bytes={0} expected={1} match={2} elapsed={3:N1} min" -f $fin,$total,($fin -eq $total),$sw.Elapsed.TotalMinutes)
if ($fin -ne $total) { exit 1 }
