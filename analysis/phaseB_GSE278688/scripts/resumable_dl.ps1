# resumable_dl.ps1 -- resumable single-stream HTTP download with progress log.
# Usage: resumable_dl.ps1 -Url <u> -Out <path> -Log <logfile>
param(
  [Parameter(Mandatory=$true)][string]$Url,
  [Parameter(Mandatory=$true)][string]$Out,
  [Parameter(Mandatory=$true)][string]$Log
)
function Log($m) { "$(Get-Date -Format 'yyyy-MM-dd HH:mm:ss') $m" | Out-File -FilePath $Log -Append -Encoding utf8 }

# total size
$req = [System.Net.HttpWebRequest]::Create($Url); $req.Method='HEAD'; $req.UserAgent='Mozilla/5.0'
$resp = $req.GetResponse(); $total = [int64]$resp.ContentLength; $resp.Close()
Log "URL=$Url TOTAL=$total OUT=$Out"

# normalize: if a .tmp exists keep it, else if Out exists and shorter than total, rename to .tmp
$tmp = "$Out.tmp"
if (-not (Test-Path $tmp)) {
  if ((Test-Path $Out) -and ((Get-Item $Out).Length -lt $total)) { Move-Item -Force $Out $tmp }
}
$have = 0
if (Test-Path $tmp) { $have = [int64](Get-Item $tmp).Length }
Log "resume from byte $have"

$sw = [System.Diagnostics.Stopwatch]::StartNew()
$lastLog = 0
while ($have -lt $total) {
  try {
    $hreq = [System.Net.HttpWebRequest]::Create($Url)
    if ($have -gt 0) { $hreq.AddRange($have, $total - 1) }
    $hreq.UserAgent = 'Mozilla/5.0'; $hreq.Timeout = 900000; $hreq.ReadWriteTimeout = 900000
    $hresp = $hreq.GetResponse()
    $in = $hresp.GetResponseStream()
    $fs = [System.IO.File]::Open($tmp, [System.IO.FileMode]::Append, [System.IO.FileAccess]::Write)
    $buf = New-Object byte[] 1048576
    while (($n = $in.Read($buf, 0, $buf.Length)) -gt 0) {
      $fs.Write($buf, 0, $n); $have += $n
      if (($sw.Elapsed.TotalSeconds - $lastLog) -ge 30) {
        $lastLog = $sw.Elapsed.TotalSeconds
        $rate = 0.0
        if ($sw.Elapsed.TotalSeconds -gt 0) { $rate = $have / 1KB / $sw.Elapsed.TotalSeconds }
        Log ("progress {0:N0}/{1:N0} ({2:N1}%) rate={3:N1} KB/s elapsed={4:N1} min" -f $have,$total,(100*$have/$total),$rate,$sw.Elapsed.TotalMinutes)
      }
    }
    $fs.Close(); $in.Close(); $hresp.Close()
    if ($have -ge $total) { break }
    Log "stream ended early at $have; reconnecting"
  } catch {
    Log "ERROR: $($_.Exception.Message) -- retrying"
  }
  Start-Sleep -Seconds 3
}
Move-Item -Force $tmp $Out
Log ("DONE bytes={0} elapsed={1:N1} min avg={2:N1} KB/s" -f (Get-Item $Out).Length, $sw.Elapsed.TotalMinutes, ((Get-Item $Out).Length/1KB/$sw.Elapsed.TotalSeconds))
