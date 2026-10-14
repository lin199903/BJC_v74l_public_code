# Release gate for the SIGLEC15 cross-scale measurement code archive.
#
# Order of operations, deliberately:
#   1. check that every required metadata file is present
#   2. re-verify SHA256_MANIFEST.txt            -- nothing is published unless this passes
#   3. check the tree the manifest locks against the tree Git will actually publish
#   4. only then create the git tag and the GitHub release
#
# On exit status: `$ErrorActionPreference = "Stop"` does NOT turn a non-zero exit
# from a native program (git, gh, python) into a terminating error in Windows
# PowerShell 5.1. Every native call below is therefore routed through
# Invoke-Native, which checks $LASTEXITCODE and throws. Without that, a failed
# `git push` or `gh release create` would fall through to the "Done" line and
# report success for a release that never happened.
#
# Zenodo's GitHub integration only archives releases published *after* the
# integration is switched on. Enable it at https://zenodo.org/account/settings/github/
# before running this script without -DryRun, or the release will not be minted
# a DOI and will have to be deleted and re-issued.
#
# Two metadata fields describe the release itself and are NOT checked here,
# because they are only correct once the release exists -- reconcile them by hand
# immediately before running this script. See the "Before publishing" section of
# RELEASE_NOTES.md:
#   * CITATION.cff  date-released  -> must equal the tag date
#   * .zenodo.json  version        -> must equal $Tag exactly
#
# Usage:
#   .\release.ps1 -DryRun        # verify only, change nothing
#   .\release.ps1                # verify, then tag and publish v1.0

param(
    [string]$Repo    = "lin199903/BJC_v74l_public_code",
    [string]$Tag     = "v1.0",
    [string]$Title   = "SIGLEC15 cross-scale measurement - v1.0 frozen analysis code",
    [switch]$DryRun
)

$ErrorActionPreference = "Stop"

function Invoke-Native {
    <#
      Run a native command and throw if it exits non-zero. Returns its output.
      $LASTEXITCODE is set by the native call itself, so it must be read before
      anything else runs.
    #>
    param(
        [Parameter(Mandatory)][string]$What,
        [Parameter(Mandatory)][scriptblock]$Command
    )
    $output = & $Command 2>&1
    $code = $LASTEXITCODE
    if ($code -ne 0) {
        throw "$What failed (exit $code)`n$($output | Out-String)"
    }
    return $output
}

function Invoke-NativeLines {
    <#
      Run a native command and return its non-empty output lines as an ARRAY.

      Always an array, never a bare string, so that callers can index [0] and use
      -contains safely. `@(Invoke-Native ...) | Where-Object {...}` does NOT
      guarantee that: the pipeline unwraps a single surviving element back to a
      scalar, after which [0] silently returns the first *character*.
    #>
    param(
        [Parameter(Mandatory)][string]$What,
        [Parameter(Mandatory)][scriptblock]$Command
    )
    $output = Invoke-Native -What $What -Command $Command
    # The leading comma is load-bearing: PowerShell unwraps a single-element array
    # on return, so `return @(...)` would hand back a bare string and the caller's
    # [0] would silently index the first character. `return ,@(...)` wraps once so
    # that the unwrap yields the array itself.
    return ,@($output | Where-Object { $_ -and "$_".Trim() } | ForEach-Object { "$_".Trim() })
}

$Root = $PSScriptRoot
if (-not $Root) { $Root = (Get-Location).Path }
Push-Location $Root
try {

    # ---- 1. required metadata ------------------------------------------------
    $required = @("README.md", "CITATION.cff", "LICENSE", ".zenodo.json",
                  "SHA256_MANIFEST.txt", "verify_manifest.py", "RELEASE_NOTES.md")
    foreach ($f in $required) {
        if (-not (Test-Path (Join-Path $Root $f))) { throw "Missing required file: $f" }
    }
    Write-Host "Metadata files present: $($required.Count)/$($required.Count)"

    # ---- 2. manifest verification -------------------------------------------
    Write-Host "== manifest verification =="
    Invoke-Native -What "verify_manifest.py" -Command { python verify_manifest.py } | Out-Host

    # No @() around the pipeline: Measure-Object emits a single object, so wrapping
    # the pipeline in an array makes .Count report 1 regardless of how many lines matched.
    $locked = (Select-String -Path (Join-Path $Root "SHA256_MANIFEST.txt") `
                              -Pattern '^[0-9a-f]{64}  ' | Measure-Object).Count
    Write-Host "Manifest verification PASS ($locked files locked)"

    # ---- 3. the tree Git will actually publish -------------------------------
    # Scope first. `git ls-files` run from a subdirectory reports only the tracked
    # files below that subdirectory, while a tag, a GitHub Release and the Zenodo
    # source archive all cover the WHOLE repository. If this directory were a
    # subdirectory of a larger repository, the file-set comparison below would pass
    # and the release would still carry everything else in that repository. Require
    # the archive root and the Git root to be the same directory.
    Write-Host "== repository scope check =="
    $toplevel = Invoke-NativeLines -What "git rev-parse --show-toplevel" `
                                   -Command { git rev-parse --show-toplevel }
    if ($toplevel.Count -ne 1) {
        throw "Cannot determine a unique Git repository root (got $($toplevel.Count) lines)"
    }
    $gitRoot  = [System.IO.Path]::GetFullPath([string]$toplevel[0]).TrimEnd('\', '/')
    $codeRoot = [System.IO.Path]::GetFullPath($Root).TrimEnd('\', '/')
    if (-not [string]::Equals($gitRoot, $codeRoot, [System.StringComparison]::OrdinalIgnoreCase)) {
        throw ("The code archive must be the Git repository root, because the tag and the" + "`n" +
               "Zenodo source archive cover the whole repository and not just this directory." + "`n" +
               "  archive root : $codeRoot" + "`n" +
               "  git root     : $gitRoot")
    }
    Write-Host "Repository root is this directory"

    # verify_manifest.py ignores archive and cache paths by design, so a committed
    # .zip would pass it and still reach GitHub and Zenodo. Compare against what
    # Git tracks instead of trusting the working directory.
    Write-Host "== tracked-file check =="
    $tracked = @(Invoke-NativeLines -What "git ls-files" -Command { git ls-files } |
        ForEach-Object { $_ -replace '\\', '/' })

    $listed = @(Select-String -Path (Join-Path $Root "SHA256_MANIFEST.txt") `
                              -Pattern '^[0-9a-f]{64}  (.+)$' |
                ForEach-Object { $_.Matches[0].Groups[1].Value.Trim() })
    $allowed = $listed + @("SHA256_MANIFEST.txt")

    $extra = @($tracked | Where-Object { $allowed -notcontains $_ })
    $absent = @($allowed | Where-Object { $tracked -notcontains $_ })
    if ($extra.Count -gt 0) {
        throw "Git tracks $($extra.Count) file(s) the manifest does not lock:`n  $($extra -join "`n  ")"
    }
    if ($absent.Count -gt 0) {
        throw "$($absent.Count) locked file(s) are not tracked by Git:`n  $($absent -join "`n  ")"
    }
    Write-Host "Tracked files match the manifest exactly ($($tracked.Count) files)"

    if ($DryRun) {
        Write-Host "`n-DryRun: stopping before tag creation. Nothing was published."
        return
    }

    # ---- 4. tag and release --------------------------------------------------
    $dirty = @(Invoke-Native -What "git status" -Command { git status --porcelain })
    if ($dirty.Count -gt 0) {
        throw "Working tree is dirty; commit everything before releasing:`n  $($dirty -join "`n  ")"
    }

    Write-Host "== creating tag $Tag =="
    Invoke-Native -What "git tag -a $Tag" -Command { git tag -a $Tag -m $Title } | Out-Null

    Invoke-Native -What "git push origin $Tag" -Command { git push origin $Tag } | Out-Host

    Write-Host "== creating GitHub release =="
    Invoke-Native -What "gh release create" -Command {
        gh release create $Tag --repo $Repo --title $Title --notes-file RELEASE_NOTES.md
    } | Out-Host

    Write-Host "`nDone. Confirm the Zenodo deposit and copy the *version* DOI"
    Write-Host "(10.5281/zenodo.XXXXXXX, not the concept DOI) into the manuscript."
}
finally {
    Pop-Location
}
