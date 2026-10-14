"""Verify SHA256_MANIFEST.txt against the working tree.

Manifest format, one line per file:

    <64 hex chars><two spaces><relative path with forward slashes>

Lines that are blank or begin with '#' are ignored. Any other line that is not a
well-formed record is a **failure**, not a warning.

Exit status is 0 only when all four of these hold:

  * every listed file is present,
  * every listed file hashes to its recorded digest,
  * every manifest line parses,
  * no file is present that the manifest does not account for.

That last condition is why this script is the release gate: a stray file added
after the manifest was generated must stop a release, otherwise the published
archive would not match its own inventory.

Two classes of path are deliberately outside the inventory and are reported
without failing the run, because the manifest generator excludes them by
construction: archive/temporary files, and Python bytecode caches. They are
listed under "not locked (ignored by design)" so their presence is still visible.

Nothing is written.
"""
from __future__ import annotations

import hashlib
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent
MANIFEST = ROOT / "SHA256_MANIFEST.txt"
CHUNK = 1 << 20

# Must stay in step with the exclusions in the manifest generator. If these two
# lists ever diverge, a clean tree starts failing — or a dirty one starts passing.
IGNORED_SUFFIXES = {".zip", ".7z", ".tar", ".gz", ".tgz", ".rar", ".bz2", ".xz",
                    ".tmp", ".bak", ".orig", ".swp", ".pyc"}
IGNORED_DIRS = {"__pycache__", ".git"}
IGNORED_NAMES = {MANIFEST.name}


def sha256_of(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        while True:
            block = fh.read(CHUNK)
            if not block:
                break
            h.update(block)
    return h.hexdigest()


def ignored_by_design(rel: Path) -> bool:
    return (rel.suffix.lower() in IGNORED_SUFFIXES
            or any(part in IGNORED_DIRS for part in rel.parts)
            or rel.name in IGNORED_NAMES)


def in_skipped_dir(rel: Path) -> bool:
    """True for paths excluded because of a whole directory, not a filename rule.

    Pruned from the walk entirely. Inside a Git working tree .git/ holds thousands
    of objects; listing them one by one buries the real findings.
    """
    return any(part in IGNORED_DIRS for part in rel.parts)


def parse_manifest(path: Path) -> tuple[list[tuple[str, str]], list[str]]:
    """Return (entries, malformed_lines). Malformed lines are returned, not skipped silently."""
    entries: list[tuple[str, str]] = []
    malformed: list[str] = []
    for raw in path.read_text(encoding="utf-8").splitlines():
        if not raw.strip() or raw.lstrip().startswith("#"):
            continue
        digest, sep, rel = raw.partition("  ")
        digest, rel = digest.strip().lower(), rel.strip()
        if sep != "  " or len(digest) != 64 or not rel:
            malformed.append(raw)
            continue
        try:
            int(digest, 16)
        except ValueError:
            malformed.append(raw)
            continue
        entries.append((digest, rel))
    return entries, malformed


def main() -> int:
    if not MANIFEST.exists():
        print(f"FAIL: {MANIFEST.name} not found next to this script")
        return 1

    entries, malformed = parse_manifest(MANIFEST)
    print(f"SHA256_MANIFEST.txt lists {len(entries)} files")

    missing: list[str] = []
    mismatched: list[tuple[str, str, str]] = []
    for digest, rel in entries:
        target = ROOT / rel
        if not target.is_file():
            missing.append(rel)
            continue
        got = sha256_of(target)
        if got != digest:
            mismatched.append((rel, digest, got))

    listed = {rel for _, rel in entries}
    unlisted: list[str] = []
    ignored: list[str] = []
    pruned = 0
    for p in ROOT.rglob("*"):
        if not p.is_file():
            continue
        rel = p.relative_to(ROOT)
        if in_skipped_dir(rel):            # .git/, __pycache__/ -- counted, not listed
            pruned += 1
            continue
        rel_s = rel.as_posix()
        if rel_s in listed:
            continue
        (ignored if ignored_by_design(rel) else unlisted).append(rel_s)
    unlisted.sort()
    ignored.sort()

    for rel in missing:
        print(f"  [missing]    {rel}")
    for rel, want, got in mismatched:
        print(f"  [mismatch]   {rel}\n               want {want}\n               got  {got}")
    for raw in malformed:
        print(f"  [malformed]  {raw!r}   (not a valid manifest record)")
    for rel in unlisted:
        print(f"  [unlisted]   {rel}   (present but not locked)")
    for rel in ignored:
        print(f"  [ignored]    {rel}   (not locked, excluded by design)")

    ok = not missing and not mismatched and not malformed and not unlisted
    print(
        f"\n{'PASS' if ok else 'FAIL'}: "
        f"{len(entries) - len(missing) - len(mismatched)}/{len(entries)} locked files verified"
        + (f"; {len(missing)} missing" if missing else "")
        + (f"; {len(mismatched)} mismatched" if mismatched else "")
        + (f"; {len(malformed)} malformed manifest lines" if malformed else "")
        + (f"; {len(unlisted)} present but unlisted" if unlisted else "")
        + (f"; {len(ignored)} ignored by design" if ignored else "")
        + (f"; {pruned} files under skipped directories" if pruned else "")
    )
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
