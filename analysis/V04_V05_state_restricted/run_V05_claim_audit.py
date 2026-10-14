"""V0.5 Point 5: claim-language audit.

Scan the main manuscript text and supplementary material for wording that tends to cross the evidence boundary, and output line numbers, context and suggested replacements.
Performs text localisation and classification only; does not modify any file.
"""

from __future__ import annotations

import json
import re
import sys
from pathlib import Path

BASE = Path(r".")
OUT = Path(r".\data\analysis\V04_discriminating_analysis")
MAIN = BASE / "Main_Manuscript_v34_BMC_Medicine.tex"
SI = BASE / "Additional_file_1_Supplementary_Material_v34.tex"

# (category, regex, risk note, suggestion)
RULES: list[tuple[str, str, str, str]] = [
    ("causal", r"\bcaus(e|al|ally|ation)\b", "causal verb", "change to associated with / consistent with"),
    ("causal", r"\bdriv(es|en|ing)\b", "driving wording", "change to associated with or delete"),
    ("causal", r"\binduction of\b", "induction wording", "restrict to reported / shown in the source study"),
    ("validation", r"\bvalidat(e|ed|es|ion|ing)\b", "validation-level wording", "change to qualification / concordance testing / orthogonal measurement"),
    ("validation", r"\bqualif(y|ied|ication)\b", "qualification-level wording", "verify that it is used only in the negative sense (does not qualify)"),
    ("validation", r"\bconfirm(s|ed|ing|ation)?\b", "confirmation-level wording", "change to did not confirm / not established"),
    ("replication", r"\breplicat(e|ed|ion|ing)\b", "replication-level wording", "after v0.4 still must be restricted to protocol-based recovery"),
    ("mechanism", r"\bmechanis(m|ms|tic)\b", "mechanism-level wording", "check whether it exceeds the RNA evidence; public RNA cannot support claims of pathway causality"),
    ("mechanism", r"\bcontext-dependent\b", "context-dependent mechanism", "change to context-associated / varies across samples"),
    ("proxy", r"\bprox(y|ies)\b", "proxy wording", "check whether it implies a substitution relationship already validated at patient level"),
    ("strength", r"\bexclud(e|ed|es|ing|ion)\b", "exclusionary wording", "RE precision is limited with three cohorts; change to did not support"),
    ("strength", r"\b(prove|proves|proven|demonstrat(e|ed|es|ing|ion))\b", "proof-level wording", "check whether the object is a directly observable quantity"),
    ("strength", r"\bsufficien(t|tly)\b", "sufficiency wording", "keep only in contexts where discrimination analysis has been performed"),
    ("strength", r"\bnoise\b", "noise attribution", "change to estimation error / sampling variability"),
    ("strength", r"\buninformative\b", "uninformative verdict", "an established operational definition; retain"),
    ("scope", r"\bgenerali[sz]ab(le|ility)\b", "generalisation wording", "check whether it exceeds the PDAC / single-target scope"),
    ("scope", r"\bclinical(ly)? (useful|actionable|valid)\b", "clinical-level wording", "change to hypothesis-generating"),
    ("scope", r"\bpatient selection\b", "patient selection", "check whether it implies validated stratification capability"),
]


def scan(path: Path, stop_at_references: bool) -> list[dict]:
    lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
    hits = []
    for i, line in enumerate(lines, start=1):
        if stop_at_references and re.match(r"\s*\\begin\{thebibliography\}", line):
            break
        if line.strip().startswith("%"):
            continue
        for cat, pat, risk, fix in RULES:
            for m in re.finditer(pat, line, flags=re.I):
                start = max(0, m.start() - 130)
                end = min(len(line), m.end() + 130)
                hits.append({
                    "file": path.name, "line": i, "category": cat,
                    "matched": m.group(0), "risk": risk, "suggested": fix,
                    "context": line[start:end].replace("\n", " "),
                })
    return hits


def main() -> None:
    hits = scan(MAIN, True) + scan(SI, False)
    (OUT / "_v05_claim_language_hits.json").write_text(
        json.dumps(hits, indent=2, ensure_ascii=False), encoding="utf-8")

    by_cat: dict[str, list[dict]] = {}
    for h in hits:
        by_cat.setdefault(h["category"], []).append(h)

    print(f"total hits: {len(hits)}")
    for cat, items in sorted(by_cat.items(), key=lambda kv: -len(kv[1])):
        print(f"\n### {cat}  ({len(items)})")
        for h in items:
            print(f"  {h['file']}:L{h['line']:<4} [{h['matched']}] ...{h['context'][:200]}...")


if __name__ == "__main__":
    sys.stdout.reconfigure(encoding="utf-8")
    main()
