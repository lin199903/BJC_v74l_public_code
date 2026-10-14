#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
extension_meta / 05_close_run_log.py

Appends the logs of scripts 02, 03 and 04 into run_log.txt, so that run_log.txt is a single
chronological record of the whole step and not only of script 01.
"""
from __future__ import annotations
import os

import sys
from datetime import datetime
from pathlib import Path

EXT = Path(os.environ.get("SIGLEC15_EXT_META_ROOT", r"E:\Tex\HMlong\BJC_EXT_v1_20260927\extension_meta"))

SECTIONS = [
    ("02_crosscheck_metafor.R (independent engine: R 4.5.1 / metafor 4.8.0)",
     EXT / "provenance" / "crosscheck_metafor_log.txt"),
    ("04_verify_summary.py (write-up vs data)",
     EXT / "provenance" / "verify_summary_log.txt"),
]


def main() -> int:
    parts = [
        "\n" + "#" * 100,
        f"# RUN {datetime.now().isoformat(timespec='seconds')}  05_close_run_log.py "
        f"-- collected logs of scripts 02/03/04",
        "#" * 100,
    ]
    for title, path in SECTIONS:
        parts.append("")
        parts.append(f"### {title}")
        parts.append(f"### source: {path}")
        parts.append("")
        if path.exists():
            parts.append(path.read_text(encoding="utf-8").rstrip())
        else:
            parts.append(f"[MISSING: {path} -- not produced]")
        parts.append("")
    parts.append("### sessionInfo.txt is written separately by 03_session_info.py "
                 f"({EXT / 'sessionInfo.txt'})")
    with open(EXT / "run_log.txt", "a", encoding="utf-8") as fh:
        fh.write("\n".join(parts) + "\n")
    print(f"appended 02/03/04 logs to {EXT / 'run_log.txt'}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
