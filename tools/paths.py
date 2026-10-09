# -*- coding: utf-8 -*-
"""Shared path resolution for the four-layer analysis code.

The individual-level inputs (the de-identified Chinese screening extract and the
NHANES working set) are **not** distributed with this repository. Every script
that needs them resolves them through this module, so the whole pipeline can be
repointed in one place instead of by editing absolute paths in dozens of files.

Resolution order for the analysis source directory
--------------------------------------------------
1. environment variable ``LDF_SRC``
2. the ``_SRC`` default below

The raw FinnGen release used by the control-exclusion sensitivity analysis is
resolved the same way through ``LDF_EXMORE_GZ`` / the ``_EXMORE_GZ`` default.

Nothing here reads data at import time; missing inputs are reported only when a
script asserts that a specific file is required.
"""

from __future__ import annotations

import os
from pathlib import Path

# --------------------------------------------------------------------------
# defaults (author's local analysis directory; override with env vars)
# --------------------------------------------------------------------------
DEFAULT_SRC = r"E:/workbuddy/2026-09-04-18-19-16"
DEFAULT_EXMORE_GZ = (
    r"E:/workbuddy/BMC Genomics投稿资料/DR，DN，DM芬兰原始数据/"
    r"finngen_R13_DM_RETINOPATHY_EXMORE.gz"
)

ENV_SRC = "LDF_SRC"
ENV_EXMORE_GZ = "LDF_EXMORE_GZ"


def repo_root() -> Path:
    """Repository root (the directory that contains ``output_tables/``)."""
    return Path(__file__).resolve().parent.parent


def out_dir() -> Path:
    d = repo_root() / "output_tables"
    d.mkdir(parents=True, exist_ok=True)
    return d


def src_dir() -> Path:
    """Analysis source directory holding the non-distributed inputs."""
    return Path(os.environ.get(ENV_SRC, DEFAULT_SRC)).expanduser()


def exmore_gz() -> Path:
    return Path(os.environ.get(ENV_EXMORE_GZ, DEFAULT_EXMORE_GZ)).expanduser()


class MissingInput(FileNotFoundError):
    """Raised when a required, non-distributed input is absent.

    The message states exactly which environment variable to set, so a failure
    is actionable rather than a bare ``FileNotFoundError``.
    """


#: logical name -> (sub-path relative to ``src_dir()``, environment variable hint)
_INPUTS = {
    "dev_set": "phase2_dev/dev_set.csv",
    "dev_set_v2": "phase2_dev/dev_set_v2.csv",
    "dev_set_v2_dedup": "phase2_dev/dev_set_v2_dedup.csv",
    "severity_rebuild_audit": "phase2_dev/severity_rebuild_audit.csv",
    "phase1_out": "phase1_out",
    "phase2_ext_nhanes_b": "phase2_ext/NHANES_B.csv",
    "raw_screening": "raw/",
}


def input_path(name: str) -> Path:
    """Resolve a named input. Does not check existence."""
    try:
        rel = _INPUTS[name]
    except KeyError:
        raise KeyError(
            "unknown input %r; known inputs: %s"
            % (name, ", ".join(sorted(_INPUTS)))
        )
    return src_dir() / rel


def require(name: str) -> Path:
    """Resolve a named input and raise ``MissingInput`` if it is absent."""
    p = input_path(name)
    if not p.exists():
        raise MissingInput(
            "required input %r not found at %s\n"
            "These inputs are not distributed with the repository "
            "(individual-level / licence-restricted data).\n"
            "Point the pipeline at a local analysis directory:\n"
            "    export %s=/path/to/analysis/dir      # bash / Git Bash\n"
            "    $env:%s = 'X:\\path\\to\\analysis'    # PowerShell\n"
            "See README.md ('Running the scripts') and INPUTS.md."
            % (name, p, ENV_SRC, ENV_SRC)
        )
    return p


def describe() -> str:
    """One-line human-readable description of the active configuration."""
    return (
        "LDF_SRC=%s (%s) | LDF_EXMORE_GZ=%s (%s)"
        % (
            src_dir(),
            "present" if src_dir().exists() else "ABSENT",
            exmore_gz(),
            "present" if exmore_gz().exists() else "ABSENT",
        )
    )


if __name__ == "__main__":  # pragma: no cover - manual probe
    print("repo_root :", repo_root())
    print("out_dir   :", out_dir())
    for _n in sorted(_INPUTS):
        _p = input_path(_n)
        print("%-22s %-6s %s" % (_n, "OK" if _p.exists() else "MISS", _p))
