#!/usr/bin/env python3
"""Install an executable local runtime without changing upstream EviAnn sources."""

import argparse
import hashlib
import json
import os
from pathlib import Path
import re
import shlex
import shutil
import subprocess


DEPENDENCIES = (
    "TransDecoder.Predict", "TransDecoder.LongOrfs", "miniprot",
    "snap", "fathom", "forge", "hmm-assembler.pl",
)


def manifest_entries(source):
    makefile = (source / "Makefile.am").read_text()
    makefile = re.sub(r"\\\r?\n", " ", makefile)
    match = re.search(r"^dist_bin_SCRIPTS\s*=\s*(.*)$", makefile, re.MULTILINE)
    if not match:
        raise ValueError("Makefile.am has no dist_bin_SCRIPTS installation list")
    return shlex.split(match.group(1), comments=True)


def install_runtime(source, outdir):
    source, outdir = Path(source).resolve(), Path(outdir).resolve()
    plan = {}
    for entry in manifest_entries(source):
        original = source / entry
        if entry.startswith("src/"):
            if not original.is_file():
                raise ValueError(f"missing upstream script: {original}")
            plan[original.name] = ("copy", original)
        elif not entry.endswith((".tar.gz", ".tgz", ".zip")):
            if not original.is_file() or not os.access(original, os.X_OK):
                raise ValueError(f"bundled program is missing or not executable: {original}")
            plan[original.name] = ("symlink", original.resolve())
    for name in DEPENDENCIES:
        executable = shutil.which(name)
        if not executable:
            raise ValueError(f"{name} is missing from PATH; activate the existing annotation environment")
        plan[name] = ("symlink", Path(executable).resolve())
    ufasta = shutil.which("ufasta")
    if ufasta:
        plan["ufasta"] = ("symlink", Path(ufasta).resolve())
    elif "ufasta.pl" not in plan:
        fallback = shutil.which("ufasta.pl")
        if not fallback:
            raise ValueError("ufasta or ufasta.pl is required")
        plan["ufasta.pl"] = ("symlink", Path(fallback).resolve())

    # Refuse destinations that would overwrite upstream files or directories.
    for name, (_, original) in plan.items():
        target = outdir / name
        if target == original or (target.is_dir() and not target.is_symlink()):
            raise ValueError(f"refusing to replace source or directory: {target}")
    manifest = outdir / "installation.json"
    if manifest.is_dir():
        raise ValueError(f"refusing to replace directory: {manifest}")
    try:
        commit = subprocess.check_output(
            ["git", "-C", str(source), "rev-parse", "HEAD"], text=True,
            stderr=subprocess.DEVNULL,
        ).strip()
    except (OSError, subprocess.CalledProcessError):
        commit = None
    outdir.mkdir(parents=True, exist_ok=True)
    installed = {}
    for name, (action, original) in plan.items():
        target = outdir / name
        if target.is_symlink() or target.exists():
            target.unlink()
        if action == "copy":
            shutil.copyfile(original, target)
            target.chmod(0o755)
        else:
            target.symlink_to(original)
        installed[name] = {"kind": action, "source": str(original)}
        if action == "copy":
            installed[name]["sha256"] = hashlib.sha256(original.read_bytes()).hexdigest()
    if manifest.is_symlink():
        manifest.unlink()
    manifest.write_text(json.dumps({
        "source": str(source), "source_commit": commit, "tools": installed,
    }, indent=2) + "\n")
    return manifest


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, required=True, help="upstream EviAnn checkout")
    parser.add_argument("--outdir", type=Path, required=True, help="private installation directory")
    args = parser.parse_args()
    try:
        print(install_runtime(args.source, args.outdir))
    except (ValueError, FileNotFoundError) as error:
        parser.error(str(error))


if __name__ == "__main__":
    main()
