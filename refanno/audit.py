import hashlib
import json
import os
import platform
import subprocess
import sys
from datetime import datetime, timezone
from functools import lru_cache
from importlib.metadata import PackageNotFoundError, version

from . import OUTPUT_SCHEMA_VERSION, __version__


@lru_cache(maxsize=128)
def cached_sha256(path, size, mtime_ns):
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def file_sha256(path):
    stat = os.stat(path)
    return cached_sha256(os.path.abspath(path), stat.st_size, stat.st_mtime_ns)


def read_json(path):
    with open(path) as handle:
        return json.load(handle)


def write_json(path, data):
    temporary = path + ".tmp"
    with open(temporary, "w") as handle:
        json.dump(data, handle, ensure_ascii=False, indent=2)
        handle.write("\n")
    os.replace(temporary, path)


def executable_version(executable):
    try:
        result = subprocess.run([executable, "--version"], text=True, capture_output=True, timeout=15)
        text = (result.stdout or result.stderr).strip()
        return text.splitlines()[0] if text else "unknown"
    except (OSError, subprocess.SubprocessError):
        return "unknown"


def package_version(name):
    try:
        return version(name)
    except PackageNotFoundError:
        return "not installed"


def bam_index(path):
    candidates = [path + ".bai", path + ".csi"]
    if path.lower().endswith(".bam"):
        base = path[:-4]
        candidates.extend([base + ".bai", base + ".csi"])
    return next((candidate for candidate in candidates if os.path.isfile(candidate)), "")


def input_records(args):
    records = []

    def add(role, path, name=""):
        record = {
            "role": role, "name": name, "path": os.path.abspath(path),
            "bytes": os.path.getsize(path), "sha256": file_sha256(path),
        }
        if role == "short_read_bam":
            index = bam_index(path)
            record["index_path"] = os.path.abspath(index) if index else ""
            record["index_sha256"] = file_sha256(index) if index else ""
        records.append(record)

    for spec in args.set:
        name, path = spec.split("=", 1)
        add("annotation", path, name)
    add("genome", args.genome)
    for path in args.bam:
        add("short_read_bam", path)
    for path in args.isoseq:
        add("isoseq", path)
    for path in args.sj:
        add("junction_table", path)
    if args.proteins:
        add("proteins", args.proteins)
    if args.miniprot_gff:
        add("protein_alignment", args.miniprot_gff)
    return records


def write_run_manifest(args, argv, started, summary):
    command = [sys.argv[0]] + sys.argv[1:] if argv is None else ["refanno"] + list(argv)
    ended = datetime.now(timezone.utc)
    dependencies = {"python": platform.python_version(), "pysam": package_version("pysam")}
    if args.proteins and not args.miniprot_gff:
        dependencies["miniprot"] = executable_version(args.miniprot)
    outputs = []
    for name in sorted(os.listdir(args.outdir)):
        path = os.path.join(args.outdir, name)
        if name == "run_manifest.json" or not os.path.isfile(path) or name.endswith(".tmp"):
            continue
        outputs.append({"path": name, "bytes": os.path.getsize(path), "sha256": file_sha256(path)})
    manifest = {
        "refanno_version": __version__, "output_schema_version": OUTPUT_SCHEMA_VERSION,
        "status": "completed", "command": command, "parameters": vars(args),
        "started_utc": started.isoformat(), "ended_utc": ended.isoformat(),
        "duration_seconds": (ended - started).total_seconds(),
        "dependencies": dependencies, "inputs": input_records(args), "outputs": outputs,
        "summary": summary,
    }
    write_json(os.path.join(args.outdir, "run_manifest.json"), manifest)
