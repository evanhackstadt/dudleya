# helper-scripts/rename_files_with_popcode.py

'''
Usage:

python helper-scripts/rename_files_with_popcode.py /path/to/fastq-files
python helper-scripts/rename_files_with_popcode.py /path/to/fastq-files --tag su26
python helper-scripts/rename_files_with_popcode.py /path/to/fastq-files --dictionary /path/to/du_pop_dictionary.json

The script reads `du_pop_dictionary.json` beside itself, 
validates every candidate before changing files, 
and restores original filenames if a later rename fails.
'''


import sys
import argparse
import json
import re
import uuid
from pathlib import Path


def load_dictionary(path):
    with path.open() as file:
        entries = json.load(file)

    if not isinstance(entries, list):
        raise ValueError("Dictionary must be a JSON array of objects.")

    dictionary = {}
    for entry in entries:
        if not isinstance(entry, dict):
            raise ValueError(f"Invalid dictionary entry: {entry!r}")

        du = entry.get("DU#")
        code = entry.get("Code")
        if du is None or code is None or not str(code).strip():
            raise ValueError(f"Missing or empty DU# / Code in entry: {entry!r}")

        du = str(du).strip()
        if du in dictionary:
            raise ValueError(f"Duplicate DU# in dictionary: {du}")

        dictionary[du] = str(code).strip()

    return dictionary


def is_fastq(path):
    name = path.name.lower()
    return name.endswith((".fastq", ".fastq.gz", ".fq", ".fq.gz"))


def get_lp_filename(filename):
    match = re.search(r"(?:^|_)LP_", filename)
    if not match:
        return None

    start = match.start()
    if filename[start] == "_":
        start += 1
    return filename[start:]


def plan_renames(directory, dictionary, tag):
    planned = []
    for path in sorted(directory.iterdir()):
        if not path.is_file() or not is_fastq(path):
            continue

        lp_filename = get_lp_filename(path.name)
        if lp_filename is None:
            print(f"Skipping (no LP_ marker): {path.name}")
            continue

        match = re.search(r"(?i)(?:^|_)Du-(?P<du>[^_./]+)", lp_filename)
        if not match:
            raise ValueError(f"Could not extract DU# from {path.name}")

        du = match.group("du")
        if du not in dictionary:
            raise ValueError(f"No population code found for DU#{du} in {path.name}")

        prefix = dictionary[du]
        if tag:
            prefix = f"{prefix}_{tag}"
        target = path.with_name(f"{prefix}_{lp_filename}")

        if target != path:
            planned.append((path, target))

    sources = {source for source, _ in planned}
    targets = set()
    for source, target in planned:
        if target in targets:
            raise ValueError(f"Multiple files would produce {target.name}")
        targets.add(target)

        if target.exists() and target not in sources:
            raise ValueError(f"Target already exists and will not be renamed: {target.name}")

    return planned


def rename_files(planned):
    # Stage all files first, allowing safe swaps and preventing partial collisions.
    staged = []
    committed = []

    try:
        for source, target in planned:
            temp = source.with_name(f".du-rename-{uuid.uuid4().hex}.tmp")
            source.rename(temp)
            staged.append((source, temp, target))

        for source, temp, target in staged:
            temp.rename(target)
            committed.append((source, temp, target))

    except Exception:
        # Return committed targets to temporary paths, then restore original names.
        for source, temp, target in reversed(committed):
            if target.exists():
                target.rename(temp)

        for source, temp, target in reversed(staged):
            if temp.exists():
                temp.rename(source)

        raise


def main():
    parser = argparse.ArgumentParser(
        description="Set population-code and optional tag prefixes on FASTQ filenames."
    )
    parser.add_argument("directory", type=Path, help="Directory containing FASTQ files")
    parser.add_argument(
        "--dictionary",
        type=Path,
        default=Path(__file__).with_name("du_pop_dictionary.json"),
        help="Path to the population dictionary JSON",
    )
    parser.add_argument(
        "--tag",
        help="Optional tag added after the population code, e.g. su26",
    )
    args = parser.parse_args()

    if args.tag is not None and (
        not args.tag.strip()
        or "_" in args.tag
        or "/" in args.tag
        or "\\" in args.tag
    ):
        parser.error("--tag must be non-empty and contain no underscores or path separators")

    directory = args.directory.expanduser().resolve()
    dictionary_path = args.dictionary.expanduser().resolve()

    if not directory.is_dir():
        raise SystemExit(f"Directory does not exist: {directory}")
    if not dictionary_path.is_file():
        raise SystemExit(f"Dictionary not found: {dictionary_path}")

    dictionary = load_dictionary(dictionary_path)
    planned = plan_renames(directory, dictionary, args.tag)

    if not planned:
        print("No filenames need changing.")
        return

    rename_files(planned)
    for source, target in planned:
        print(f"{source.name} -> {target.name}")


if __name__ == "__main__":
    try:
        main()
    except Exception as error:
        print(f"Error: {error}", file=sys.stderr)
        sys.exit(1)