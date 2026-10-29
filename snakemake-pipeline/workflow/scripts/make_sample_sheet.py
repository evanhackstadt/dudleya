#!/usr/bin/env python3
"""
make_sample_sheet.py -- build config/samples.tsv from the legacy samples.yaml.

This is the ONLY place the FASTQ naming convention is encoded. Everything
downstream reads columns from the sheet and never parses a sample name.

Design notes
  * Metadata is parsed from the FASTQ *filename* (the primary record), not from
    the legacy YAML key. The key is cross-checked against it and kept as
    `legacy_id`, so existing BAMs/results can be mapped to new ids.
  * Unparseable or inconsistent entries are collected and reported together;
    nothing is written unless every sample parses.
  * Row order follows the legacy YAML, so bam.filelist / population.info order
    is unchanged.
  * `batch` (which pooled sequencing run) comes from config/batches.yaml, which
    maps batch tags to raw-reads folders. It is unrelated to `field_tag`, the
    per-sample collection tag that appears in some filenames after LP_<n>_.

Usage (from the repo root)
  python workflow/scripts/make_sample_sheet.py --check-files
"""

import argparse
import csv
import os
import re
import sys
from collections import Counter, defaultdict

import yaml

from sampleinfo import Taxonomy, check_sheet

# <popcode>_LP_<lp>_<tag>_S<n>_L<lane>_R<1|2>_<chunk>.fastq.gz
#   popcode: 'DAN' or 'CY_SIE' (optional taxon prefix + population)
#   tag:     'Du-590' (main batch) | 'EAMF6' | 'SSRJDu-791' (biol173 batch)
FASTQ_RE = re.compile(
    r"^(?P<popcode>[A-Z0-9]+(?:_[A-Z0-9]+)?)"
    r"_LP_(?P<lp>[A-Za-z0-9-]+)"
    r"_(?P<tag>[A-Za-z0-9-]+)"
    r"_S(?P<snum>\d+)_L(?P<lane>\d{3})"
    r"_R(?P<read>[12])_(?P<chunk>\d{3})\.fastq\.gz$"
)
DU_SUFFIX_RE = re.compile(r"(Du-[0-9A-Za-z]+)$")

COLUMNS = [
    "sample", "legacy_id", "taxon", "species", "population",
    "lp_id", "du_id", "field_tag", "batch", "lane", "r1", "r2",
]


# --- legacy YAML loading ----------------------------------------------------

class UniqueKeyLoader(yaml.SafeLoader):
    """SafeLoader that rejects duplicate mapping keys (PyYAML silently keeps the last)."""


def _construct_unique_mapping(loader, node, deep=False):
    seen = set()
    for key_node, _ in node.value:
        key = loader.construct_object(key_node, deep=deep)
        if key in seen:
            raise yaml.constructor.ConstructorError(
                None, None, f"duplicate key {key!r}", key_node.start_mark
            )
        seen.add(key)
    return yaml.SafeLoader.construct_mapping(loader, node, deep)


UniqueKeyLoader.add_constructor(
    yaml.resolver.BaseResolver.DEFAULT_MAPPING_TAG, _construct_unique_mapping
)


def load_legacy(path):
    with open(path) as fh:
        return yaml.load(fh, Loader=UniqueKeyLoader)


# --- parsing ----------------------------------------------------------------

def parse_fastq(path):
    m = FASTQ_RE.match(os.path.basename(path))
    if not m:
        raise ValueError(f"unparseable FASTQ filename: {os.path.basename(path)}")
    return m.groupdict()


def load_batches(path):
    """Load batches.yaml -> {normalized folder: tag}. Validates tags and uniqueness."""
    with open(path) as fh:
        raw = yaml.load(fh, Loader=UniqueKeyLoader)
    entries = (raw or {}).get("batches")
    if not entries:
        raise ValueError(f"{path} must contain a non-empty 'batches:' mapping (tag: folder)")
    by_dir, errs = {}, []
    for tag, folder in entries.items():
        if not re.fullmatch(r"[A-Za-z0-9_-]+", str(tag)):
            errs.append(f"batch tag {tag!r} may only contain letters, digits, '-' and '_'")
        key = os.path.normpath(str(folder))
        if key in by_dir:
            errs.append(f"folder {key} is registered for both {by_dir[key]!r} and {tag!r}")
        by_dir[key] = str(tag)
    if errs:
        raise ValueError(f"Invalid {path}:\n  " + "\n  ".join(errs))
    return by_dir


def batch_for(path, batches):
    folder = os.path.normpath(os.path.dirname(path))
    try:
        return batches[folder]
    except KeyError:
        raise ValueError(
            f"FASTQ folder {folder} is not registered in batches.yaml "
            f"(registered: {', '.join(sorted(batches))})"
        ) from None


def build_row(legacy_id, entry, taxonomy, batches):
    """Return (row_dict, warnings) for one legacy sample entry; raise ValueError if invalid."""
    if not isinstance(entry, dict) or "r1" not in entry or "r2" not in entry:
        raise ValueError("entry must provide both r1 and r2")
    r1, r2 = entry["r1"], entry["r2"]
    f1, f2 = parse_fastq(r1), parse_fastq(r2)
    if os.path.dirname(r1) != os.path.dirname(r2):
        raise ValueError("r1 and r2 are in different folders")

    if f1["read"] != "1" or f2["read"] != "2":
        raise ValueError("r1/r2 are not R1/R2 respectively (swapped?)")
    for key in ("popcode", "lp", "tag", "snum", "lane"):
        if f1[key] != f2[key]:
            raise ValueError(f"r1/r2 disagree on {key}: {f1[key]!r} vs {f2[key]!r}")

    warnings = []
    if f1["chunk"] != f2["chunk"]:
        warnings.append(f"r1 chunk {f1['chunk']} != r2 chunk {f2['chunk']} (typo?)")

    sample = f"{f1['popcode']}_LP_{f1['lp']}_{f1['tag']}"
    if sample != legacy_id:
        if legacy_id.endswith("_") and sample.startswith(legacy_id):
            warnings.append(f"legacy id {legacy_id!r} repaired to {sample!r}")
        else:
            raise ValueError(f"legacy key disagrees with FASTQ name (-> {sample!r})")

    du = DU_SUFFIX_RE.search(f1["tag"])
    taxon = taxonomy.resolve(f1["popcode"])
    row = {
        "sample": sample,
        "legacy_id": legacy_id,
        "taxon": taxon,
        "species": taxonomy.species_of(taxon),
        "population": f1["popcode"],
        "lp_id": f1["lp"],
        "du_id": du.group(1) if du else "",
        "field_tag": f1["tag"][: du.start()] if du else f1["tag"],
        "batch": batch_for(r1, batches),
        "lane": f"L{f1['lane']}",
        "r1": r1,
        "r2": r2,
    }
    return row, warnings


# --- driver -----------------------------------------------------------------

def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--legacy-yaml", default="config/samples.legacy.yaml")
    ap.add_argument("--taxa", default="config/taxa.yaml")
    ap.add_argument("--batches", default="config/batches.yaml")
    ap.add_argument("--out", default="config/samples.tsv")
    ap.add_argument("--check-files", action="store_true", help="verify every FASTQ exists on disk")
    args = ap.parse_args(argv)

    taxonomy = Taxonomy.load(args.taxa)
    batches = load_batches(args.batches)
    legacy = load_legacy(args.legacy_yaml)

    rows, errors, warnings = [], [], []
    for legacy_id, entry in legacy["samples"].items():
        try:
            row, warns = build_row(legacy_id, entry, taxonomy, batches)
        except (ValueError, TypeError) as exc:
            errors.append(f"{legacy_id}: {exc}")
            continue
        rows.append(row)
        warnings += [f"{row['sample']}: {w}" for w in warns]
        if args.check_files:
            errors += [
                f"{row['sample']}: missing file {row[k]}"
                for k in ("r1", "r2") if not os.path.isfile(row[k])
            ]

    # cross-row checks
    counts = Counter(r["sample"] for r in rows)
    errors += [f"duplicate sample id after parsing: {s}" for s, n in counts.items() if n > 1]
    by_du = defaultdict(list)
    by_tag = defaultdict(list)
    for r in rows:
        if r["du_id"]:
            by_du[r["du_id"]].append(r["sample"])
        if r["field_tag"]:
            by_tag[(r["population"], r["field_tag"])].append(r["sample"])
    warnings += [f"{du} appears in multiple samples: {', '.join(s)}" for du, s in by_du.items() if len(s) > 1]
    warnings += [
        f"field tag {tag!r} reused within population {pop}: {', '.join(s)}"
        for (pop, tag), s in by_tag.items() if len(s) > 1
    ]
    used = {r["batch"] for r in rows}
    warnings += [
        f"batch {tag!r} is registered ({folder}) but has no samples"
        for folder, tag in batches.items() if tag not in used
    ]
    try:
        check_sheet(rows, taxonomy)
    except ValueError as exc:
        errors.append(str(exc))

    for w in warnings:
        print(f"WARN  {w}", file=sys.stderr)
    if errors:
        for e in errors:
            print(f"ERROR {e}", file=sys.stderr)
        print(f"\n{len(errors)} error(s); {args.out} not written.", file=sys.stderr)
        return 1

    with open(args.out, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=COLUMNS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)

    print(f"Wrote {len(rows)} samples to {args.out} ({len(warnings)} warning(s))", file=sys.stderr)
    for title, key in (("taxon", "taxon"), ("batch", "batch"), ("lane", "lane")):
        summary = ", ".join(f"{k}={n}" for k, n in sorted(Counter(r[key] for r in rows).items()))
        print(f"  by {title}: {summary}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())