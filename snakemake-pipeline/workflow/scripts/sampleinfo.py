"""
sampleinfo.py -- the single home for Dudleya sample-metadata semantics.

Anything that needs to know what a taxon code means, which color/marker a
species gets, or how a population code maps to a taxon imports it from here.
No other script should parse sample names or hard-code taxa.

Consumers (Snakefile, make_sample_sheet.py, pcangsd_visualize.py, ...) read
columns from config/samples.tsv and use Taxonomy for labels and styling.
"""

import yaml

RANKS = ("species", "subspecies")


class TaxonomyError(ValueError):
    """Invalid taxa.yaml, or a popcode that cannot be resolved to a taxon."""


class Taxonomy:
    """Taxon tree loaded from taxa.yaml (two ranks: species > subspecies)."""

    def __init__(self, raw):
        if not isinstance(raw, dict) or not raw.get("taxa"):
            raise TaxonomyError("taxa.yaml must contain a non-empty 'taxa:' mapping")
        self.taxa = raw["taxa"]
        self.implicit = raw.get("implicit_taxon")
        self._prefix_to_taxon = {}
        self._validate()

    @classmethod
    def load(cls, path):
        with open(path) as fh:
            return cls(yaml.safe_load(fh))

    # -- validation ---------------------------------------------------------

    def _validate(self):
        errs = []
        for code, t in self.taxa.items():
            for key in ("label", "rank"):
                if key not in t:
                    errs.append(f"{code}: missing '{key}'")
            rank = t.get("rank")
            if rank not in RANKS:
                errs.append(f"{code}: rank must be one of {RANKS}, got {rank!r}")
                continue
            parent = t.get("parent")
            if rank == "species":
                if parent:
                    errs.append(f"{code}: a species cannot have a parent")
                for key in ("color", "marker"):
                    if key not in t:
                        errs.append(f"{code}: species requires '{key}'")
            elif parent not in self.taxa or self.taxa[parent].get("rank") != "species":
                errs.append(f"{code}: subspecies requires 'parent' to be a species code")
            for prefix in t.get("prefixes") or []:
                if prefix in self._prefix_to_taxon:
                    errs.append(
                        f"{code}: prefix {prefix!r} already used by "
                        f"{self._prefix_to_taxon[prefix]}"
                    )
                self._prefix_to_taxon[prefix] = code
        if self.implicit not in self.taxa:
            errs.append(f"implicit_taxon {self.implicit!r} is not a defined taxon")
        if errs:
            raise TaxonomyError("Invalid taxa.yaml:\n  " + "\n  ".join(errs))

    # -- queries ------------------------------------------------------------

    def resolve(self, popcode):
        """Map a popcode ('CY_SIE', 'DAN', 'CYPU_XYZ') to its most specific taxon code."""
        parts = popcode.split("_")
        if len(parts) == 1:
            return self.implicit
        if len(parts) == 2 and parts[0] in self._prefix_to_taxon:
            return self._prefix_to_taxon[parts[0]]
        known = ", ".join(sorted(self._prefix_to_taxon))
        raise TaxonomyError(
            f"Cannot resolve popcode {popcode!r}: prefix {parts[0]!r} is not in "
            f"taxa.yaml (known prefixes: {known})"
        )

    def species_of(self, code):
        """Species-level code for any taxon (itself, or its parent)."""
        return self.taxa[code].get("parent") or code

    def label(self, code):
        return self.taxa[code]["label"]

    def style(self, code):
        """{'color', 'marker'} for a taxon; subspecies inherit from the species."""
        t = self.taxa[code]
        sp = self.taxa[self.species_of(code)]
        return {
            "color": t.get("color", sp["color"]),
            "marker": t.get("marker", sp["marker"]),
        }

    def species_codes(self):
        return [c for c, t in self.taxa.items() if t["rank"] == "species"]


def check_sheet(rows, taxonomy):
    """Cross-row consistency checks that a per-row JSON schema cannot express.

    rows: iterable of dicts (e.g. DataFrame.to_dict('records')).
    Raises ValueError listing every problem found.
    """
    problems, seen, pop_taxon = [], set(), {}
    for r in rows:
        sample = r["sample"]
        if sample in seen:
            problems.append(f"duplicate sample id: {sample}")
        seen.add(sample)

        taxon = r["taxon"]
        if taxon not in taxonomy.taxa:
            problems.append(f"{sample}: taxon {taxon!r} not in taxa.yaml")
            continue
        if r["species"] != taxonomy.species_of(taxon):
            problems.append(
                f"{sample}: species {r['species']!r} disagrees with taxa.yaml "
                f"({taxonomy.species_of(taxon)!r})"
            )
        prev = pop_taxon.setdefault(r["population"], taxon)
        if prev != taxon:
            problems.append(
                f"{sample}: population {r['population']!r} maps to both {prev} and {taxon}"
            )
    if problems:
        raise ValueError("Sample sheet failed consistency checks:\n  " + "\n  ".join(problems))


def select_samples(rows, subset=None, exclude=None):
    """Filter sheet records for a partial run.

    subset:  {column: value-or-list}; a row is kept if it matches ALL columns
             (any listed value within a column), e.g. {"batch": ["su26"], "species": "DUSE"}
    exclude: list of sample ids to drop afterwards.
    Row order is preserved. Typos fail loudly: unknown columns, values that match
    nothing, unknown sample ids, and empty selections all raise ValueError.
    """
    rows = list(rows)
    columns = set(rows[0]) if rows else set()
    wanted, problems = {}, []
    for col, vals in (subset or {}).items():
        if col not in columns:
            problems.append(f"subset column {col!r} is not in the sample sheet")
            continue
        vals = [str(v) for v in ([vals] if isinstance(vals, (str, int)) else vals)]
        present = {r[col] for r in rows}
        missing = [v for v in vals if v not in present]
        if missing:
            problems.append(
                f"subset {col}: no samples match {missing} "
                f"(available: {', '.join(sorted(present)[:15])})"
            )
        wanted[col] = set(vals)
    exclude = [str(s) for s in (exclude or [])]
    unknown = sorted(set(exclude) - {r["sample"] for r in rows})
    if unknown:
        problems.append(f"exclude: unknown sample id(s) {unknown}")
    if problems:
        raise ValueError("Invalid subset selection:\n  " + "\n  ".join(problems))

    keep = [
        r for r in rows
        if all(r[c] in vs for c, vs in wanted.items()) and r["sample"] not in exclude
    ]
    if not keep:
        raise ValueError("Subset selects zero samples")
    return keep