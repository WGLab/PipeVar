#!/usr/bin/env python3
"""Compare inheritance-aware PipeVar and Exomiser gene rankings."""

from __future__ import annotations

import argparse
import csv
import re
import string
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Sequence


CUTOFFS = (1, 5, 10, 20)
VALID_STATUSES = {"ok", "empty"}
INHERITANCE_ALIASES = {
    "AD": "AD",
    "AUTOSOMAL DOMINANT": "AD",
    "DOMINANT": "AD",
    "AR": "AR",
    "AUTOSOMAL RECESSIVE": "AR",
    "RECESSIVE": "AR",
    "XLD": "XLD",
    "XD": "XLD",
    "X LINKED DOMINANT": "XLD",
    "XLR": "XLR",
    "XR": "XLR",
    "X RECESSIVE": "XLR",
    "X LINKED RECESSIVE": "XLR",
    "XL": "XL",
    "X LINKED": "XL",
    "UNKNOWN": "UNSPECIFIED",
    "UNSPECIFIED": "UNSPECIFIED",
    "NOT SPECIFIED": "UNSPECIFIED",
    "NOT KNOWN": "UNSPECIFIED",
    "NA": "UNSPECIFIED",
    "N A": "UNSPECIFIED",
    "NONE": "UNSPECIFIED",
    ".": "UNSPECIFIED",
}


class ComparisonError(ValueError):
    """Raised for invalid truth, result, or output contracts."""


class ResultFormatError(ValueError):
    """Raised when a ranking file is malformed."""


@dataclass(frozen=True)
class TruthCandidate:
    gene: str
    inheritance: str


@dataclass(frozen=True)
class TruthSample:
    sample_id: str
    candidates: tuple[TruthCandidate, ...]


@dataclass(frozen=True)
class InheritanceRanking:
    gene: str
    inheritance: str
    rank: int


@dataclass(frozen=True)
class ParsedResult:
    status: str
    rankings: tuple[InheritanceRanking, ...]
    error: str = ""


@dataclass(frozen=True)
class ToolMatch:
    exact_gene: str = ""
    exact_inheritance: str = ""
    exact_rank: int | None = None
    gene_only_gene: str = ""
    gene_only_rank: int | None = None
    gene_only_predicted_models: tuple[str, ...] = ()


@dataclass(frozen=True)
class SampleResult:
    truth: TruthSample
    pipevar_path: Path
    exomiser_path: Path
    pipevar_result: ParsedResult
    exomiser_result: ParsedResult
    pipevar_match: ToolMatch
    exomiser_match: ToolMatch

    @property
    def paired(self) -> bool:
        return (
            self.pipevar_result.status in VALID_STATUSES
            and self.exomiser_result.status in VALID_STATUSES
        )

    @property
    def inheritance_evaluable(self) -> bool:
        return any(
            candidate.inheritance != "UNSPECIFIED"
            for candidate in self.truth.candidates
        )


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Compare whether PipeVar and Exomiser rank the correct gene under "
            "the correct inheritance model at top 1, 5, 10, and 20."
        )
    )
    parser.add_argument("--truth-tsv", required=True, type=Path)
    parser.add_argument("--pipevar-dir", required=True, type=Path)
    parser.add_argument("--exomiser-dir", required=True, type=Path)
    parser.add_argument("--output-prefix", required=True, type=Path)
    parser.add_argument(
        "--pipevar-template",
        default="{sample}.norm.alt.prio_gene.vcf",
        help="Path template relative to --pipevar-dir (default: %(default)s)",
    )
    parser.add_argument(
        "--exomiser-template",
        default="{sample}.exomiser.input.vcf.gz-results.genes.tsv",
        help="Path template relative to --exomiser-dir (default: %(default)s)",
    )
    return parser.parse_args(argv)


def normalize_gene(value: str) -> str:
    return value.strip().upper()


def normalize_inheritance(value: str, context: str) -> str:
    key = re.sub(r"[_\-/]+", " ", value.strip().upper())
    key = re.sub(r"\s+", " ", key)
    if key in INHERITANCE_ALIASES:
        return INHERITANCE_ALIASES[key]
    accepted = "AD, AR, XLD, XLR, XL, or UNSPECIFIED"
    raise ComparisonError(
        f"{context} has unsupported inheritance '{value.strip()}'; expected {accepted}"
    )


def parse_genes(value: str, context: str) -> tuple[str, ...]:
    genes: list[str] = []
    seen: set[str] = set()
    for raw_gene in re.split(r"[;,]", value):
        gene = normalize_gene(raw_gene)
        if gene and gene not in seen:
            genes.append(gene)
            seen.add(gene)
    if not genes:
        raise ComparisonError(f"{context} has no gene symbols")
    return tuple(genes)


def read_truth_tsv(path: Path) -> list[TruthSample]:
    if not path.is_file():
        raise ComparisonError(f"truth TSV not found: {path}")
    try:
        with path.open("r", encoding="utf-8-sig", newline="") as handle:
            rows = [
                (line_number, row)
                for line_number, row in enumerate(
                    csv.reader(handle, delimiter="\t", strict=True), 1
                )
                if any(row)
            ]
    except (OSError, UnicodeError, csv.Error) as exc:
        raise ComparisonError(f"could not parse truth TSV '{path}': {exc}") from exc
    if not rows:
        raise ComparisonError(f"truth TSV is empty: {path}")
    header = [value.strip() for value in rows[0][1]]
    if len(header) != len(set(header)):
        raise ComparisonError("truth TSV contains duplicate column names")
    required = ("sample_id", "genes", "inheritance")
    missing = [column for column in required if column not in header]
    if missing:
        raise ComparisonError(
            "truth TSV is missing required column(s): " + ", ".join(missing)
        )
    indices = {column: header.index(column) for column in required}

    sample_candidates: dict[str, list[TruthCandidate]] = {}
    pair_rows: dict[str, dict[TruthCandidate, int]] = {}
    for line_number, row in rows[1:]:
        if len(row) != len(header):
            raise ComparisonError(
                f"truth TSV row {line_number} has {len(row)} columns; "
                f"expected {len(header)}"
            )
        sample_id = row[indices["sample_id"]].strip()
        if not sample_id:
            raise ComparisonError(f"truth TSV row {line_number} has an empty sample_id")
        genes = parse_genes(
            row[indices["genes"]], f"truth sample '{sample_id}' at row {line_number}"
        )
        inheritance = normalize_inheritance(
            row[indices["inheritance"]],
            f"truth sample '{sample_id}' at row {line_number}",
        )
        candidates = sample_candidates.setdefault(sample_id, [])
        seen_pairs = pair_rows.setdefault(sample_id, {})
        for gene in genes:
            candidate = TruthCandidate(gene, inheritance)
            if candidate in seen_pairs:
                raise ComparisonError(
                    f"truth TSV rows {seen_pairs[candidate]} and {line_number} "
                    f"duplicate {gene}/{inheritance} for sample_id '{sample_id}'"
                )
            seen_pairs[candidate] = line_number
            candidates.append(candidate)
    if not sample_candidates:
        raise ComparisonError("truth TSV has no sample rows")
    return [
        TruthSample(sample_id, tuple(candidates))
        for sample_id, candidates in sample_candidates.items()
    ]


def parse_positive_rank(value: str, context: str) -> int:
    try:
        rank = int(value.strip())
    except ValueError as exc:
        raise ResultFormatError(
            f"{context} has non-integer rank '{value.strip()}'"
        ) from exc
    if rank < 1:
        raise ResultFormatError(f"{context} has non-positive rank {rank}")
    return rank


def canonical_prediction(value: str) -> str | None:
    key = re.sub(r"[_\-/]+", " ", value.strip().upper())
    key = re.sub(r"\s+", " ", key)
    aliases = {
        "AD": "AD",
        "AUTOSOMAL DOMINANT": "AD",
        "DOMINANT": "AD",
        "AR": "AR",
        "AUTOSOMAL RECESSIVE": "AR",
        "RECESSIVE": "AR",
        "XLD": "XLD",
        "XD": "XLD",
        "X DOMINANT": "XLD",
        "X LINKED DOMINANT": "XLD",
        "XLR": "XLR",
        "XR": "XLR",
        "X RECESSIVE": "XLR",
        "X LINKED RECESSIVE": "XLR",
    }
    return aliases.get(key)


def inheritance_matches(truth: str, predicted: str) -> bool:
    if truth == "UNSPECIFIED":
        return False
    if truth == "XL":
        return predicted in {"XLD", "XLR"}
    return truth == predicted


def parse_info_models(vcf_tail: str) -> tuple[str, ...]:
    models: list[str] = []
    for record in vcf_tail.split(" ||| "):
        fields = record.split("\t")
        if len(fields) < 8:
            continue
        info = fields[7]
        for item in info.split(";"):
            if item.startswith("MODEL="):
                model_value = item.split("=", 1)[1]
                model_value = re.sub(
                    r"\\x(?:3b|2c|7c)|%(?:3b|2c|7c)",
                    ",",
                    model_value,
                    flags=re.IGNORECASE,
                )
                raw_models = re.split(r"[,|]+", model_value)
            else:
                # Recover legacy malformed MODEL=Dominant;Recessive values,
                # where the second model is parsed as a bare INFO flag.
                raw_models = [item] if canonical_prediction(item) else []
            for raw_model in raw_models:
                model = canonical_prediction(raw_model)
                if model and model not in models:
                    models.append(model)
    return tuple(models)


def has_standalone_duplication_model(vcf_tail: str) -> bool:
    for record in vcf_tail.split(" ||| "):
        fields = record.split("\t")
        if len(fields) < 8:
            continue
        for item in fields[7].split(";"):
            if item.startswith("MODEL="):
                values = re.split(r"[,|]+", item.split("=", 1)[1])
                if any(value.strip().upper() == "DUPLICATION" for value in values):
                    return True
    return False


def infer_pipevar_tier_model(tier: str) -> str | None:
    normalized = tier.strip().upper()
    if "HEMIZYGOUS" in normalized:
        return "XLR"
    if normalized.endswith("_AD"):
        return "AD"
    if "_AR" in normalized:
        return "AR"
    return None


def read_pipevar(path: Path) -> tuple[InheritanceRanking, ...]:
    try:
        lines = path.read_text(encoding="utf-8-sig").splitlines()
    except (OSError, UnicodeError) as exc:
        raise ResultFormatError(f"could not read file: {exc}") from exc
    nonblank = [(number, line) for number, line in enumerate(lines, 1) if line.strip()]
    if not nonblank:
        raise ResultFormatError("file is empty and has no header")
    header = [value.strip() for value in nonblank[0][1].split("\t")]
    required = ("#RANK", "GENE", "PRIORITY_TIER")
    missing = [column for column in required if column not in header]
    if missing:
        raise ResultFormatError("missing required column(s): " + ", ".join(missing))
    indices = {column: header.index(column) for column in required}
    vcf_index = header.index("VCF_LINE") if "VCF_LINE" in header else None
    rankings: list[InheritanceRanking] = []
    seen: set[tuple[str, str, int]] = set()
    for line_number, line in nonblank[1:]:
        fields = line.split("\t", maxsplit=vcf_index if vcf_index is not None else len(header) - 1)
        if len(fields) <= max(indices.values()):
            raise ResultFormatError(f"row {line_number} has too few columns")
        rank = parse_positive_rank(fields[indices["#RANK"]], f"row {line_number}")
        gene = normalize_gene(fields[indices["GENE"]])
        if not gene:
            raise ResultFormatError(f"row {line_number} has an empty GENE")
        vcf_tail = fields[vcf_index] if vcf_index is not None and len(fields) > vcf_index else ""
        models = parse_info_models(vcf_tail) if vcf_tail else ()
        if not models:
            fallback = infer_pipevar_tier_model(fields[indices["PRIORITY_TIER"]])
            models = (fallback,) if fallback else ()
        if not models:
            if has_standalone_duplication_model(vcf_tail):
                continue
            raise ResultFormatError(
                f"row {line_number} has no recognized inheritance prediction"
            )
        for model in models:
            key = (gene, model, rank)
            if key not in seen:
                rankings.append(InheritanceRanking(*key))
                seen.add(key)
    rankings.sort(key=lambda item: item.rank)
    return tuple(rankings)


def read_exomiser(path: Path) -> tuple[InheritanceRanking, ...]:
    try:
        with path.open("r", encoding="utf-8-sig", newline="") as handle:
            rows = [
                (number, row)
                for number, row in enumerate(
                    csv.reader(handle, delimiter="\t", strict=True), 1
                )
                if any(row)
            ]
    except (OSError, UnicodeError, csv.Error) as exc:
        raise ResultFormatError(f"could not parse file: {exc}") from exc
    if not rows:
        raise ResultFormatError("file is empty and has no header")
    header = [value.strip() for value in rows[0][1]]
    required = ("#RANK", "GENE_SYMBOL", "MOI")
    missing = [column for column in required if column not in header]
    if missing:
        raise ResultFormatError("missing required column(s): " + ", ".join(missing))
    indices = {column: header.index(column) for column in required}
    rankings: list[InheritanceRanking] = []
    seen_pairs: set[tuple[str, str]] = set()
    for line_number, row in rows[1:]:
        if len(row) <= max(indices.values()):
            raise ResultFormatError(f"row {line_number} has too few columns")
        rank = parse_positive_rank(row[indices["#RANK"]], f"row {line_number}")
        gene = normalize_gene(row[indices["GENE_SYMBOL"]])
        model = canonical_prediction(row[indices["MOI"]])
        if not gene:
            raise ResultFormatError(f"row {line_number} has an empty GENE_SYMBOL")
        if not model:
            raise ResultFormatError(
                f"row {line_number} has unsupported MOI '{row[indices['MOI']].strip()}'"
            )
        pair = (gene, model)
        if pair in seen_pairs:
            continue
        seen_pairs.add(pair)
        rankings.append(InheritanceRanking(gene, model, rank))
    rankings.sort(key=lambda item: item.rank)
    return tuple(rankings)


def safe_read(path: Path, parser) -> ParsedResult:
    if not path.exists():
        return ParsedResult("missing", (), "file not found")
    if not path.is_file():
        return ParsedResult("malformed", (), "path is not a regular file")
    try:
        rankings = parser(path)
    except ResultFormatError as exc:
        return ParsedResult("malformed", (), str(exc))
    return ParsedResult("ok" if rankings else "empty", rankings)


def find_match(truth: TruthSample, rankings: Iterable[InheritanceRanking]) -> ToolMatch:
    rankings = tuple(rankings)
    truth_genes = {candidate.gene for candidate in truth.candidates}
    exact = [
        ranking
        for ranking in rankings
        for candidate in truth.candidates
        if ranking.gene == candidate.gene
        and inheritance_matches(candidate.inheritance, ranking.inheritance)
    ]
    gene_matches = [ranking for ranking in rankings if ranking.gene in truth_genes]
    best_exact = min(exact, key=lambda item: (item.rank, item.gene)) if exact else None
    best_gene = min(gene_matches, key=lambda item: (item.rank, item.gene)) if gene_matches else None
    predicted_models: tuple[str, ...] = ()
    if best_gene:
        predicted_models = tuple(
            dict.fromkeys(
                ranking.inheritance
                for ranking in gene_matches
                if ranking.gene == best_gene.gene
            )
        )
    return ToolMatch(
        exact_gene=best_exact.gene if best_exact else "",
        exact_inheritance=best_exact.inheritance if best_exact else "",
        exact_rank=best_exact.rank if best_exact else None,
        gene_only_gene=best_gene.gene if best_gene else "",
        gene_only_rank=best_gene.rank if best_gene else None,
        gene_only_predicted_models=predicted_models,
    )


def validate_template(template: str, option: str) -> None:
    fields: list[str] = []
    try:
        for _, field, spec, conversion in string.Formatter().parse(template):
            if field is None:
                continue
            fields.append(field)
            if spec or conversion:
                raise ComparisonError(f"{option} does not support formatting modifiers")
    except ValueError as exc:
        raise ComparisonError(f"invalid {option}: {exc}") from exc
    if fields != ["sample"]:
        raise ComparisonError(f"{option} must contain exactly one '{{sample}}' placeholder")


def resolve_path(root: Path, template: str, sample_id: str) -> Path:
    root = root.resolve()
    path = (root / template.format(sample=sample_id)).resolve()
    try:
        path.relative_to(root)
    except ValueError as exc:
        raise ComparisonError(f"sample '{sample_id}' resolves outside {root}") from exc
    return path


def evaluate(
    samples: Iterable[TruthSample],
    pipevar_dir: Path,
    exomiser_dir: Path,
    pipevar_template: str,
    exomiser_template: str,
) -> list[SampleResult]:
    results: list[SampleResult] = []
    for truth in samples:
        pipevar_path = resolve_path(pipevar_dir, pipevar_template, truth.sample_id)
        exomiser_path = resolve_path(exomiser_dir, exomiser_template, truth.sample_id)
        pipevar_result = safe_read(pipevar_path, read_pipevar)
        exomiser_result = safe_read(exomiser_path, read_exomiser)
        results.append(
            SampleResult(
                truth,
                pipevar_path,
                exomiser_path,
                pipevar_result,
                exomiser_result,
                find_match(truth, pipevar_result.rankings),
                find_match(truth, exomiser_result.rankings),
            )
        )
    return results


def hit(match: ToolMatch, cutoff: int, exact: bool = True) -> bool:
    rank = match.exact_rank if exact else match.gene_only_rank
    return rank is not None and rank <= cutoff


def optional(value: object | None) -> object:
    return "" if value is None else value


def write_details(path: Path, results: Iterable[SampleResult]) -> None:
    fixed = [
        "sample_id", "truth_gene_inheritance", "inheritance_evaluable",
        "pipevar_path", "pipevar_status", "pipevar_error",
        "pipevar_exact_gene", "pipevar_exact_inheritance", "pipevar_exact_rank",
        "pipevar_gene_only_gene", "pipevar_gene_only_rank", "pipevar_gene_predicted_models",
        "exomiser_path", "exomiser_status", "exomiser_error",
        "exomiser_exact_gene", "exomiser_exact_inheritance", "exomiser_exact_rank",
        "exomiser_gene_only_gene", "exomiser_gene_only_rank", "exomiser_gene_predicted_models",
        "included_in_paired_inheritance_summary",
    ]
    cutoff_fields = [
        f"{tool}_{kind}_top{cutoff}"
        for tool in ("pipevar", "exomiser")
        for kind in ("exact", "gene_only")
        for cutoff in CUTOFFS
    ]
    fieldnames = fixed + cutoff_fields
    try:
        with path.open("x", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
            writer.writeheader()
            for result in results:
                row: dict[str, object] = {
                    "sample_id": result.truth.sample_id,
                    "truth_gene_inheritance": ";".join(
                        f"{candidate.gene}:{candidate.inheritance}"
                        for candidate in result.truth.candidates
                    ),
                    "inheritance_evaluable": int(result.inheritance_evaluable),
                    "pipevar_path": result.pipevar_path,
                    "pipevar_status": result.pipevar_result.status,
                    "pipevar_error": result.pipevar_result.error,
                    "pipevar_exact_gene": result.pipevar_match.exact_gene,
                    "pipevar_exact_inheritance": result.pipevar_match.exact_inheritance,
                    "pipevar_exact_rank": optional(result.pipevar_match.exact_rank),
                    "pipevar_gene_only_gene": result.pipevar_match.gene_only_gene,
                    "pipevar_gene_only_rank": optional(result.pipevar_match.gene_only_rank),
                    "pipevar_gene_predicted_models": ";".join(result.pipevar_match.gene_only_predicted_models),
                    "exomiser_path": result.exomiser_path,
                    "exomiser_status": result.exomiser_result.status,
                    "exomiser_error": result.exomiser_result.error,
                    "exomiser_exact_gene": result.exomiser_match.exact_gene,
                    "exomiser_exact_inheritance": result.exomiser_match.exact_inheritance,
                    "exomiser_exact_rank": optional(result.exomiser_match.exact_rank),
                    "exomiser_gene_only_gene": result.exomiser_match.gene_only_gene,
                    "exomiser_gene_only_rank": optional(result.exomiser_match.gene_only_rank),
                    "exomiser_gene_predicted_models": ";".join(result.exomiser_match.gene_only_predicted_models),
                    "included_in_paired_inheritance_summary": int(result.paired and result.inheritance_evaluable),
                }
                for tool, parsed, match in (
                    ("pipevar", result.pipevar_result, result.pipevar_match),
                    ("exomiser", result.exomiser_result, result.exomiser_match),
                ):
                    valid = parsed.status in VALID_STATUSES
                    for cutoff in CUTOFFS:
                        row[f"{tool}_exact_top{cutoff}"] = int(hit(match, cutoff)) if valid and result.inheritance_evaluable else ""
                        row[f"{tool}_gene_only_top{cutoff}"] = int(hit(match, cutoff, exact=False)) if valid else ""
                writer.writerow(row)
    except OSError as exc:
        raise ComparisonError(f"could not write detail report '{path}': {exc}") from exc


def build_summary(results: Iterable[SampleResult]) -> list[dict[str, object]]:
    paired = [result for result in results if result.paired and result.inheritance_evaluable]
    denominator = len(paired)
    rows: list[dict[str, object]] = []
    for cutoff in CUTOFFS:
        pipevar = [hit(result.pipevar_match, cutoff) for result in paired]
        exomiser = [hit(result.exomiser_match, cutoff) for result in paired]
        pipevar_gene = [hit(result.pipevar_match, cutoff, exact=False) for result in paired]
        exomiser_gene = [hit(result.exomiser_match, cutoff, exact=False) for result in paired]
        p_hits, e_hits = sum(pipevar), sum(exomiser)
        rows.append({
            "cutoff": cutoff,
            "paired_inheritance_evaluable_samples": denominator,
            "pipevar_correct_gene_and_inheritance_hits": p_hits,
            "pipevar_correct_gene_and_inheritance_rate": f"{p_hits / denominator:.6f}" if denominator else "NA",
            "pipevar_gene_only_hits": sum(pipevar_gene),
            "exomiser_correct_gene_and_inheritance_hits": e_hits,
            "exomiser_correct_gene_and_inheritance_rate": f"{e_hits / denominator:.6f}" if denominator else "NA",
            "exomiser_gene_only_hits": sum(exomiser_gene),
            "hit_difference_pipevar_minus_exomiser": p_hits - e_hits,
            "rate_difference_percentage_points": f"{(p_hits - e_hits) * 100 / denominator:.6f}" if denominator else "NA",
            "both_correct": sum(p and e for p, e in zip(pipevar, exomiser)),
            "pipevar_only_correct": sum(p and not e for p, e in zip(pipevar, exomiser)),
            "exomiser_only_correct": sum(e and not p for p, e in zip(pipevar, exomiser)),
            "neither_correct": sum(not p and not e for p, e in zip(pipevar, exomiser)),
        })
    return rows


def write_summary(path: Path, rows: list[dict[str, object]]) -> None:
    try:
        with path.open("x", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
            writer.writeheader()
            writer.writerows(rows)
    except OSError as exc:
        raise ComparisonError(f"could not write summary report '{path}': {exc}") from exc


def run(args: argparse.Namespace) -> int:
    for root, label in ((args.pipevar_dir, "PipeVar"), (args.exomiser_dir, "Exomiser")):
        if not root.is_dir():
            raise ComparisonError(f"{label} result directory not found: {root}")
    validate_template(args.pipevar_template, "--pipevar-template")
    validate_template(args.exomiser_template, "--exomiser-template")
    detail_path = Path(f"{args.output_prefix}.inheritance.details.tsv")
    summary_path = Path(f"{args.output_prefix}.inheritance.summary.tsv")
    if not detail_path.parent.is_dir():
        raise ComparisonError(f"output directory not found: {detail_path.parent}")
    existing = [str(path) for path in (detail_path, summary_path) if path.exists()]
    if existing:
        raise ComparisonError("output file(s) already exist: " + ", ".join(existing))
    samples = read_truth_tsv(args.truth_tsv)
    results = evaluate(samples, args.pipevar_dir, args.exomiser_dir, args.pipevar_template, args.exomiser_template)
    summary = build_summary(results)
    write_details(detail_path, results)
    write_summary(summary_path, summary)
    for result in results:
        for tool, parsed, path in (
            ("PipeVar", result.pipevar_result, result.pipevar_path),
            ("Exomiser", result.exomiser_result, result.exomiser_path),
        ):
            if parsed.status not in VALID_STATUSES:
                print(f"WARNING: {result.truth.sample_id} {tool} {parsed.status}: {path}: {parsed.error}", file=sys.stderr)
    denominator = summary[0]["paired_inheritance_evaluable_samples"]
    print(f"Paired inheritance-evaluable samples: {denominator}/{len(results)}")
    for row in summary:
        print(
            f"Top {row['cutoff']}: PipeVar correct gene+inheritance "
            f"{row['pipevar_correct_gene_and_inheritance_hits']}/{denominator}; "
            f"Exomiser {row['exomiser_correct_gene_and_inheritance_hits']}/{denominator}"
        )
    print(f"Detail report: {detail_path}")
    print(f"Summary report: {summary_path}")
    return 0 if denominator else 1


def main(argv: Sequence[str] | None = None) -> int:
    try:
        return run(parse_args(argv))
    except ComparisonError as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
