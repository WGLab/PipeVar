#!/usr/bin/env python3
"""Compare PipeVar, RankVar, and Exomiser rankings against causal genes."""

from __future__ import annotations

import argparse
import csv
import re
import string
import sys
from dataclasses import dataclass
from decimal import Decimal, InvalidOperation
from pathlib import Path
from typing import Callable, Iterable, Sequence


CUTOFFS = (1, 5, 10, 20)
VALID_RESULT_STATUSES = {"ok", "empty"}


class BenchmarkError(ValueError):
    """Raised when benchmark inputs or outputs are invalid."""


class RankingError(ValueError):
    """Raised when a ranking file does not satisfy its format contract."""


@dataclass(frozen=True)
class TruthSample:
    sample_id: str
    genes: tuple[str, ...]


@dataclass(frozen=True)
class RankedGene:
    gene: str
    rank: int
    raw_rank: int


@dataclass(frozen=True)
class RankingResult:
    status: str
    rankings: tuple[RankedGene, ...]
    error: str = ""


@dataclass(frozen=True)
class MatchResult:
    matched_gene: str = ""
    rank: int | None = None
    raw_rank: int | None = None


@dataclass(frozen=True)
class SampleResult:
    truth: TruthSample
    pipevar_path: Path
    rankvar_path: Path
    exomiser_path: Path
    pipevar_result: RankingResult
    rankvar_result: RankingResult
    exomiser_result: RankingResult
    pipevar_match: MatchResult
    rankvar_match: MatchResult
    exomiser_match: MatchResult

    @property
    def paired(self) -> bool:
        return (
            self.pipevar_result.status in VALID_RESULT_STATUSES
            and self.rankvar_result.status in VALID_RESULT_STATUSES
            and self.exomiser_result.status in VALID_RESULT_STATUSES
        )


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Compare known causal genes with PipeVar, RankVar, and Exomiser "
            "gene "
            "rankings at top 1, top 5, top 10, and top 20."
        )
    )
    parser.add_argument("--truth-tsv", required=True, type=Path)
    parser.add_argument("--pipevar-dir", required=True, type=Path)
    parser.add_argument(
        "--rankvar-dir",
        type=Path,
        help="RankVar result directory (default: --pipevar-dir)",
    )
    parser.add_argument("--exomiser-dir", required=True, type=Path)
    parser.add_argument("--output-prefix", required=True, type=Path)
    parser.add_argument(
        "--pipevar-template",
        default="{sample}.norm.alt.prio_gene.vcf",
        help="Path template relative to --pipevar-dir (default: %(default)s)",
    )
    parser.add_argument(
        "--rankvar-template",
        default="{sample}.norm.alt.rank_var.tsv",
        help="Path template relative to --rankvar-dir (default: %(default)s)",
    )
    parser.add_argument(
        "--exomiser-template",
        default="{sample}.exomiser.input.vcf.gz-results.genes.tsv",
        help="Path template relative to --exomiser-dir (default: %(default)s)",
    )
    return parser.parse_args(argv)


def normalize_gene(raw_gene: str) -> str:
    return raw_gene.strip().upper()


def parse_truth_genes(raw_genes: str, context: str) -> tuple[str, ...]:
    genes: list[str] = []
    seen: set[str] = set()
    for value in re.split(r"[;,]", raw_genes):
        gene = normalize_gene(value)
        if not gene or gene in seen:
            continue
        genes.append(gene)
        seen.add(gene)
    if not genes:
        raise BenchmarkError(f"{context} has no gene symbols")
    return tuple(genes)


def read_truth_tsv(path: Path) -> list[TruthSample]:
    if not path.is_file():
        raise BenchmarkError(f"truth TSV not found: {path}")
    try:
        with path.open("r", encoding="utf-8-sig", newline="") as handle:
            reader = csv.reader(handle, delimiter="\t", strict=True)
            rows = [(line_number, row) for line_number, row in enumerate(reader, 1) if any(row)]
    except (OSError, UnicodeError, csv.Error) as exc:
        raise BenchmarkError(f"could not parse truth TSV '{path}': {exc}") from exc

    if not rows:
        raise BenchmarkError(f"truth TSV is empty: {path}")
    header = [value.strip() for value in rows[0][1]]
    if len(header) != len(set(header)):
        raise BenchmarkError("truth TSV contains duplicate column names")
    missing = [name for name in ("sample_id", "genes") if name not in header]
    if missing:
        raise BenchmarkError(
            "truth TSV is missing required column(s): " + ", ".join(missing)
        )
    sample_index = header.index("sample_id")
    genes_index = header.index("genes")

    sample_genes: dict[str, list[str]] = {}
    gene_rows: dict[str, dict[str, int]] = {}
    for line_number, row in rows[1:]:
        if len(row) != len(header):
            raise BenchmarkError(
                f"truth TSV row {line_number} has {len(row)} columns; "
                f"expected {len(header)}"
            )
        sample_id = row[sample_index].strip()
        if not sample_id:
            raise BenchmarkError(f"truth TSV row {line_number} has an empty sample_id")
        genes = parse_truth_genes(
            row[genes_index], f"truth TSV sample '{sample_id}' at row {line_number}"
        )
        combined_genes = sample_genes.setdefault(sample_id, [])
        sample_gene_rows = gene_rows.setdefault(sample_id, {})
        for gene in genes:
            if gene in sample_gene_rows:
                raise BenchmarkError(
                    f"truth TSV rows {sample_gene_rows[gene]} and {line_number} "
                    f"have duplicate gene '{gene}' for sample_id '{sample_id}'"
                )
            combined_genes.append(gene)
            sample_gene_rows[gene] = line_number

    if not sample_genes:
        raise BenchmarkError("truth TSV has no sample rows")
    return [
        TruthSample(sample_id, tuple(genes))
        for sample_id, genes in sample_genes.items()
    ]


def validate_template(template: str, option_name: str) -> None:
    fields: list[str] = []
    try:
        for _, field_name, format_spec, conversion in string.Formatter().parse(template):
            if field_name is None:
                continue
            fields.append(field_name)
            if format_spec or conversion:
                raise BenchmarkError(
                    f"{option_name} does not support format specifications or conversions"
                )
    except ValueError as exc:
        raise BenchmarkError(f"invalid {option_name}: {exc}") from exc
    if fields != ["sample"]:
        raise BenchmarkError(
            f"{option_name} must contain exactly one '{{sample}}' placeholder"
        )


def resolve_result_path(root: Path, template: str, sample_id: str) -> Path:
    root_resolved = root.resolve()
    candidate = (root_resolved / template.format(sample=sample_id)).resolve()
    try:
        candidate.relative_to(root_resolved)
    except ValueError as exc:
        raise BenchmarkError(
            f"sample '{sample_id}' resolves outside result directory using template "
            f"'{template}'"
        ) from exc
    return candidate


def parse_positive_rank(raw_rank: str, context: str) -> int:
    try:
        rank = int(raw_rank.strip())
    except ValueError as exc:
        raise RankingError(f"{context} has non-integer rank '{raw_rank.strip()}'") from exc
    if rank < 1:
        raise RankingError(f"{context} has non-positive rank {rank}")
    return rank


def read_pipevar_ranking(path: Path) -> tuple[RankedGene, ...]:
    try:
        with path.open("r", encoding="utf-8-sig", newline="") as handle:
            lines = [(line_number, line.rstrip("\r\n")) for line_number, line in enumerate(handle, 1)]
    except (OSError, UnicodeError) as exc:
        raise RankingError(f"could not read file: {exc}") from exc

    nonblank = [(line_number, line) for line_number, line in lines if line.strip()]
    if not nonblank:
        raise RankingError("file is empty and has no header")
    header_line_number, header_line = nonblank[0]
    header = [value.strip() for value in header_line.split("\t")]
    missing = [name for name in ("#RANK", "GENE") if name not in header]
    if missing:
        raise RankingError("missing required column(s): " + ", ".join(missing))
    rank_index = header.index("#RANK")
    gene_index = header.index("GENE")
    vcf_line_index = header.index("VCF_LINE") if "VCF_LINE" in header else None
    max_required_index = max(rank_index, gene_index)

    rankings: list[RankedGene] = []
    seen_ranks: dict[int, int] = {}
    seen_genes: dict[str, int] = {}
    for line_number, line in nonblank[1:]:
        maxsplit = vcf_line_index if vcf_line_index is not None else len(header) - 1
        fields = line.split("\t", maxsplit=maxsplit)
        if len(fields) <= max_required_index:
            raise RankingError(
                f"row {line_number} has too few columns for #RANK and GENE"
            )
        rank = parse_positive_rank(fields[rank_index], f"row {line_number}")
        gene = normalize_gene(fields[gene_index])
        if not gene:
            raise RankingError(f"row {line_number} has an empty GENE")
        if rank in seen_ranks:
            raise RankingError(
                f"rows {seen_ranks[rank]} and {line_number} have duplicate rank {rank}"
            )
        if gene in seen_genes:
            raise RankingError(
                f"rows {seen_genes[gene]} and {line_number} have duplicate gene '{gene}'"
            )
        seen_ranks[rank] = line_number
        seen_genes[gene] = line_number
        rankings.append(RankedGene(gene, rank, rank))

    rankings.sort(key=lambda item: item.rank)
    return tuple(rankings)


def read_exomiser_ranking(path: Path) -> tuple[RankedGene, ...]:
    try:
        with path.open("r", encoding="utf-8-sig", newline="") as handle:
            reader = csv.reader(handle, delimiter="\t", strict=True)
            rows = [(line_number, row) for line_number, row in enumerate(reader, 1) if any(row)]
    except (OSError, UnicodeError, csv.Error) as exc:
        raise RankingError(f"could not parse file: {exc}") from exc

    if not rows:
        raise RankingError("file is empty and has no header")
    header = [value.strip() for value in rows[0][1]]
    missing = [name for name in ("#RANK", "GENE_SYMBOL") if name not in header]
    if missing:
        raise RankingError("missing required column(s): " + ", ".join(missing))
    rank_index = header.index("#RANK")
    gene_index = header.index("GENE_SYMBOL")
    max_required_index = max(rank_index, gene_index)

    raw_rows: list[tuple[int, str, int]] = []
    for line_number, row in rows[1:]:
        if len(row) <= max_required_index:
            raise RankingError(
                f"row {line_number} has too few columns for #RANK and GENE_SYMBOL"
            )
        rank = parse_positive_rank(row[rank_index], f"row {line_number}")
        gene = normalize_gene(row[gene_index])
        if not gene:
            raise RankingError(f"row {line_number} has an empty GENE_SYMBOL")
        raw_rows.append((rank, gene, line_number))

    raw_rows.sort(key=lambda item: item[0])
    rankings: list[RankedGene] = []
    seen_genes: set[str] = set()
    for raw_rank, gene, _ in raw_rows:
        if gene in seen_genes:
            continue
        seen_genes.add(gene)
        # Exomiser #RANK is the authoritative GeneScore rank. A gene can
        # appear more than once for different MOIs, so retain its first/best
        # native rank without compressing the remaining genes into a new list.
        rankings.append(RankedGene(gene, raw_rank, raw_rank))
    return tuple(rankings)


def parse_rankvar_rank(raw_rank: str, context: str) -> int:
    text = raw_rank.strip()
    try:
        numeric_rank = Decimal(text)
    except InvalidOperation as exc:
        raise RankingError(f"{context} has non-numeric rank '{text}'") from exc
    if not numeric_rank.is_finite() or numeric_rank != numeric_rank.to_integral_value():
        raise RankingError(f"{context} has non-integral rank '{text}'")
    rank = int(numeric_rank)
    if rank < 1:
        raise RankingError(f"{context} has non-positive rank {rank}")
    return rank


def parse_rankvar_genes(raw_genes: str, context: str) -> tuple[str, ...]:
    genes: list[str] = []
    seen: set[str] = set()
    for value in re.split(r"[,|/&;]+", raw_genes):
        gene = normalize_gene(value)
        if not gene or gene in {".", "NA", "N/A", "NONE", "NULL", "UNKNOWN"}:
            continue
        if gene not in seen:
            genes.append(gene)
            seen.add(gene)
    if not genes:
        raise RankingError(f"{context} has no usable Gene.refGene symbols")
    return tuple(genes)


def read_rankvar_ranking(path: Path) -> tuple[RankedGene, ...]:
    try:
        with path.open("r", encoding="utf-8-sig", newline="") as handle:
            reader = csv.reader(handle, delimiter="\t", strict=True)
            rows = [
                (line_number, row)
                for line_number, row in enumerate(reader, 1)
                if any(row)
            ]
    except (OSError, UnicodeError, csv.Error) as exc:
        raise RankingError(f"could not parse file: {exc}") from exc

    if not rows:
        raise RankingError("file is empty and has no header")
    header = [value.strip() for value in rows[0][1]]
    missing = [name for name in ("Gene.refGene", "rank") if name not in header]
    if missing:
        raise RankingError("missing required column(s): " + ", ".join(missing))
    gene_index = header.index("Gene.refGene")
    rank_index = header.index("rank")
    max_required_index = max(gene_index, rank_index)

    raw_rows: list[tuple[int, str, int]] = []
    for line_number, row in rows[1:]:
        if len(row) <= max_required_index:
            raise RankingError(
                f"row {line_number} has too few columns for Gene.refGene and rank"
            )
        raw_rank = parse_rankvar_rank(row[rank_index], f"row {line_number}")
        genes = parse_rankvar_genes(row[gene_index], f"row {line_number}")
        for gene in genes:
            raw_rows.append((raw_rank, gene, line_number))

    raw_rows.sort(key=lambda item: (item[0], item[2]))
    rankings: list[RankedGene] = []
    seen_genes: set[str] = set()
    for raw_rank, gene, _ in raw_rows:
        if gene in seen_genes:
            continue
        seen_genes.add(gene)
        rankings.append(RankedGene(gene, len(rankings) + 1, raw_rank))
    return tuple(rankings)


def safe_read_ranking(
    path: Path, parser: Callable[[Path], tuple[RankedGene, ...]]
) -> RankingResult:
    if not path.exists():
        return RankingResult("missing", (), "file not found")
    if not path.is_file():
        return RankingResult("malformed", (), "path is not a regular file")
    try:
        rankings = parser(path)
    except RankingError as exc:
        return RankingResult("malformed", (), str(exc))
    return RankingResult("ok" if rankings else "empty", rankings)


def find_best_match(
    truth_genes: Iterable[str], rankings: Iterable[RankedGene]
) -> MatchResult:
    truth_set = set(truth_genes)
    matches = [item for item in rankings if item.gene in truth_set]
    if not matches:
        return MatchResult()
    best = min(matches, key=lambda item: (item.rank, item.raw_rank, item.gene))
    return MatchResult(best.gene, best.rank, best.raw_rank)


def evaluate_samples(
    samples: Iterable[TruthSample],
    pipevar_dir: Path,
    exomiser_dir: Path,
    pipevar_template: str,
    exomiser_template: str,
    rankvar_dir: Path | None = None,
    rankvar_template: str = "{sample}.norm.alt.rank_var.tsv",
) -> list[SampleResult]:
    effective_rankvar_dir = rankvar_dir or pipevar_dir
    results: list[SampleResult] = []
    for truth in samples:
        pipevar_path = resolve_result_path(pipevar_dir, pipevar_template, truth.sample_id)
        rankvar_path = resolve_result_path(
            effective_rankvar_dir, rankvar_template, truth.sample_id
        )
        exomiser_path = resolve_result_path(
            exomiser_dir, exomiser_template, truth.sample_id
        )
        pipevar_result = safe_read_ranking(pipevar_path, read_pipevar_ranking)
        rankvar_result = safe_read_ranking(rankvar_path, read_rankvar_ranking)
        exomiser_result = safe_read_ranking(exomiser_path, read_exomiser_ranking)
        results.append(
            SampleResult(
                truth=truth,
                pipevar_path=pipevar_path,
                rankvar_path=rankvar_path,
                exomiser_path=exomiser_path,
                pipevar_result=pipevar_result,
                rankvar_result=rankvar_result,
                exomiser_result=exomiser_result,
                pipevar_match=find_best_match(truth.genes, pipevar_result.rankings),
                rankvar_match=find_best_match(truth.genes, rankvar_result.rankings),
                exomiser_match=find_best_match(truth.genes, exomiser_result.rankings),
            )
        )
    return results


def hit_at_cutoff(match: MatchResult, cutoff: int) -> bool:
    return match.rank is not None and match.rank <= cutoff


def display_optional(value: object | None) -> object:
    return "" if value is None else value


def write_detail_report(path: Path, results: Iterable[SampleResult]) -> None:
    fieldnames = [
        "sample_id",
        "truth_genes",
        "pipevar_path",
        "rankvar_path",
        "exomiser_path",
        "pipevar_status",
        "pipevar_error",
        "pipevar_matched_gene",
        "pipevar_best_rank",
        "pipevar_top1",
        "pipevar_top5",
        "pipevar_top10",
        "pipevar_top20",
        "rankvar_status",
        "rankvar_error",
        "rankvar_matched_gene",
        "rankvar_best_unique_gene_rank",
        "rankvar_best_raw_variant_rank",
        "rankvar_top1",
        "rankvar_top5",
        "rankvar_top10",
        "rankvar_top20",
        "exomiser_status",
        "exomiser_error",
        "exomiser_matched_gene",
        "exomiser_best_rank",
        "exomiser_best_raw_rank",
        "exomiser_top1",
        "exomiser_top5",
        "exomiser_top10",
        "exomiser_top20",
        "included_in_paired_summary",
    ]
    try:
        with path.open("x", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(
                handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n"
            )
            writer.writeheader()
            for result in results:
                pipevar_valid = result.pipevar_result.status in VALID_RESULT_STATUSES
                rankvar_valid = result.rankvar_result.status in VALID_RESULT_STATUSES
                exomiser_valid = result.exomiser_result.status in VALID_RESULT_STATUSES
                row: dict[str, object] = {
                    "sample_id": result.truth.sample_id,
                    "truth_genes": ";".join(result.truth.genes),
                    "pipevar_path": result.pipevar_path,
                    "rankvar_path": result.rankvar_path,
                    "exomiser_path": result.exomiser_path,
                    "pipevar_status": result.pipevar_result.status,
                    "pipevar_error": result.pipevar_result.error,
                    "pipevar_matched_gene": result.pipevar_match.matched_gene,
                    "pipevar_best_rank": display_optional(result.pipevar_match.rank),
                    "rankvar_status": result.rankvar_result.status,
                    "rankvar_error": result.rankvar_result.error,
                    "rankvar_matched_gene": result.rankvar_match.matched_gene,
                    "rankvar_best_unique_gene_rank": display_optional(
                        result.rankvar_match.rank
                    ),
                    "rankvar_best_raw_variant_rank": display_optional(
                        result.rankvar_match.raw_rank
                    ),
                    "exomiser_status": result.exomiser_result.status,
                    "exomiser_error": result.exomiser_result.error,
                    "exomiser_matched_gene": result.exomiser_match.matched_gene,
                    "exomiser_best_rank": display_optional(
                        result.exomiser_match.rank
                    ),
                    "exomiser_best_raw_rank": display_optional(
                        result.exomiser_match.raw_rank
                    ),
                    "included_in_paired_summary": int(result.paired),
                }
                for cutoff in CUTOFFS:
                    row[f"pipevar_top{cutoff}"] = (
                        int(hit_at_cutoff(result.pipevar_match, cutoff))
                        if pipevar_valid
                        else ""
                    )
                    row[f"rankvar_top{cutoff}"] = (
                        int(hit_at_cutoff(result.rankvar_match, cutoff))
                        if rankvar_valid
                        else ""
                    )
                    row[f"exomiser_top{cutoff}"] = (
                        int(hit_at_cutoff(result.exomiser_match, cutoff))
                        if exomiser_valid
                        else ""
                    )
                writer.writerow(row)
    except OSError as exc:
        raise BenchmarkError(f"could not write detail report '{path}': {exc}") from exc


def format_top_rankings(
    rankings: Iterable[RankedGene], cutoff: int = 20
) -> str:
    return ";".join(
        f"{ranking.rank}:{ranking.gene}"
        for ranking in rankings
        if ranking.rank <= cutoff
    )


def not_found_category(result: SampleResult, cutoff: int = CUTOFFS[-1]) -> str:
    """Return which primary tool missed the truth gene at the requested cutoff."""
    pipevar_missing = not hit_at_cutoff(result.pipevar_match, cutoff)
    exomiser_missing = not hit_at_cutoff(result.exomiser_match, cutoff)
    if pipevar_missing and exomiser_missing:
        return "both"
    if pipevar_missing:
        return "pipevar"
    if exomiser_missing:
        return "exomiser"
    return ""


def write_not_found_report(path: Path, results: Iterable[SampleResult]) -> int:
    fieldnames = [
        "sample_id",
        "truth_genes",
        "not_in_top20",
        "pipevar_truth_gene_top20",
        "exomiser_truth_gene_top20",
        "rankvar_truth_gene_top20",
        "pipevar_truth_gene_found_anywhere",
        "exomiser_truth_gene_found_anywhere",
        "rankvar_truth_gene_found_anywhere",
        "pipevar_status",
        "pipevar_error",
        "pipevar_path",
        "pipevar_matched_gene",
        "pipevar_best_rank",
        "pipevar_top20_genes",
        "rankvar_status",
        "rankvar_error",
        "rankvar_path",
        "rankvar_matched_gene",
        "rankvar_best_unique_gene_rank",
        "rankvar_top20_genes",
        "exomiser_status",
        "exomiser_error",
        "exomiser_path",
        "exomiser_matched_gene",
        "exomiser_best_rank",
        "exomiser_top20_genes",
    ]
    diagnostic_results = [
        result for result in results if not_found_category(result)
    ]
    try:
        with path.open("x", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(
                handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n"
            )
            writer.writeheader()
            for result in diagnostic_results:
                writer.writerow(
                    {
                        "sample_id": result.truth.sample_id,
                        "truth_genes": ";".join(result.truth.genes),
                        "not_in_top20": not_found_category(result),
                        "pipevar_truth_gene_top20": int(
                            hit_at_cutoff(result.pipevar_match, 20)
                        ),
                        "exomiser_truth_gene_top20": int(
                            hit_at_cutoff(result.exomiser_match, 20)
                        ),
                        "rankvar_truth_gene_top20": int(
                            hit_at_cutoff(result.rankvar_match, 20)
                        ),
                        "pipevar_truth_gene_found_anywhere": int(
                            result.pipevar_match.rank is not None
                        ),
                        "exomiser_truth_gene_found_anywhere": int(
                            result.exomiser_match.rank is not None
                        ),
                        "rankvar_truth_gene_found_anywhere": int(
                            result.rankvar_match.rank is not None
                        ),
                        "pipevar_status": result.pipevar_result.status,
                        "pipevar_error": result.pipevar_result.error,
                        "pipevar_path": result.pipevar_path,
                        "pipevar_matched_gene": result.pipevar_match.matched_gene,
                        "pipevar_best_rank": display_optional(
                            result.pipevar_match.rank
                        ),
                        "pipevar_top20_genes": format_top_rankings(
                            result.pipevar_result.rankings
                        ),
                        "rankvar_status": result.rankvar_result.status,
                        "rankvar_error": result.rankvar_result.error,
                        "rankvar_path": result.rankvar_path,
                        "rankvar_matched_gene": result.rankvar_match.matched_gene,
                        "rankvar_best_unique_gene_rank": display_optional(
                            result.rankvar_match.rank
                        ),
                        "rankvar_top20_genes": format_top_rankings(
                            result.rankvar_result.rankings
                        ),
                        "exomiser_status": result.exomiser_result.status,
                        "exomiser_error": result.exomiser_result.error,
                        "exomiser_path": result.exomiser_path,
                        "exomiser_matched_gene": result.exomiser_match.matched_gene,
                        "exomiser_best_rank": display_optional(
                            result.exomiser_match.rank
                        ),
                        "exomiser_top20_genes": format_top_rankings(
                            result.exomiser_result.rankings
                        ),
                    }
                )
    except OSError as exc:
        raise BenchmarkError(
            f"could not write not-found report '{path}': {exc}"
        ) from exc
    return len(diagnostic_results)


def build_summary_rows(results: Iterable[SampleResult]) -> list[dict[str, object]]:
    paired = [result for result in results if result.paired]
    denominator = len(paired)
    rows: list[dict[str, object]] = []
    for cutoff in CUTOFFS:
        pipevar_flags = [hit_at_cutoff(result.pipevar_match, cutoff) for result in paired]
        rankvar_flags = [hit_at_cutoff(result.rankvar_match, cutoff) for result in paired]
        exomiser_flags = [
            hit_at_cutoff(result.exomiser_match, cutoff) for result in paired
        ]
        pipevar_hits = sum(pipevar_flags)
        rankvar_hits = sum(rankvar_flags)
        exomiser_hits = sum(exomiser_flags)
        if denominator:
            pipevar_rate: object = f"{pipevar_hits / denominator:.6f}"
            rankvar_rate: object = f"{rankvar_hits / denominator:.6f}"
            exomiser_rate: object = f"{exomiser_hits / denominator:.6f}"
            rankvar_rate_difference: object = (
                f"{(pipevar_hits - rankvar_hits) * 100.0 / denominator:.6f}"
            )
            exomiser_rate_difference: object = (
                f"{(pipevar_hits - exomiser_hits) * 100.0 / denominator:.6f}"
            )
        else:
            pipevar_rate = "NA"
            rankvar_rate = "NA"
            exomiser_rate = "NA"
            rankvar_rate_difference = "NA"
            exomiser_rate_difference = "NA"
        rows.append(
            {
                "cutoff": cutoff,
                "paired_samples": denominator,
                "pipevar_hits": pipevar_hits,
                "pipevar_rate": pipevar_rate,
                "rankvar_hits": rankvar_hits,
                "rankvar_rate": rankvar_rate,
                "exomiser_hits": exomiser_hits,
                "exomiser_rate": exomiser_rate,
                "hit_difference_pipevar_minus_rankvar": pipevar_hits
                - rankvar_hits,
                "rate_difference_pipevar_minus_rankvar_percentage_points": rankvar_rate_difference,
                "hit_difference_pipevar_minus_exomiser": pipevar_hits
                - exomiser_hits,
                "rate_difference_percentage_points": exomiser_rate_difference,
                "both_hit": sum(p and e for p, e in zip(pipevar_flags, exomiser_flags)),
                "pipevar_only": sum(
                    p and not e for p, e in zip(pipevar_flags, exomiser_flags)
                ),
                "exomiser_only": sum(
                    e and not p for p, e in zip(pipevar_flags, exomiser_flags)
                ),
                "neither_hit": sum(
                    not p and not e for p, e in zip(pipevar_flags, exomiser_flags)
                ),
                "pipevar_rankvar_both_hit": sum(
                    p and r for p, r in zip(pipevar_flags, rankvar_flags)
                ),
                "pipevar_only_vs_rankvar": sum(
                    p and not r for p, r in zip(pipevar_flags, rankvar_flags)
                ),
                "rankvar_only_vs_pipevar": sum(
                    r and not p for p, r in zip(pipevar_flags, rankvar_flags)
                ),
                "pipevar_rankvar_neither_hit": sum(
                    not p and not r for p, r in zip(pipevar_flags, rankvar_flags)
                ),
            }
        )
    return rows


def write_summary_report(path: Path, rows: Iterable[dict[str, object]]) -> None:
    fieldnames = [
        "cutoff",
        "paired_samples",
        "pipevar_hits",
        "pipevar_rate",
        "rankvar_hits",
        "rankvar_rate",
        "exomiser_hits",
        "exomiser_rate",
        "hit_difference_pipevar_minus_rankvar",
        "rate_difference_pipevar_minus_rankvar_percentage_points",
        "hit_difference_pipevar_minus_exomiser",
        "rate_difference_percentage_points",
        "both_hit",
        "pipevar_only",
        "exomiser_only",
        "neither_hit",
        "pipevar_rankvar_both_hit",
        "pipevar_only_vs_rankvar",
        "rankvar_only_vs_pipevar",
        "pipevar_rankvar_neither_hit",
    ]
    try:
        with path.open("x", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(
                handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n"
            )
            writer.writeheader()
            writer.writerows(rows)
    except OSError as exc:
        raise BenchmarkError(f"could not write summary report '{path}': {exc}") from exc


def validate_run_inputs(args: argparse.Namespace) -> tuple[Path, Path, Path]:
    rankvar_dir = args.rankvar_dir or args.pipevar_dir
    for root, label in (
        (args.pipevar_dir, "PipeVar result directory"),
        (rankvar_dir, "RankVar result directory"),
        (args.exomiser_dir, "Exomiser result directory"),
    ):
        if not root.is_dir():
            raise BenchmarkError(f"{label} not found: {root}")
    validate_template(args.pipevar_template, "--pipevar-template")
    validate_template(args.rankvar_template, "--rankvar-template")
    validate_template(args.exomiser_template, "--exomiser-template")

    detail_path = Path(f"{args.output_prefix}.details.tsv")
    summary_path = Path(f"{args.output_prefix}.summary.tsv")
    not_found_path = Path(f"{args.output_prefix}.not_found.tsv")
    if not detail_path.parent.is_dir():
        raise BenchmarkError(f"output directory not found: {detail_path.parent}")
    existing = [
        str(path)
        for path in (detail_path, summary_path, not_found_path)
        if path.exists()
    ]
    if existing:
        raise BenchmarkError("output file(s) already exist: " + ", ".join(existing))
    return detail_path, summary_path, not_found_path


def print_results(
    results: Iterable[SampleResult], summary_rows: Iterable[dict[str, object]]
) -> None:
    results = list(results)
    for result in results:
        for tool, ranking_result, path in (
            ("PipeVar", result.pipevar_result, result.pipevar_path),
            ("RankVar", result.rankvar_result, result.rankvar_path),
            ("Exomiser", result.exomiser_result, result.exomiser_path),
        ):
            if ranking_result.status not in VALID_RESULT_STATUSES:
                print(
                    f"WARNING: {result.truth.sample_id} {tool} "
                    f"{ranking_result.status}: {path}: {ranking_result.error}",
                    file=sys.stderr,
                )
    paired_count = sum(result.paired for result in results)
    print(f"Paired evaluable samples: {paired_count}/{len(results)}")
    for row in summary_rows:
        denominator = int(row["paired_samples"])
        if denominator:
            pipevar_percent = 100.0 * int(row["pipevar_hits"]) / denominator
            rankvar_percent = 100.0 * int(row["rankvar_hits"]) / denominator
            exomiser_percent = 100.0 * int(row["exomiser_hits"]) / denominator
            print(
                f"Top {row['cutoff']}: PipeVar {row['pipevar_hits']}/{denominator} "
                f"({pipevar_percent:.1f}%); RankVar "
                f"{row['rankvar_hits']}/{denominator} ({rankvar_percent:.1f}%); "
                f"Exomiser "
                f"{row['exomiser_hits']}/{denominator} ({exomiser_percent:.1f}%)"
            )
        else:
            print(f"Top {row['cutoff']}: no paired evaluable samples")


def run(args: argparse.Namespace) -> int:
    detail_path, summary_path, not_found_path = validate_run_inputs(args)
    samples = read_truth_tsv(args.truth_tsv)
    results = evaluate_samples(
        samples,
        args.pipevar_dir,
        args.exomiser_dir,
        args.pipevar_template,
        args.exomiser_template,
        args.rankvar_dir,
        args.rankvar_template,
    )
    summary_rows = build_summary_rows(results)
    write_detail_report(detail_path, results)
    write_summary_report(summary_path, summary_rows)
    not_found_count = write_not_found_report(not_found_path, results)
    print_results(results, summary_rows)
    print(f"Detail report: {detail_path}")
    print(f"Summary report: {summary_path}")
    print(f"Not-found diagnostics: {not_found_path} ({not_found_count} sample(s))")
    return 0 if any(result.paired for result in results) else 1


def main(argv: Sequence[str] | None = None) -> int:
    try:
        return run(parse_args(argv))
    except BenchmarkError as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
