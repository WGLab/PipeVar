#!/usr/bin/env python3

import argparse
import gzip
import re
import sys
from collections import Counter, defaultdict
from dataclasses import dataclass
from difflib import SequenceMatcher
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple


DEFAULT_COMMON_VCF = "/opt/pipevar/resources/common_svs.vcf"
INTERVAL_SVTYPES = {"DEL", "DUP", "INV", "CNV"}
SKIP_SVTYPES = {"BND", "TRA"}


def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)


def open_maybe_gzip(path: str, mode: str = "rt"):
    if path.endswith(".gz"):
        return gzip.open(path, mode)
    return open(path, mode)


def clean_field(value: Optional[str]) -> str:
    if value is None:
        return "."
    value = str(value).strip()
    return value if value else "."


def normalize_func(value: str) -> str:
    value = clean_field(value)
    return value.lower() if value != "." else "."


def parse_info(info_str: str) -> Dict[str, str]:
    info = {}
    for part in info_str.split(";"):
        if not part:
            continue
        if "=" in part:
            key, value = part.split("=", 1)
            info[key] = value
        else:
            info[part] = "True"
    return info


def pick_first(info: Dict[str, str], keys: List[str], default: str = ".") -> str:
    for key in keys:
        value = info.get(key)
        if value not in (None, ""):
            return value
    return default


def parse_int_or_none(value: Optional[str]) -> Optional[int]:
    if value is None or value in ("", ".", "NA"):
        return None
    try:
        return int(str(value).split(",")[0])
    except ValueError:
        return None


def parse_float_max(value: Optional[str]) -> Optional[float]:
    if value is None or value in ("", ".", "NA"):
        return None
    vals = []
    for item in str(value).split(","):
        try:
            vals.append(float(item))
        except ValueError:
            continue
    return max(vals) if vals else None


def normalize_chrom(chrom: str) -> str:
    chrom = chrom.strip()
    lower = chrom.lower()
    if lower.startswith("chr"):
        lower = lower[3:]
    if lower in {"m", "mt"}:
        return "mt"
    return lower


def infer_svtype(alt: str, info: Dict[str, str]) -> str:
    svtype = clean_field(info.get("SVTYPE"))
    if svtype != ".":
        return svtype.upper()
    first_alt = alt.split(",", 1)[0]
    if first_alt.startswith("<") and first_alt.endswith(">"):
        return first_alt[1:-1].upper()
    if "[" in first_alt or "]" in first_alt:
        return "BND"
    return "UNK"


def infer_end(pos: int, svtype: str, info: Dict[str, str]) -> int:
    end = parse_int_or_none(info.get("END"))
    if end is not None:
        return end
    svlen = parse_int_or_none(info.get("SVLEN"))
    if svlen is not None:
        if svtype == "INS":
            return pos
        return pos + max(0, abs(svlen) - 1)
    return pos


def infer_svlen(pos: int, end: int, svtype: str, info: Dict[str, str]) -> int:
    svlen = parse_int_or_none(info.get("SVLEN"))
    if svlen is not None:
        return abs(svlen)
    if svtype == "INS":
        return 1
    return max(1, end - pos + 1)


def extract_annovar_fields(info: Dict[str, str], gene_model: str) -> Tuple[str, str, str, str]:
    func = pick_first(
        info,
        ["Func.{0}".format(gene_model), "Func.refGene", "Func.ensGene", "Func.knownGene", "FUNC"],
    )
    gene = pick_first(
        info,
        ["Gene.{0}".format(gene_model), "Gene.refGene", "Gene.ensGene", "Gene.knownGene", "GENE"],
    )
    gene_detail = pick_first(
        info,
        [
            "GeneDetail.{0}".format(gene_model),
            "GeneDetail.refGene",
            "GeneDetail.ensGene",
            "GeneDetail.knownGene",
            "GENEDETAIL",
        ],
    )
    exonic_func = pick_first(
        info,
        [
            "ExonicFunc.{0}".format(gene_model),
            "ExonicFunc.refGene",
            "ExonicFunc.ensGene",
            "ExonicFunc.knownGene",
            "EXONICFUNC",
        ],
    )
    return normalize_func(func), clean_field(gene), clean_field(gene_detail), normalize_func(exonic_func)


def reciprocal_overlap(start1: int, end1: int, start2: int, end2: int) -> float:
    left = max(start1, start2)
    right = min(end1, end2)
    if right < left:
        return 0.0
    overlap = right - left + 1
    len1 = max(1, end1 - start1 + 1)
    len2 = max(1, end2 - start2 + 1)
    return min(overlap / len1, overlap / len2)


def inserted_sequence(ref: str, alt: str) -> Optional[str]:
    first_alt = alt.split(",", 1)[0].upper()
    ref = ref.upper()
    if not re.fullmatch(r"[ACGTN]+", first_alt):
        return None
    if first_alt.startswith(ref) and len(first_alt) > len(ref):
        return first_alt[len(ref):]
    if len(first_alt) > len(ref):
        return first_alt
    return None


def sequence_identity(seq1: Optional[str], seq2: Optional[str]) -> Optional[float]:
    if not seq1 or not seq2:
        return None
    return SequenceMatcher(None, seq1, seq2, autojunk=False).ratio()


def safe_info_value(value: object) -> str:
    value = clean_field(str(value))
    return (
        value.replace("%", "%25")
        .replace(";", "%3B")
        .replace(",", "%2C")
        .replace(" ", "_")
        .replace("\t", "_")
    )


def append_info(info_str: str, values: Dict[str, object]) -> str:
    existing = "." if info_str in ("", ".") else info_str
    parts = []
    for key, value in values.items():
        if value is True:
            parts.append(key)
        else:
            parts.append("{0}={1}".format(key, safe_info_value(value)))
    extra = ";".join(parts)
    if existing == ".":
        return extra
    return existing + ";" + extra


@dataclass
class SVRecord:
    chrom: str
    norm_chrom: str
    pos: int
    end: int
    svtype: str
    svlen: int
    ref: str
    alt: str
    record_id: str
    info: Dict[str, str]
    func: str
    gene: str
    gene_detail: str
    exonic_func: str
    raw_cols: List[str]
    raw_line: str


@dataclass
class CommonSV:
    chrom: str
    norm_chrom: str
    pos: int
    end: int
    svtype: str
    svlen: int
    ref: str
    alt: str
    record_id: str
    af: float
    ac: str
    an: str
    support: str
    start_min: int
    start_max: int
    end_min: int
    end_max: int
    func: str
    gene: str
    gene_detail: str
    exonic_func: str


@dataclass
class MatchResult:
    common: CommonSV
    method: str
    ro: float
    start_diff: int
    end_diff: int
    len_diff: int
    ins_identity: Optional[float]


def record_from_cols(cols: List[str], gene_model: str, raw_line: str) -> Optional[SVRecord]:
    if len(cols) < 8:
        return None
    try:
        pos = int(cols[1])
    except ValueError:
        return None
    info = parse_info(cols[7])
    svtype = infer_svtype(cols[4], info)
    end = infer_end(pos, svtype, info)
    svlen = infer_svlen(pos, end, svtype, info)
    func, gene, gene_detail, exonic_func = extract_annovar_fields(info, gene_model)
    return SVRecord(
        chrom=cols[0],
        norm_chrom=normalize_chrom(cols[0]),
        pos=pos,
        end=end,
        svtype=svtype,
        svlen=svlen,
        ref=cols[3],
        alt=cols[4],
        record_id=cols[2],
        info=info,
        func=func,
        gene=gene,
        gene_detail=gene_detail,
        exonic_func=exonic_func,
        raw_cols=cols,
        raw_line=raw_line,
    )


def common_from_record(rec: SVRecord, min_af: float) -> Optional[CommonSV]:
    af = parse_float_max(rec.info.get("AF"))
    if af is None or af < min_af:
        return None
    start_min = parse_int_or_none(rec.info.get("STARTMIN")) or rec.pos
    start_max = parse_int_or_none(rec.info.get("STARTMAX")) or rec.pos
    end_min = parse_int_or_none(rec.info.get("ENDMIN")) or rec.end
    end_max = parse_int_or_none(rec.info.get("ENDMAX")) or rec.end
    return CommonSV(
        chrom=rec.chrom,
        norm_chrom=rec.norm_chrom,
        pos=rec.pos,
        end=rec.end,
        svtype=rec.svtype,
        svlen=rec.svlen,
        ref=rec.ref,
        alt=rec.alt,
        record_id=rec.record_id,
        af=af,
        ac=clean_field(rec.info.get("AC")),
        an=clean_field(rec.info.get("AN")),
        support=clean_field(rec.info.get("SUPPORT")),
        start_min=min(start_min, start_max),
        start_max=max(start_min, start_max),
        end_min=min(end_min, end_max),
        end_max=max(end_min, end_max),
        func=rec.func,
        gene=rec.gene,
        gene_detail=rec.gene_detail,
        exonic_func=rec.exonic_func,
    )


def load_common_vcf(path: str, min_af: float, gene_model: str) -> Dict[Tuple[str, str], List[CommonSV]]:
    common = defaultdict(list)
    with open_maybe_gzip(path) as fh:
        for raw in fh:
            line = raw.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            cols = line.split("\t")
            rec = record_from_cols(cols[:8], gene_model, line)
            if rec is None or rec.svtype in SKIP_SVTYPES:
                continue
            common_rec = common_from_record(rec, min_af)
            if common_rec is not None:
                common[(common_rec.norm_chrom, common_rec.svtype)].append(common_rec)
    for key in common:
        common[key].sort(key=lambda r: (r.start_min, r.end_min, r.pos))
    return common


def interval_match(
    rec: SVRecord,
    common: CommonSV,
    reciprocal_overlap_threshold: float,
    distance: int,
) -> Optional[MatchResult]:
    ro = reciprocal_overlap(rec.pos, rec.end, common.pos, common.end)
    start_diff = min(abs(rec.pos - common.pos), abs(rec.pos - common.start_min), abs(rec.pos - common.start_max))
    end_diff = min(abs(rec.end - common.end), abs(rec.end - common.end_min), abs(rec.end - common.end_max))
    len_diff = abs(rec.svlen - common.svlen)
    if ro >= reciprocal_overlap_threshold:
        return MatchResult(common, "RECIPROCAL_OVERLAP", ro, start_diff, end_diff, len_diff, None)
    if start_diff <= distance and end_diff <= distance:
        return MatchResult(common, "BREAKPOINT_DISTANCE", ro, start_diff, end_diff, len_diff, None)
    return None


def insertion_match(
    rec: SVRecord,
    common: CommonSV,
    ins_distance: int,
    ins_len_frac: float,
    min_ins_len_diff: int,
    ins_identity_threshold: float,
) -> Optional[MatchResult]:
    start_diff = min(abs(rec.pos - common.pos), abs(rec.pos - common.start_min), abs(rec.pos - common.start_max))
    end_diff = abs(rec.end - common.end)
    len_diff = abs(rec.svlen - common.svlen)
    allowed_len_diff = max(min_ins_len_diff, int(round(ins_len_frac * max(rec.svlen, common.svlen))))
    if start_diff > ins_distance or len_diff > allowed_len_diff:
        return None
    identity = sequence_identity(inserted_sequence(rec.ref, rec.alt), inserted_sequence(common.ref, common.alt))
    if identity is None:
        return MatchResult(common, "INS_POS_LEN_ONLY", 0.0, start_diff, end_diff, len_diff, None)
    if identity >= ins_identity_threshold:
        return MatchResult(common, "INS_SEQUENCE_IDENTITY", 0.0, start_diff, end_diff, len_diff, identity)
    return None


def find_match(rec: SVRecord, common_index: Dict[Tuple[str, str], List[CommonSV]], args) -> Optional[MatchResult]:
    if rec.svtype in SKIP_SVTYPES:
        return None
    candidates = common_index.get((rec.norm_chrom, rec.svtype), [])
    best = None
    best_key = None
    for common in candidates:
        if rec.svtype == "INS":
            match = insertion_match(
                rec,
                common,
                args.ins_distance,
                args.ins_len_frac,
                args.min_ins_len_diff,
                args.ins_identity,
            )
        elif rec.svtype in INTERVAL_SVTYPES or rec.end != rec.pos:
            match = interval_match(rec, common, args.reciprocal_overlap, args.distance)
        else:
            match = interval_match(rec, common, args.reciprocal_overlap, args.distance)
        if match is None:
            continue
        key = (
            match.start_diff,
            match.end_diff,
            match.len_diff,
            -match.ro,
            -(match.ins_identity if match.ins_identity is not None else 0.0),
            -match.common.af,
        )
        if best is None or key < best_key:
            best = match
            best_key = key
    return best


def write_header_with_common_info(header_lines: List[str], out) -> None:
    inserted = False
    for line in header_lines:
        if line.startswith("#CHROM") and not inserted:
            out.write('##INFO=<ID=COMMON_SV,Number=0,Type=Flag,Description="Variant matched the common SV filter database">\n')
            out.write('##INFO=<ID=COMMON_SV_ID,Number=1,Type=String,Description="Matched common SV record ID">\n')
            out.write('##INFO=<ID=COMMON_SV_AF,Number=1,Type=Float,Description="Matched common SV allele frequency">\n')
            out.write('##INFO=<ID=COMMON_SV_MATCH,Number=1,Type=String,Description="Common SV matching method">\n')
            out.write('##INFO=<ID=COMMON_SV_RO,Number=1,Type=Float,Description="Reciprocal overlap with matched common SV">\n')
            out.write('##INFO=<ID=COMMON_SV_INS_IDENTITY,Number=1,Type=Float,Description="Insertion sequence identity with matched common SV when available">\n')
            inserted = True
        out.write(line + "\n")


def write_summary(path: str, total: int, kept: int, removed: int, by_type: Counter, args) -> None:
    with open(path, "w") as out:
        out.write("metric\tvalue\n")
        out.write("total_records\t{0}\n".format(total))
        out.write("kept_records\t{0}\n".format(kept))
        out.write("removed_records\t{0}\n".format(removed))
        for svtype, count in sorted(by_type.items()):
            out.write("removed_{0}\t{1}\n".format(svtype, count))
        out.write("common_sv_af\t{0}\n".format(args.common_af))
        out.write("common_sv_reciprocal_overlap\t{0}\n".format(args.reciprocal_overlap))
        out.write("common_sv_distance\t{0}\n".format(args.distance))
        out.write("common_sv_ins_distance\t{0}\n".format(args.ins_distance))
        out.write("common_sv_ins_identity\t{0}\n".format(args.ins_identity))


def filter_vcf(args) -> None:
    common_index = load_common_vcf(args.common_vcf, args.common_af, args.gene_model)
    total_common = sum(len(records) for records in common_index.values())
    eprint("[INFO] Loaded {0} common SV records from {1}".format(total_common, args.common_vcf))

    header_lines = []
    total = kept = removed = 0
    removed_by_type = Counter()

    with open_maybe_gzip(args.input_vcf) as inp, open(args.output_vcf, "w") as kept_out, open(args.removed_vcf, "w") as removed_out:
        kept_header_written = False
        removed_header_written = False
        for raw in inp:
            line = raw.rstrip("\n")
            if line.startswith("#"):
                header_lines.append(line)
                continue
            if not kept_header_written:
                for header in header_lines:
                    kept_out.write(header + "\n")
                kept_header_written = True
            if not removed_header_written:
                write_header_with_common_info(header_lines, removed_out)
                removed_header_written = True

            if not line:
                continue
            cols = line.split("\t")
            rec = record_from_cols(cols, args.gene_model, line)
            if rec is None:
                kept_out.write(line + "\n")
                continue
            total += 1
            match = find_match(rec, common_index, args)
            if match is None:
                kept += 1
                kept_out.write(line + "\n")
                continue

            removed += 1
            removed_by_type[rec.svtype] += 1
            annotated_cols = list(cols)
            annotated_cols[7] = append_info(
                annotated_cols[7],
                {
                    "COMMON_SV": True,
                    "COMMON_SV_ID": match.common.record_id,
                    "COMMON_SV_AF": "{0:.6g}".format(match.common.af),
                    "COMMON_SV_MATCH": match.method,
                    "COMMON_SV_RO": "{0:.6g}".format(match.ro),
                    "COMMON_SV_INS_IDENTITY": "." if match.ins_identity is None else "{0:.6g}".format(match.ins_identity),
                },
            )
            removed_out.write("\t".join(annotated_cols) + "\n")

        if not kept_header_written:
            for header in header_lines:
                kept_out.write(header + "\n")
        if not removed_header_written:
            write_header_with_common_info(header_lines, removed_out)

    write_summary(args.summary_tsv, total, kept, removed, removed_by_type, args)
    eprint("[INFO] total={0} kept={1} removed={2}".format(total, kept, removed))


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Filter common structural variants from an ANNOVAR SV VCF.")
    parser.add_argument("--input-vcf", required=True, help="Input ANNOVAR SV VCF")
    parser.add_argument("--output-vcf", required=True, help="Output VCF with common SVs removed")
    parser.add_argument("--removed-vcf", required=True, help="Output VCF containing removed common SVs")
    parser.add_argument("--summary-tsv", required=True, help="Summary TSV")
    parser.add_argument("--common-vcf", default=DEFAULT_COMMON_VCF, help="Common SV database VCF")
    parser.add_argument("--common-af", type=float, default=0.01, help="Minimum AF for common database records")
    parser.add_argument("--gene-model", default="refGene", help="Preferred ANNOVAR gene model")
    parser.add_argument("--reciprocal-overlap", type=float, default=0.5, help="Minimum reciprocal overlap for interval SVs")
    parser.add_argument("--distance", type=int, default=1000, help="Breakpoint fallback distance for interval SVs")
    parser.add_argument("--ins-distance", type=int, default=500, help="Insertion position window")
    parser.add_argument("--ins-len-frac", type=float, default=0.5, help="Insertion length tolerance fraction")
    parser.add_argument("--min-ins-len-diff", type=int, default=50, help="Minimum insertion length tolerance")
    parser.add_argument("--ins-identity", type=float, default=0.8, help="Insertion sequence identity threshold")
    parser.add_argument("--version", action="version", version="filter_common_svs.py 0.1")
    return parser


def main() -> None:
    parser = build_parser()
    args = parser.parse_args()
    if not Path(args.common_vcf).exists():
        parser.error("Common SV VCF not found: {0}".format(args.common_vcf))
    filter_vcf(args)


if __name__ == "__main__":
    main()
