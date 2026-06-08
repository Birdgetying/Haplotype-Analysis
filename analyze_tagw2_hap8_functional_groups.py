#!/usr/bin/env python3
"""Check whether TaGW2-A1 Hap8 matches published functional variants.

This script is intentionally separate from the core HaplotypeScorer workflow.
It diagnoses the biological interpretation of the existing TaGW2-A1 result by
splitting Hap8 into literature-marker and background-marker components.
"""

import argparse
import io
import json
import math
import sys
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence

if sys.platform == "win32":
    try:
        if not isinstance(sys.stdout, io.TextIOWrapper) or sys.stdout.encoding != "utf-8":
            sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")
        if not isinstance(sys.stderr, io.TextIOWrapper) or sys.stderr.encoding != "utf-8":
            sys.stderr = io.TextIOWrapper(sys.stderr.buffer, encoding="utf-8", errors="replace")
    except (ValueError, AttributeError):
        pass

import numpy as np
import pandas as pd
from scipy import stats


DEFAULT_DB_DIR = Path("star_gene_database/wheat_nature_2024/TaGW2-A1")
DEFAULT_RESULT_DIR = Path("star_gene_results/wheat_nature_2024/TaGW2-A1")
DEFAULT_MARKER_MATRIX = Path("external_data/wheat_nature_2024/tagw2_ad/TaGW2-A1/marker_matrix.tsv")
DEFAULT_PHENOTYPE_XLSX = Path(
    "external_data/wheat_nature_2024/watseq/"
    "Watkins_Collection_WGIN_WISP_DFW_watseq_phenotype_data_JIC.xlsx"
)
DEFAULT_GW2_DATA_DIR = Path(r"D:\Desktop\data\GW2")

TRAIT_COLUMNS = ["TGW_mean", "TGW_CFLN10", "TGW_CFLN14"]
GENE_START = 237_734_651
CDS_START_OFFSET_1BASED = 185

SU_QIN_PAIR = {
    "chr6A_237734096": "G",  # promoter offset -739, carried by local Hap8
    "chr6A_237734242": "A",  # promoter offset -593, carried by local Hap8
}

PUBLISHED_VARIANTS = [
    {
        "variant_name": "SNP-988",
        "variant_class": "promoter_snp",
        "promoter_offset_from_atg": -988,
        "literature_alleles": "G/A",
        "jaiswal_hap5_allele": "G",
        "tested_group": "coverage_only",
        "source_note": "Jaiswal 2015 reported four promoter SNPs; Hap5 is G_A_G_A.",
    },
    {
        "variant_name": "SNP-739",
        "variant_class": "promoter_snp",
        "promoter_offset_from_atg": -739,
        "literature_alleles": "G/A",
        "jaiswal_hap5_allele": "A",
        "tested_group": "su_qin_two_snp_pair",
        "source_note": "Previously reported by Su et al.; part of the two-SNP state carried by local Hap8.",
    },
    {
        "variant_name": "SNP-593",
        "variant_class": "promoter_snp",
        "promoter_offset_from_atg": -593,
        "literature_alleles": "A/G",
        "jaiswal_hap5_allele": "G",
        "tested_group": "su_qin_two_snp_pair",
        "source_note": "Previously reported by Su et al.; part of the two-SNP state carried by local Hap8.",
    },
    {
        "variant_name": "SNP-494",
        "variant_class": "promoter_snp",
        "promoter_offset_from_atg": -494,
        "literature_alleles": "G/A",
        "jaiswal_hap5_allele": "A",
        "tested_group": "not_tested_missing_marker",
        "source_note": "Jaiswal 2015 nominated SNP-494 as causal and expression-regulating.",
    },
    {
        "variant_name": "G2373A",
        "variant_class": "splice_acceptor_snp",
        "promoter_offset_from_atg": None,
        "literature_alleles": "G/A",
        "jaiswal_hap5_allele": "",
        "tested_group": "not_tested_coordinate_not_lifted",
        "source_note": "Simmonds 2016 reported G2373A on an earlier TaGW2-A1 reference context.",
    },
    {
        "variant_name": "exon8_1bp_insertion",
        "variant_class": "coding_indel",
        "promoter_offset_from_atg": None,
        "literature_alleles": "1-bp T insertion",
        "jaiswal_hap5_allele": "",
        "tested_group": "not_tested_indel_data_missing",
        "source_note": "Reported coding indel; the current database was built from SNP VCF only.",
    },
]

CHR6A_INDEL_URL = (
    "https://opendata.earlham.ac.uk/wheat/under_license/toronto/"
    "WatSeq_2023-09-15_landrace_modern_Variation_Data/"
    "WatSeq_VCF_Raw_ChineseSpringRefSeqv1.0/chr6A/"
    "chr6A.HARD.INDEL.Missing-unphasing.ID.ann.finalSID.vcf.gz"
)
CHR6A_INDEL_INDEX_URL = CHR6A_INDEL_URL + ".csi"


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Split TaGW2-A1 Hap8 into published promoter SNP and background SNP components.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--db-dir", default=str(DEFAULT_DB_DIR), help="TaGW2-A1 precomputed database")
    parser.add_argument("--result-dir", default=str(DEFAULT_RESULT_DIR), help="TaGW2-A1 result directory")
    parser.add_argument("--marker-matrix", default=str(DEFAULT_MARKER_MATRIX), help="TaGW2-A1 marker matrix TSV")
    parser.add_argument("--phenotype-xlsx", default=str(DEFAULT_PHENOTYPE_XLSX), help="Watkins phenotype workbook")
    parser.add_argument("--gw2-data-dir", default=str(DEFAULT_GW2_DATA_DIR), help="Local GW2 data directory")
    parser.add_argument("--full-haplotype", default="Hap8", help="Full haplotype to diagnose")
    return parser


def map_promoter_offset_to_position(
    gene_start: int,
    cds_start_offset_1based: int,
    promoter_offset_from_atg: int,
) -> int:
    """Convert an ATG-relative promoter offset to a 1-based RefSeq position."""
    atg_position = int(gene_start) + int(cds_start_offset_1based) - 1
    return atg_position + int(promoter_offset_from_atg)


def _pattern_match(df: pd.DataFrame, alleles: Dict[str, str]) -> pd.Series:
    if not alleles:
        return pd.Series(False, index=df.index)
    missing = [marker for marker in alleles if marker not in df.columns]
    if missing:
        return pd.Series(False, index=df.index)
    expected = pd.Series(alleles)
    return df[list(alleles)].astype(str).eq(expected).all(axis=1)


def classify_literature_background_groups(
    df: pd.DataFrame,
    literature_alleles: Dict[str, str],
    background_alleles: Dict[str, str],
    full_haplotype_name: str = "Hap8",
) -> pd.DataFrame:
    """Add published-pair, background-pattern, and split-group labels."""
    out = df.copy()
    out["PublishedPair"] = _pattern_match(out, literature_alleles)
    out["Hap8Background"] = _pattern_match(out, background_alleles)
    out["FullHap8"] = out["Hap_Name"].astype(str).eq(full_haplotype_name)
    out["Hap8SplitGroup"] = "Neither"
    out.loc[out["PublishedPair"] & ~out["Hap8Background"], "Hap8SplitGroup"] = "PublishedOnly"
    out.loc[out["PublishedPair"] & out["Hap8Background"], "Hap8SplitGroup"] = "Published+Background"
    out.loc[~out["PublishedPair"] & out["Hap8Background"], "Hap8SplitGroup"] = "BackgroundOnly"
    return out


def summarize_binary_trait(
    df: pd.DataFrame,
    group_col: str,
    trait_col: str,
    case_label: str,
    control_label: str,
) -> Dict[str, object]:
    case = pd.to_numeric(df.loc[df[group_col].astype(bool), trait_col], errors="coerce").dropna()
    control = pd.to_numeric(df.loc[~df[group_col].astype(bool), trait_col], errors="coerce").dropna()
    pvalue = math.nan
    statistic = math.nan
    if len(case) >= 2 and len(control) >= 2:
        test = stats.ttest_ind(case, control, equal_var=False, nan_policy="omit")
        statistic = float(test.statistic)
        pvalue = float(test.pvalue)
    return {
        "group_col": group_col,
        "trait": trait_col,
        "case_label": case_label,
        "control_label": control_label,
        "case_n": int(len(case)),
        "control_n": int(len(control)),
        "case_mean": float(case.mean()) if len(case) else math.nan,
        "control_mean": float(control.mean()) if len(control) else math.nan,
        "case_sd": float(case.std(ddof=1)) if len(case) > 1 else math.nan,
        "control_sd": float(control.std(ddof=1)) if len(control) > 1 else math.nan,
        "effect_case_minus_control": (
            float(case.mean() - control.mean()) if len(case) and len(control) else math.nan
        ),
        "welch_t": statistic,
        "welch_pvalue": pvalue,
    }


def _bh_fdr(pvalues: Sequence[float]) -> List[float]:
    valid = [(idx, float(p)) for idx, p in enumerate(pvalues) if not pd.isna(p)]
    adjusted = [math.nan] * len(pvalues)
    if not valid:
        return adjusted
    valid_sorted = sorted(valid, key=lambda item: item[1])
    m = len(valid_sorted)
    running = 1.0
    for rank_from_end, (idx, pvalue) in enumerate(reversed(valid_sorted), start=1):
        rank = m - rank_from_end + 1
        running = min(running, pvalue * m / rank)
        adjusted[idx] = min(running, 1.0)
    return adjusted


def _read_required_csv(path: Path) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(path)
    return pd.read_csv(path)


def load_watkins_country_table(phenotype_xlsx: Path) -> pd.DataFrame:
    if not phenotype_xlsx.exists():
        raise FileNotFoundError(phenotype_xlsx)
    pbi = pd.read_excel(phenotype_xlsx, sheet_name="PBI data (historic)")
    required = ["StoreCode", "COUNTRY of origin", "acc_name"]
    missing = [col for col in required if col not in pbi.columns]
    if missing:
        raise ValueError(f"PBI data sheet missing columns: {', '.join(missing)}")
    out = pbi[required].copy()
    out = out.rename(columns={"StoreCode": "SampleID", "COUNTRY of origin": "Country", "acc_name": "AccessionName"})
    out["SampleID"] = out["SampleID"].astype(str).str.strip()
    out["Country"] = out["Country"].fillna("Unknown").astype(str).str.strip()
    out.loc[out["Country"].eq(""), "Country"] = "Unknown"
    return out.drop_duplicates(subset=["SampleID"])


def build_hap8_defining_snps(db_dir: Path, result_dir: Path, full_haplotype: str) -> pd.DataFrame:
    variant_info = _read_required_csv(db_dir / "variant_info.csv")
    haplotypes = _read_required_csv(db_dir / "haplotype_data.csv")
    reference = haplotypes.loc[haplotypes["Hap_Name"].eq("Hap1")]
    target = haplotypes.loc[haplotypes["Hap_Name"].eq(full_haplotype)]
    if reference.empty:
        raise ValueError("Hap1 is missing from haplotype_data.csv")
    if target.empty:
        raise ValueError(f"{full_haplotype} is missing from haplotype_data.csv")

    ref_alleles = str(reference.iloc[0]["Haplotype_Seq"]).split("|")
    target_alleles = str(target.iloc[0]["Haplotype_Seq"]).split("|")
    rows = []
    for index, (ref_allele, target_allele) in enumerate(zip(ref_alleles, target_alleles)):
        if ref_allele == target_allele:
            continue
        row = variant_info.iloc[index].to_dict()
        row["index_1based"] = index + 1
        row["Hap1_allele"] = ref_allele
        row[f"{full_haplotype}_allele"] = target_allele
        row["Hap8_allele"] = target_allele
        rows.append(row)
    out = pd.DataFrame(rows)
    out_path = result_dir / f"{full_haplotype}_vs_Hap1_defining_snps.csv"
    out.to_csv(out_path, index=False)
    return out


def _hap_allele(haplotypes: pd.DataFrame, variant_info: pd.DataFrame, hap_name: str, marker_id: str) -> str:
    hit = haplotypes.loc[haplotypes["Hap_Name"].eq(hap_name)]
    if hit.empty or marker_id not in set(variant_info["marker_id"].astype(str)):
        return ""
    index = int(variant_info.index[variant_info["marker_id"].astype(str).eq(marker_id)][0])
    alleles = str(hit.iloc[0]["Haplotype_Seq"]).split("|")
    return alleles[index] if index < len(alleles) else ""


def build_literature_coverage_table(
    db_dir: Path,
    defining_snps: pd.DataFrame,
    full_haplotype: str,
) -> pd.DataFrame:
    variant_info = _read_required_csv(db_dir / "variant_info.csv")
    haplotypes = _read_required_csv(db_dir / "haplotype_data.csv")
    marker_ids = set(variant_info["marker_id"].astype(str))
    defining_ids = set(defining_snps["marker_id"].astype(str)) if not defining_snps.empty else set()
    rows = []
    for item in PUBLISHED_VARIANTS:
        offset = item["promoter_offset_from_atg"]
        position = ""
        marker_id = ""
        if offset is not None:
            position = map_promoter_offset_to_position(GENE_START, CDS_START_OFFSET_1BASED, int(offset))
            marker_id = f"chr6A_{position}"
        covered = bool(marker_id and marker_id in marker_ids)
        row = dict(item)
        row.update({
            "refseqv1_position": position,
            "marker_id": marker_id,
            "covered_in_current_snp_markers": covered,
            "hap8_defining_snp": bool(marker_id and marker_id in defining_ids),
            "Hap1_allele_current": _hap_allele(haplotypes, variant_info, "Hap1", marker_id) if covered else "",
            f"{full_haplotype}_allele_current": _hap_allele(haplotypes, variant_info, full_haplotype, marker_id) if covered else "",
        })
        rows.append(row)
    return pd.DataFrame(rows)


def background_alleles_from_defining_snps(
    defining_snps: pd.DataFrame,
    literature_alleles: Dict[str, str],
) -> Dict[str, str]:
    out: Dict[str, str] = {}
    for row in defining_snps.to_dict("records"):
        marker_id = str(row["marker_id"])
        if marker_id in literature_alleles:
            continue
        out[marker_id] = str(row.get("Hap8_allele", ""))
    return out


def summarize_split_groups(df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for group, sub in df.groupby("Hap8SplitGroup", dropna=False):
        row = {"Hap8SplitGroup": group, "n": int(len(sub))}
        row["haplotypes"] = ";".join(f"{k}:{v}" for k, v in sub["Hap_Name"].value_counts().sort_index().items())
        for trait in TRAIT_COLUMNS:
            values = pd.to_numeric(sub[trait], errors="coerce").dropna()
            row[f"{trait}_n"] = int(len(values))
            row[f"{trait}_mean"] = float(values.mean()) if len(values) else math.nan
            row[f"{trait}_sd"] = float(values.std(ddof=1)) if len(values) > 1 else math.nan
        rows.append(row)
    order = {"Published+Background": 0, "PublishedOnly": 1, "BackgroundOnly": 2, "Neither": 3}
    return pd.DataFrame(rows).sort_values(
        by="Hap8SplitGroup",
        key=lambda s: s.map(order).fillna(99),
    )


def build_trait_tests(df: pd.DataFrame) -> pd.DataFrame:
    tests = []
    binary_specs = [
        ("PublishedPair", "PublishedPair", "NonPublishedPair"),
        ("Hap8Background", "Hap8Background", "NonHap8Background"),
        ("FullHap8", "FullHap8", "NonHap8"),
    ]
    for group_col, case_label, control_label in binary_specs:
        for trait in TRAIT_COLUMNS:
            tests.append(summarize_binary_trait(df, group_col, trait, case_label, control_label))

    for trait in TRAIT_COLUMNS:
        subset = df[df["Hap8SplitGroup"].isin(["Published+Background", "PublishedOnly"])].copy()
        subset["_pub_bg"] = subset["Hap8SplitGroup"].eq("Published+Background")
        tests.append(
            summarize_binary_trait(
                subset,
                "_pub_bg",
                trait,
                "Published+Background",
                "PublishedOnly",
            )
        )
    return pd.DataFrame(tests)


def build_country_counts(df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for (split_group, country), sub in df.groupby(["Hap8SplitGroup", "Country"], dropna=False):
        rows.append({
            "Hap8SplitGroup": split_group,
            "Country": country,
            "n": int(len(sub)),
            "haplotypes": ";".join(f"{k}:{v}" for k, v in sub["Hap_Name"].value_counts().sort_index().items()),
            "TGW_mean": float(pd.to_numeric(sub["TGW_mean"], errors="coerce").mean()),
        })
    return pd.DataFrame(rows).sort_values(["Hap8SplitGroup", "n", "Country"], ascending=[True, False, True])


def build_country_enrichment(df: pd.DataFrame, group_col: str, group_label: str) -> pd.DataFrame:
    rows = []
    countries = sorted(country for country in df["Country"].dropna().unique() if country)
    group_mask = df[group_col].astype(bool)
    for country in countries:
        country_mask = df["Country"].eq(country)
        table = np.array([
            [int((group_mask & country_mask).sum()), int((group_mask & ~country_mask).sum())],
            [int((~group_mask & country_mask).sum()), int((~group_mask & ~country_mask).sum())],
        ])
        if table.sum() == 0 or table[0].sum() == 0 or table[:, 0].sum() == 0:
            odds_ratio = math.nan
            pvalue = math.nan
        else:
            test = stats.fisher_exact(table)
            odds_ratio = float(test.statistic)
            pvalue = float(test.pvalue)
        rows.append({
            "group": group_label,
            "Country": country,
            "case_in_country": int(table[0, 0]),
            "case_total": int(group_mask.sum()),
            "control_in_country": int(table[1, 0]),
            "control_total": int((~group_mask).sum()),
            "odds_ratio": odds_ratio,
            "fisher_pvalue": pvalue,
        })
    out = pd.DataFrame(rows)
    out["fdr_bh"] = _bh_fdr(out["fisher_pvalue"].tolist())
    return out.sort_values(["fdr_bh", "fisher_pvalue", "Country"], na_position="last")


def build_country_adjusted_tests(df: pd.DataFrame) -> pd.DataFrame:
    adjusted = df.copy()
    rows = []
    for trait in TRAIT_COLUMNS:
        values = pd.to_numeric(adjusted[trait], errors="coerce")
        country_means = values.groupby(adjusted["Country"]).transform("mean")
        resid_col = f"{trait}_country_residual"
        adjusted[resid_col] = values - country_means
        for group_col, case_label, control_label in [
            ("PublishedPair", "PublishedPair", "NonPublishedPair"),
            ("Hap8Background", "Hap8Background", "NonHap8Background"),
            ("FullHap8", "FullHap8", "NonHap8"),
        ]:
            row = summarize_binary_trait(adjusted, group_col, resid_col, case_label, control_label)
            row["trait"] = trait
            row["adjustment"] = "within_country_centered"
            rows.append(row)
    return pd.DataFrame(rows)


def build_haplotype_literature_status(df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for hap_name, sub in df.groupby("Hap_Name"):
        row = {
            "Hap_Name": hap_name,
            "n": int(len(sub)),
            "PublishedPair_n": int(sub["PublishedPair"].sum()),
            "Hap8Background_n": int(sub["Hap8Background"].sum()),
            "FullHap8_n": int(sub["FullHap8"].sum()),
            "TGW_mean": float(pd.to_numeric(sub["TGW_mean"], errors="coerce").mean()),
            "countries": ";".join(f"{k}:{v}" for k, v in sub["Country"].value_counts().sort_index().items()),
        }
        rows.append(row)
    return pd.DataFrame(rows).sort_values(["PublishedPair_n", "Hap8Background_n", "n"], ascending=False)


def build_indel_sv_status(gw2_data_dir: Path) -> pd.DataFrame:
    files = list(gw2_data_dir.rglob("*")) if gw2_data_dir.exists() else []
    file_names = [path.name for path in files if path.is_file()]
    local_indel = [path for path in files if path.is_file() and "INDEL" in path.name.upper()]
    local_sv = [
        path for path in files
        if path.is_file() and any(token in path.name.upper() for token in ["SV", "PAV", "CNV"])
    ]
    return pd.DataFrame([
        {
            "item": "local_chr6A_indel_or_sv_vcf",
            "status": "present" if local_indel or local_sv else "missing",
            "local_match_count": len(local_indel) + len(local_sv),
            "details": ";".join(str(path) for path in [*local_indel, *local_sv]) if (local_indel or local_sv) else "",
        },
        {
            "item": "earlham_chr6A_indel_vcf",
            "status": "available_remote_large_not_downloaded",
            "local_match_count": int(any(name == Path(CHR6A_INDEL_URL).name for name in file_names)),
            "details": CHR6A_INDEL_URL,
        },
        {
            "item": "earlham_chr6A_indel_csi",
            "status": "available_remote_small_index",
            "local_match_count": int(any(name == Path(CHR6A_INDEL_INDEX_URL).name for name in file_names)),
            "details": CHR6A_INDEL_INDEX_URL,
        },
        {
            "item": "local_extraction_tools",
            "status": "requires_bcftools_or_pysam_for_indexed_region_extraction",
            "local_match_count": 0,
            "details": "Current Windows environment lacks pysam/bcftools in the checked run.",
        },
    ])


def write_report(
    output_path: Path,
    coverage: pd.DataFrame,
    split_summary: pd.DataFrame,
    trait_tests: pd.DataFrame,
    country_adjusted: pd.DataFrame,
    indel_status: pd.DataFrame,
) -> None:
    def _fmt(value: object, digits: int = 4) -> str:
        if pd.isna(value):
            return "NA"
        if isinstance(value, (float, np.floating)):
            return f"{float(value):.{digits}g}"
        return str(value)

    coverage_pair = coverage[coverage["tested_group"].eq("su_qin_two_snp_pair")]
    pair_covered = bool(coverage_pair["covered_in_current_snp_markers"].all()) if not coverage_pair.empty else False
    causal_494 = coverage[coverage["variant_name"].eq("SNP-494")].iloc[0]
    g2373a = coverage[coverage["variant_name"].eq("G2373A")].iloc[0]

    key_tests = trait_tests[
        (trait_tests["trait"].eq("TGW_mean"))
        & (trait_tests["group_col"].isin(["PublishedPair", "Hap8Background", "FullHap8", "_pub_bg"]))
    ].copy()
    residual_tests = country_adjusted[
        (country_adjusted["trait"].eq("TGW_mean"))
        & (country_adjusted["group_col"].isin(["PublishedPair", "Hap8Background", "FullHap8"]))
    ].copy()

    lines = [
        "# TaGW2-A1 Hap8 functional-marker diagnosis",
        "",
        "## Interpretation",
        "",
        "- Hap8 is not a one-to-one proxy for the published TaGW2-A1 functional haplotype.",
        "- The two local promoter SNPs carried by Hap8 are covered, but the same two-SNP state also occurs outside Hap8.",
        "- The later Jaiswal four-SNP promoter Hap5 cannot be confirmed because SNP-494 is absent from the current SNP-only database.",
        "- The G2373A splice mutation and coding indel/SV variants are not tested by the current SNP-only A1 database.",
        "",
        "## Marker Coverage",
        "",
        f"- Su/Qin two-SNP pair covered in current SNP markers: {pair_covered}.",
        f"- Jaiswal causal SNP-494 covered: {bool(causal_494['covered_in_current_snp_markers'])}.",
        f"- G2373A tested: {bool(g2373a['covered_in_current_snp_markers'])}.",
        "",
        "## Group Summary",
        "",
    ]
    for row in split_summary.to_dict("records"):
        lines.append(
            f"- {row['Hap8SplitGroup']}: n={row['n']}, "
            f"TGW_mean={_fmt(row.get('TGW_mean_mean'))}, haplotypes={row.get('haplotypes', '')}"
        )
    lines.extend(["", "## TGW Mean Tests", ""])
    for row in key_tests.to_dict("records"):
        lines.append(
            f"- {row['case_label']} vs {row['control_label']}: "
            f"n={row['case_n']}/{row['control_n']}, "
            f"effect={_fmt(row['effect_case_minus_control'])}, "
            f"p={_fmt(row['welch_pvalue'])}"
        )
    lines.extend(["", "## Country-Centered TGW Tests", ""])
    for row in residual_tests.to_dict("records"):
        lines.append(
            f"- {row['case_label']} vs {row['control_label']}: "
            f"effect={_fmt(row['effect_case_minus_control'])}, "
            f"p={_fmt(row['welch_pvalue'])}"
        )
    lines.extend(["", "## INDEL/SV Status", ""])
    for row in indel_status.to_dict("records"):
        lines.append(f"- {row['item']}: {row['status']}; {row['details']}")

    output_path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_indel_sv_download_instructions(output_path: Path) -> None:
    output_path.write_text(
        "\n".join([
            "# TaGW2-A1 chr6A INDEL/SV supplement instructions",
            "",
            "The current TaGW2-A1 database was built from the WatSeq chr6A SNP VCF only.",
            "To supplement INDEL calls, download the matching chr6A INDEL VCF and CSI index.",
            "",
            "PowerShell download commands:",
            "",
            "```powershell",
            "New-Item -ItemType Directory -Force -Path 'D:/Desktop/data/GW2' | Out-Null",
            f"Invoke-WebRequest -Uri '{CHR6A_INDEL_URL}' -OutFile 'D:/Desktop/data/GW2/{Path(CHR6A_INDEL_URL).name}'",
            f"Invoke-WebRequest -Uri '{CHR6A_INDEL_INDEX_URL}' -OutFile 'D:/Desktop/data/GW2/{Path(CHR6A_INDEL_INDEX_URL).name}'",
            "```",
            "",
            "Current checked remote metadata:",
            "",
            "- chr6A INDEL VCF Content-Length: 8519486608 bytes, about 7.94 GiB.",
            "- chr6A INDEL CSI Content-Length: 462684 bytes.",
            "- A test `curl -I -r 0-1023` returned HTTP 200 rather than HTTP 206, so simple range download is not enough for remote subsetting on this machine.",
            "- This Windows environment currently lacks `pysam`, `bcftools`, and `tabix`; indexed regional extraction should be run on an environment with one of those tools.",
            "",
            "Suggested extraction after download on Linux/HPC with bcftools:",
            "",
            "```bash",
            "bcftools view -r 6A:237732651-237760058 \\",
            f"  {Path(CHR6A_INDEL_URL).name} \\",
            "  -Oz -o TaGW2-A1.chr6A.indel.region.vcf.gz",
            "tabix -p vcf TaGW2-A1.chr6A.indel.region.vcf.gz",
            "```",
            "",
            "After region extraction, the preparation script should be extended to merge SNP and INDEL markers for the same SampleID set before rerunning the Hap8 diagnosis.",
        ]) + "\n",
        encoding="utf-8",
    )


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = build_arg_parser().parse_args(argv)
    db_dir = Path(args.db_dir)
    result_dir = Path(args.result_dir)
    marker_matrix = Path(args.marker_matrix)
    phenotype_xlsx = Path(args.phenotype_xlsx)
    gw2_data_dir = Path(args.gw2_data_dir)
    full_haplotype = str(args.full_haplotype)
    result_dir.mkdir(parents=True, exist_ok=True)

    phenotype = _read_required_csv(db_dir / "phenotype_data.csv")
    marker_df = pd.read_csv(marker_matrix, sep="\t")
    country = load_watkins_country_table(phenotype_xlsx)
    defining_snps = build_hap8_defining_snps(db_dir, result_dir, full_haplotype)
    coverage = build_literature_coverage_table(db_dir, defining_snps, full_haplotype)

    background_alleles = background_alleles_from_defining_snps(defining_snps, SU_QIN_PAIR)
    required_marker_cols = sorted(set(SU_QIN_PAIR) | set(background_alleles))
    merged = phenotype.merge(marker_df[["SampleID", *required_marker_cols]], on="SampleID", how="left")
    merged = merged.merge(country, on="SampleID", how="left")
    merged["Country"] = merged["Country"].fillna("Unknown")
    classified = classify_literature_background_groups(
        merged,
        literature_alleles=SU_QIN_PAIR,
        background_alleles=background_alleles,
        full_haplotype_name=full_haplotype,
    )

    split_summary = summarize_split_groups(classified)
    trait_tests = build_trait_tests(classified)
    country_counts = build_country_counts(classified)
    enrichment = pd.concat([
        build_country_enrichment(classified, "PublishedPair", "PublishedPair"),
        build_country_enrichment(classified, "Hap8Background", "Hap8Background"),
        build_country_enrichment(classified, "FullHap8", "FullHap8"),
    ], ignore_index=True)
    country_adjusted = build_country_adjusted_tests(classified)
    hap_status = build_haplotype_literature_status(classified)
    indel_status = build_indel_sv_status(gw2_data_dir)

    outputs = {
        "literature_variant_coverage": result_dir / "tagw2_hap8_literature_variant_coverage.csv",
        "classified_samples": result_dir / "tagw2_hap8_classified_samples.csv",
        "split_group_summary": result_dir / "tagw2_hap8_split_group_summary.csv",
        "trait_tests": result_dir / "tagw2_hap8_group_trait_tests.csv",
        "country_counts": result_dir / "tagw2_hap8_country_counts.csv",
        "country_enrichment": result_dir / "tagw2_hap8_country_enrichment.csv",
        "country_adjusted_trait_tests": result_dir / "tagw2_hap8_country_adjusted_trait_tests.csv",
        "haplotype_literature_status": result_dir / "tagw2_hap8_haplotype_literature_status.csv",
        "indel_sv_status": result_dir / "tagw2_hap8_indel_sv_status.csv",
        "indel_sv_download_instructions": result_dir / "tagw2_hap8_indel_sv_download_instructions.md",
        "report": result_dir / "tagw2_hap8_functional_report.md",
    }

    coverage.to_csv(outputs["literature_variant_coverage"], index=False)
    classified.to_csv(outputs["classified_samples"], index=False)
    split_summary.to_csv(outputs["split_group_summary"], index=False)
    trait_tests.to_csv(outputs["trait_tests"], index=False)
    country_counts.to_csv(outputs["country_counts"], index=False)
    enrichment.to_csv(outputs["country_enrichment"], index=False)
    country_adjusted.to_csv(outputs["country_adjusted_trait_tests"], index=False)
    hap_status.to_csv(outputs["haplotype_literature_status"], index=False)
    indel_status.to_csv(outputs["indel_sv_status"], index=False)
    write_report(outputs["report"], coverage, split_summary, trait_tests, country_adjusted, indel_status)
    write_indel_sv_download_instructions(outputs["indel_sv_download_instructions"])

    print("[INFO] Hap8 defining SNPs:", len(defining_snps))
    print("[INFO] Published pair samples:", int(classified["PublishedPair"].sum()))
    print("[INFO] Full Hap8 samples:", int(classified["FullHap8"].sum()))
    print("[INFO] Hap8 background markers:", len(background_alleles))
    print("[INFO] Wrote outputs:")
    for label, path in outputs.items():
        print(f"  - {label}: {path}")
    print("[INFO] Summary JSON:")
    print(json.dumps({
        "published_pair_samples": int(classified["PublishedPair"].sum()),
        "full_hap8_samples": int(classified["FullHap8"].sum()),
        "published_pair_haplotypes": classified.loc[
            classified["PublishedPair"], "Hap_Name"
        ].value_counts().sort_index().to_dict(),
        "hap8_is_one_to_one_with_published_pair": bool(
            classified["PublishedPair"].equals(classified["FullHap8"])
        ),
    }, indent=2, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
