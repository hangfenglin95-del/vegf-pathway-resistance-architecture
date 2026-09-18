#!/usr/bin/env python3
"""Rebuild Phase 10 validation summaries from bundled frozen results.

This script does not rerun differential expression, enrichment, or pathway
selection. It summarizes the 72 frozen-program results in Table S5 while
preserving GEO-series nesting.
"""

from __future__ import annotations

import csv
from collections import Counter, defaultdict
from pathlib import Path
from statistics import median


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "results" / "validation" / "TableS5_expanded_external_validation_statistics.tsv"
COHORT = ROOT / "config" / "external_validation_cohort.tsv"
OUT = ROOT / "outputs" / "validation"
PATHWAYS = [
    "HALLMARK_MTORC1_SIGNALING",
    "HALLMARK_UNFOLDED_PROTEIN_RESPONSE",
    "HALLMARK_INTERFERON_GAMMA_RESPONSE",
    "HALLMARK_INTERFERON_ALPHA_RESPONSE",
    "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
    "HALLMARK_INFLAMMATORY_RESPONSE",
]
PATHWAY_DISPLAY = {
    "HALLMARK_MTORC1_SIGNALING": "MTORC1 signaling",
    "HALLMARK_UNFOLDED_PROTEIN_RESPONSE": "Unfolded protein response",
    "HALLMARK_INTERFERON_GAMMA_RESPONSE": "Interferon-gamma response",
    "HALLMARK_INTERFERON_ALPHA_RESPONSE": "Interferon-alpha response",
    "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION": "Epithelial-mesenchymal transition",
    "HALLMARK_INFLAMMATORY_RESPONSE": "Inflammatory response",
}


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def direction(values: list[str]) -> str:
    signs = set(values)
    if signs == {"positive"}:
        return "positive"
    if signs == {"negative"}:
        return "negative"
    return "mixed"


def summarize(rows: list[dict[str, str]]) -> None:
    cohort = read_tsv(COHORT)
    dataset_order = list(dict.fromkeys(row["dataset"] for row in cohort))
    context_order = [row["context_id"] for row in cohort]
    tier_by_context = {row["context_id"]: row["evidence_tier"] for row in cohort}
    nested_by_context = {row["context_id"]: row["nested_context_group"] for row in cohort}

    by_path_dataset: dict[tuple[str, str], list[dict[str, str]]] = defaultdict(list)
    by_path_context: dict[tuple[str, str], dict[str, str]] = {}
    for row in rows:
        by_path_dataset[(row["pathway"], row["dataset"])].append(row)
        by_path_context[(row["pathway"], row["contrast"])] = row

    dataset_rows: list[dict[str, object]] = []
    program_rows: list[dict[str, object]] = []
    for pathway in PATHWAYS:
        statuses = []
        pathway_rows = [row for row in rows if row["pathway"] == pathway]
        for dataset in dataset_order:
            block = by_path_dataset[(pathway, dataset)]
            directions = [row["direction"] for row in block]
            status = direction(directions)
            nes = [float(row["NES"]) for row in block]
            fdr = [float(row["padj"]) for row in block]
            statuses.append(status)
            dataset_rows.append({
                "pathway": pathway,
                "pathway_display": PATHWAY_DISPLAY[pathway],
                "dataset": dataset,
                "n_contexts": len(block),
                "context_directions": ";".join(directions),
                "dataset_direction": status,
                "median_NES": format(median(nes), ".15g"),
                "min_NES": format(min(nes), ".15g"),
                "max_NES": format(max(nes), ".15g"),
                "n_FDR_lt_0.05": sum(value < 0.05 for value in fdr),
                "minimum_FDR": format(min(fdr), ".15g"),
                "evidence_tier": ";".join(dict.fromkeys(tier_by_context[row["contrast"]] for row in block)),
                "nested_context_group": nested_by_context[block[0]["contrast"]],
            })
        nes = [float(row["NES"]) for row in pathway_rows]
        program_rows.append({
            "pathway": pathway,
            "pathway_display": PATHWAY_DISPLAY[pathway],
            "n_independent_datasets": len(dataset_order),
            "dataset_positive": statuses.count("positive"),
            "dataset_negative": statuses.count("negative"),
            "dataset_mixed": statuses.count("mixed"),
            "n_contexts": len(pathway_rows),
            "context_positive": sum(row["direction"] == "positive" for row in pathway_rows),
            "context_negative": sum(row["direction"] == "negative" for row in pathway_rows),
            "median_context_NES": format(median(nes), ".15g"),
            "minimum_context_NES": format(min(nes), ".15g"),
            "maximum_context_NES": format(max(nes), ".15g"),
            "n_nominal_P_lt_0.05": sum(float(row["pval"]) < 0.05 for row in pathway_rows),
            "n_FDR_lt_0.05": sum(float(row["padj"]) < 0.05 for row in pathway_rows),
        })

    lodo_rows: list[dict[str, object]] = []
    for pathway in PATHWAYS:
        for omitted in dataset_order:
            kept = [row for row in dataset_rows if row["pathway"] == pathway and row["dataset"] != omitted]
            counts = Counter(row["dataset_direction"] for row in kept)
            lodo_rows.append({
                "pathway": pathway,
                "omitted_dataset": omitted,
                "remaining_datasets": len(kept),
                "dataset_positive": counts["positive"],
                "dataset_negative": counts["negative"],
                "dataset_mixed": counts["mixed"],
                "remaining_contexts": len([row for row in rows if row["pathway"] == pathway and row["dataset"] != omitted]),
                "context_positive": sum(row["direction"] == "positive" for row in rows if row["pathway"] == pathway and row["dataset"] != omitted),
                "context_negative": sum(row["direction"] == "negative" for row in rows if row["pathway"] == pathway and row["dataset"] != omitted),
                "median_context_NES": format(median(float(row["NES"]) for row in rows if row["pathway"] == pathway and row["dataset"] != omitted), ".15g"),
            })

    loco_rows: list[dict[str, object]] = []
    for pathway in PATHWAYS:
        for omitted in context_order:
            kept = [row for row in rows if row["pathway"] == pathway and row["contrast"] != omitted]
            loco_rows.append({
                "pathway": pathway,
                "omitted_context": omitted,
                "remaining_contexts": len(kept),
                "context_positive": sum(row["direction"] == "positive" for row in kept),
                "context_negative": sum(row["direction"] == "negative" for row in kept),
                "median_context_NES": format(median(float(row["NES"]) for row in kept), ".15g"),
            })

    tier_rows: list[dict[str, object]] = []
    for pathway in PATHWAYS:
        for tier_name, selector in (
            ("Tier A including caution", lambda value: value.startswith("Tier A")),
            ("Tier B", lambda value: value == "Tier B"),
        ):
            block = [row for row in rows if row["pathway"] == pathway and selector(tier_by_context[row["contrast"]])]
            datasets = list(dict.fromkeys(row["dataset"] for row in block))
            statuses = [direction([row["direction"] for row in block if row["dataset"] == dataset]) for dataset in datasets]
            counts = Counter(statuses)
            tier_rows.append({
                "pathway": pathway,
                "tier_group": tier_name,
                "n_datasets": len(datasets),
                "dataset_positive": counts["positive"],
                "dataset_negative": counts["negative"],
                "dataset_mixed": counts["mixed"],
                "n_contexts": len(block),
                "context_positive": sum(row["direction"] == "positive" for row in block),
                "context_negative": sum(row["direction"] == "negative" for row in block),
                "median_context_NES": format(median(float(row["NES"]) for row in block), ".15g"),
            })

    write_tsv(OUT / "PHASE10_VALIDATION_DATASET_SUMMARY.tsv", dataset_rows, list(dataset_rows[0]))
    write_tsv(OUT / "PHASE10_VALIDATION_PROGRAM_SUMMARY.tsv", program_rows, list(program_rows[0]))
    write_tsv(OUT / "PHASE10_VALIDATION_LEAVE_ONE_DATASET_OUT.tsv", lodo_rows, list(lodo_rows[0]))
    write_tsv(OUT / "PHASE10_VALIDATION_LEAVE_ONE_CONTEXT_OUT.tsv", loco_rows, list(loco_rows[0]))
    write_tsv(OUT / "PHASE10_VALIDATION_TIER_SUMMARY.tsv", tier_rows, list(tier_rows[0]))


def main() -> None:
    rows = [row for row in read_tsv(SOURCE) if row["final_evidence_status"] == "FINAL_EXTERNAL_VALIDATION"]
    assert len(rows) == 72, f"expected 72 formal results, found {len(rows)}"
    assert len({row["dataset"] for row in rows}) == 8
    assert len({row["contrast"] for row in rows}) == 12
    assert Counter(row["pathway"] for row in rows) == Counter({pathway: 12 for pathway in PATHWAYS})
    summarize(rows)
    print("EXPANDED_VALIDATION_SUMMARY=PASS")


if __name__ == "__main__":
    main()
