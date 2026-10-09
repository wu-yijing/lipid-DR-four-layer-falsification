#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Regenerate ``output_tables/PROVENANCE.md`` from ``tools/artefact_map.json``.

Why this exists
---------------
Most committed aggregate tables under ``output_tables/`` are **not** written under
those names by any script. The scripts write phase-named files into ``phase*_out/``
directories in the author's analysis workspace, and the ``layer*`` / ``MR_*`` names
in this repository are the result of a selection-and-rename step that used to be
undocumented. The README claimed that "the file names in ``output_tables/``
indicate which script produces which aggregate output", which was not mechanically
true.

``tools/artefact_map.json`` records, for every committed table, one of:

``direct``
    a committed script writes that exact file name into ``output_tables/``;
``collected``
    the file is byte-identical to a named output of a named script in the analysis
    workspace (the rename step is then explicit and verifiable);
``hand-assembled``
    a summary table built by hand from script outputs — the highest-risk category,
    and the reason ``tools/audit_repo_consistency.py`` exists.

This script renders that map into a human-readable register and re-checks it.

Usage
-----
    python tools/make_provenance.py            # write output_tables/PROVENANCE.md
    python tools/make_provenance.py --check    # verify the map covers every table (exit 1 otherwise)
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import paths  # noqa: E402

MAP_PATH = Path(__file__).resolve().parent / "artefact_map.json"
OUT_PATH = None  # resolved at run time (output_tables/PROVENANCE.md)

# --------------------------------------------------------------------------
# annotation overlay: analysis variant, supersession, what the paper uses it for
# variant: primary | sensitivity | legacy
# --------------------------------------------------------------------------
VARIANT = {
    # layer 2, primary analysis (all screening records)
    "layer2_cohort_repeat_screening.csv": ("primary", "repeat-screening accounting"),
    "layer2_complete_case_performance.csv": ("primary", "Table 3 performance block"),
    "layer2_complete_case_coefficients_VIF.csv": ("primary", "Table 3 coefficients + VIF"),
    "layer2_MI_pooled_coefficients_3seeds.csv": ("primary", "Supplementary Table S7"),
    "layer2_missingness_rates.csv": ("primary", "Methods 2.4 missingness rates"),
    "layer2_missingness_mechanism.csv": ("primary", "Methods 2.4 missingness mechanism"),
    "layer2_fraction_missing_information.csv": ("primary", "Supplementary Figure S1 / Table S6"),
    "layer2_imputation_observed_vs_imputed.csv": ("primary", "Supplementary Figure S1 / Table S6"),
    "layer2_selection_model.csv": ("primary", "Section 3.2 completeness model (OR 1.44)"),
    "layer2_delta_auc_bootstrap.csv": ("primary", "Section 3.2 bootstrap difference"),
    "layer2_three_row_AUC_decomposition.csv": ("primary", "Section 3.2 three-row decomposition"),
    "layer2_pmsampsize_grid.csv": ("primary", "Methods 2.4 / limitation 6 sample size"),
    "layer2_pmsampsize_grid_fine.csv": ("primary", "limitation 6 grid (C = 0.56-0.70)"),
    "layer2_severity_coding_table.csv": ("primary", "Methods 2.4 severity coding"),
    "layer2_outcome_rebuild_patient_grades.csv": ("primary", "Methods 2.4 grade distribution"),
    "layer2_outcome_rebuild_exclusions.csv": ("primary", "Methods 2.4 / limitation 12 exclusions"),
    "layer2_outcome_rebuild_eye_categories.csv": ("primary", "Methods 2.4 eye-level categories"),
    "layer2_outcome_crosstab_full.csv": ("primary", "Methods 2.4 / Supplementary Table S11"),
    "layer2_missingness_patterns.csv": ("primary", "Supplementary Table S10"),
    "layer2_parsimonious_models.csv": ("primary", "Supplementary Table S9"),
    "layer2_outcome_audit_and_bounds.csv": ("primary", "Supplementary Table S13"),
    "layer2_independent_reextraction_summary.csv": ("primary", "Supplementary Table S14a"),
    "layer2_second_grader_sampling_summary.json": ("primary", "limitation 12 re-abstraction plan"),
    "layer2_text_mode_5level_matrix_r2r3_adj.csv": ("primary", "Supplementary Table S14a (adjudicated)"),
    "layer2_text_mode_agreement_r2r3_adj.json": ("primary", "Supplementary Table S14a"),
    # layer 3 / 4, primary analysis
    "layer3_external_validation_AUC.csv": ("primary", "Table 4 core / std / LASSO rows"),
    "layer3_weighted_sensitivity_MEC.csv": ("primary", "Table 4 MEC-weighted columns"),
    "layer3_recomputed_AUC_CI_B2000.csv": ("primary", "Table 4 confidence intervals"),
    "layer3_chinese_coefficients_by_cycle.csv": ("primary", "Section 3.3 temporal subsets"),
    "layer3_MI_pooled_export_sensitivity.csv": ("primary", "Section 3.3 MI-pooled export (0.517 / 0.518)"),
    "layer4_three_row_AUC.csv": ("primary", "Section 3.4 / Figure 4 three-row comparison"),
    "layer4_weighted_full_vs_lipid.csv": ("primary", "Section 3.4 weighted re-check + Brier"),
    "layer4_weighted_calibration.csv": ("primary", "Section 3.4 weighted calibration"),
    "layer4_full_reference_model.csv": ("primary", "Table 5"),
    "layer4_temporal_validation.csv": ("primary", "Section 3.4 temporal validation"),
    "layer34_dedup_transport_and_benchmark.csv": ("primary", "Table 1 / README transport + benchmark summary"),
    # layer 2 / 3 / 4, repeat-screening sensitivity (one record per patient)
    "layer2_dedup_cohort_and_outcome.csv": ("sensitivity", "limitation 12 sensitivity cohort"),
    "layer2_dedup_complete_case_performance.csv": ("sensitivity", "limitation 12"),
    "layer2_dedup_complete_case_coefficients_VIF.csv": ("sensitivity", "limitation 12"),
    "layer2_dedup_delta_auc_bootstrap.csv": ("sensitivity", "limitation 12 (+0.005)"),
    "layer2_dedup_discrimination_and_missingness.csv": ("sensitivity", "limitation 12 summary"),
    "layer2_dedup_missingness_rates.csv": ("sensitivity", "limitation 12"),
    "layer2_dedup_pmsampsize_grid.csv": ("sensitivity", "limitation 12"),
    "layer2_dedup_selection_model.csv": ("sensitivity", "limitation 12"),
    "layer2_dedup_screening_selection_rule.csv": ("sensitivity", "limitation 12 selection rule"),
    "layer2_dedup_three_row_AUC_decomposition.csv": ("sensitivity", "limitation 12"),
    "layer3_dedup_external_validation_AUC.csv": ("sensitivity", "limitation 12 (core 0.490 / std 0.485)"),
    "layer3_dedup_weighted_sensitivity_MEC.csv": ("sensitivity", "limitation 12 (slopes -1.65 / -3.19)"),
    "layer4_dedup_three_row_AUC.csv": ("sensitivity", "limitation 12 (0.218 preserved)"),
    "layer4_dedup_weighted_calibration.csv": ("sensitivity", "limitation 12"),
    "layer4_dedup_weighted_full_vs_lipid.csv": ("sensitivity", "limitation 12"),
    # layer 1 and cis
    "MR_two_sample_all_methods.csv": ("primary", "Table 2 (all methods)"),
    "MR_Egger_intercepts.csv": ("primary", "Table 2 note / Supplementary Table S3"),
    "MR_heterogeneity_CochranQ.csv": ("primary", "Table 2 note / Supplementary Table S3"),
    "MR_PRESSO_10k_summary.csv": ("primary", "Table 2 MR-PRESSO row"),
    "MVMR_core_TG_LDL_HDL.csv": ("primary", "Table 2 MVMR block"),
    "MR_exmore_sensitivity.csv": ("primary", "Section 3.1 / Supplementary Table S12"),
    "cis_drug_target_MR_six_targets.csv": ("primary", "Figure 2B / Supplementary Table S4"),
    "cis_instrument_strength_by_target.csv": ("primary", "Methods 2.3 instrument strength"),
    "cis_marginal_F_TG_LDL_HDL.csv": ("primary", "Methods 2.3 conditional/marginal F"),
    "instrument_mean_F_by_exposure.csv": ("primary", "Methods 2.3 mean F by exposure"),
    "emerging_markers_NHANES.csv": ("primary", "Section 3.5 / Supplementary Table S1"),
    "emerging_markers_chinese_cohort.csv": ("primary", "Section 3.5 / Supplementary Table S1"),
    # added 2026-10-09 to close reporting gaps and repair the alignment defect
    "layer2_reporting_gaps.csv": ("primary", "Section 3.2 / Table 3 - intervals for AUC and IDI"),
    "layer2_design_condition_numbers.csv": ("primary", "Section 3.2 - basis of the condition number"),
    "layer2_parsimonious_markers_classified.csv": ("primary", "Supplementary Table S9 - explicit marker classes"),
    "layer2_second_grader_composition.csv": ("primary", "limitation 12 - composition of the review sample"),
    "layer2_pmsampsize_verification.csv": ("primary", "limitation 6 - verification of the quoted requirements"),
    "layer2_reporting_gap_inventory.csv": ("primary", "audit trail - statistics still without an interval"),
    "layer4_local_comparator_increment.csv": ("primary", "Section 3.4 - attributable increment with intervals"),
    "layer4_three_row_AUC_corrected.csv": ("primary", "Section 3.4 - corrected three-row comparison"),
}

# --------------------------------------------------------------------------
# superseded / conflicting artefacts. These are kept for provenance but must not
# be quoted. The audit asserts that these entries stay flagged.
# --------------------------------------------------------------------------
SUPERSEDED = {
    "layer4_increment_NRI_IDI.csv": (
        "full-model out-of-fold AUC 0.7354 and dAUC 0.2494, from the earlier phase-4 fold "
        "split; the current values are 0.7331 and 0.2470 in layer4_three_row_AUC.csv. The NRI "
        "and IDI it reports (0.410 / 0.152) are therefore computed on a different dAUC basis "
        "than the dAUC quoted alongside them. Either recompute NRI/IDI on the current split or "
        "state the provenance when citing 0.410 / 0.152."
    ),
    "layer2_outcome_rebuild_crosstab.csv": (
        "SUPERSEDED and actively misleading if read as counts. Its cells are the binary "
        "comparison 'source text states DR / does not' at an earlier pipeline stage, with "
        "analysed_as_non_case n = 2,207 rather than the current 1,691. The labels "
        "analysed_as_case = 271 + 4 and analysed_as_non_case = 1,752 + 455 were read in an "
        "earlier draft of the manuscript as '4 cases had an explicit diagnosis and 455 "
        "diagnosed records were analysed as non-cases'. That reading is WRONG. The current "
        "source of truth is layer2_outcome_crosstab_full.csv: 474 of the 475 records carrying "
        "the miscoded value are true DR cases, and none of the 275 records analysed as cases "
        "carries a gradeable DR diagnosis (50 are grade 0, 225 are not gradeable)."
    ),
}

# --------------------------------------------------------------------------
# hand-assembled summaries: state what they were assembled from
# --------------------------------------------------------------------------
# --------------------------------------------------------------------------
# files that are correct for some rows/columns and wrong for others. These stay in
# place because they are cited by the manuscript, so the defective part is named
# explicitly instead of the whole file being disowned.
# --------------------------------------------------------------------------
PARTIAL = {
    "layer4_three_row_AUC.csv": {
        "defective_rows": "local_lipid_only_std",
        "why": (
            "The `local_lipid_only_std` row reports n = 888 and AUC_unweighted = 0.4861. "
            "Both are artefacts of an alignment defect in the `oof()` helper of "
            "phase4e_local_v2.R: `model.matrix()` drops rows with missing predictors, so for "
            "the extended lipid model (terms tg_log and ldl_mmol, missing in 473 of the 888 "
            "reference records) the fold-wise replacement vector was shorter than "
            "`p[fold == k]` and R recycled it with a warning rather than failing. The "
            "quantity genuinely exists on 415 records and equals 0.4988 (95% CI "
            "0.4368-0.5575). The other three rows are unaffected and still reproduce exactly. "
            "Corrected values: layer4_three_row_AUC_corrected.csv, computed by "
            "4_external/phase10_local_comparator_increment.R, which asserts the reproduction "
            "of the three unaffected rows and asserts that the corrected std row is not 888."
        ),
        "manuscript_impact": (
            "Section 3.4's sentence about the extended variable set must drop 'AUC 0.486' and "
            "either quote 0.4988 on 415 records or remove the sentence. The conclusion - the "
            "lipid-only model performs no better - is unchanged."
        ),
    },
    "layer4_dedup_three_row_AUC.csv": {
        "defective_rows": "local_lipid_only_std",
        "why": (
            "The same alignment defect, in the repeat-screening sensitivity run "
            "(phase4e_local_v2_dedup.R). Not corrected here because no manuscript statement "
            "depends on this row; re-running the patched script regenerates it correctly."
        ),
        "manuscript_impact": "none (no manuscript statement cites this row)",
    },
}

HAND = {
    "layer2_cohort_repeat_screening.csv": "transcribed from the console output of 3_prediction/build_dedup_devset.py",
    "layer2_dedup_cohort_and_outcome.csv": "assembled by hand from the layer2_dedup_* outputs",
    "layer2_dedup_discrimination_and_missingness.csv":
        "assembled by hand: primary and repeat-screening values placed side by side",
    "layer2_dedup_screening_selection_rule.csv":
        "assembled by hand from 3_prediction/build_dedup_devset.py and the dedup outcome comparison",
    "layer2_outcome_rebuild_patient_grades.csv": "assembled by hand from severity_rebuild_audit.csv (3_prediction/rebuild_outcome.py)",
    "layer2_outcome_rebuild_exclusions.csv": "assembled by hand from severity_rebuild_audit.csv",
    "layer2_outcome_rebuild_eye_categories.csv": "assembled by hand from severity_rebuild_audit.csv",
    "layer2_outcome_rebuild_crosstab.csv": "assembled by hand at an earlier pipeline stage (see SUPERSEDED)",
    "layer2_pmsampsize_grid_fine.csv": "assembled by hand from the pmsampsize routine in 3_prediction/phase2e_tail2.R",
    "layer2_text_mode_5level_matrix_r2r3_adj.csv": "assembled by hand from 3_prediction/phase9c_text_mode.py (adjudicated matrix)",
    "layer34_dedup_transport_and_benchmark.csv":
        "assembled by hand from layer3_* and layer4_* outputs; every value in it is independently "
        "re-derived by tools/audit_repo_consistency.py, which is why it is safe to cite",
}


def load_map() -> dict:
    return json.loads(MAP_PATH.read_text(encoding="utf-8"))


def committed() -> list[str]:
    d = paths.out_dir()
    return sorted(p.name for p in d.iterdir() if p.suffix in (".csv", ".json"))


def render(m: dict) -> str:
    rows = sorted(m.items())
    by_kind: dict[str, list] = {"direct": [], "collected": [], "hand-assembled": [], "local-only": []}
    for name, e in rows:
        by_kind.setdefault(e["provenance"], []).append((name, e))

    L = []
    A = L.append
    A("# Aggregate result tables: provenance register")
    A("")
    A("Every file under `output_tables/` is listed below with the script that produces it, the")
    A("analysis variant it belongs to, and which part of the manuscript it supports. This")
    A("register exists because the `layer*` / `MR_*` file names are **not** the names the")
    A("scripts write: the scripts emit phase-named files into `phase*_out/` in the analysis")
    A("workspace, and this repository holds a documented selection of those outputs under")
    A("stable names. `tools/artefact_map.json` is the machine-readable form and")
    A("`tools/make_provenance.py` regenerates this file from it.")
    A("")
    A("**Run `tools/audit_repo_consistency.py` before quoting any number from this directory.**")
    A("The audit re-derives every quantity the manuscript reports and fails if a committed")
    A("table stops supporting it.")
    A("")
    A("## Provenance classes")
    A("")
    A("| class | n | meaning |")
    A("|---|---|---|")
    for k, label in (("direct", "a committed script writes this exact name into `output_tables/`"),
                     ("collected", "byte-identical to a named output of a named script (rename step made explicit)"),
                     ("hand-assembled", "built by hand from script outputs; values re-derived by the audit"),
                     ("local-only", "record-level file held locally only; listed in `.gitignore`, never committed")):
        A("| `%s` | %d | %s |" % (k, len(by_kind.get(k, [])), label))
    A("")

    A("## Superseded or conflicting artefacts (do not quote)")
    A("")
    if SUPERSEDED:
        for name, why in SUPERSEDED.items():
            A("### `%s`" % name)
            A("")
            A(why)
            A("")
    else:
        A("_none_")
        A("")

    A("## Rows that are known to be defective (partial)")
    A("")
    A("These files are cited by the manuscript and are correct except for the rows named")
    A("below. `tools/audit_repo_consistency.py` asserts that each of these flags is present,")
    A("so the defect cannot be forgotten or silently repaired out of the register.")
    A("")
    for name, d in PARTIAL.items():
        A("### `%s` — row `%s`" % (name, d["defective_rows"]))
        A("")
        A(d["why"])
        A("")
        A("**Manuscript impact:** %s" % d["manuscript_impact"])
        A("")

    A("## Hand-assembled summaries")
    A("")
    A("These have no producing script. They are the highest-risk class in a data-integrity")
    A("screen, so each is stated explicitly.")
    A("")
    A("| file | assembled from |")
    A("|---|---|")
    for name, _ in sorted(by_kind["hand-assembled"]):
        A("| `%s` | %s |" % (name, HAND.get(name, "**undocumented — please describe**")))
    A("")

    for kind, title in (("direct", "Written directly by a committed script"),
                        ("collected", "Collected and renamed from a script output"),
                        ("hand-assembled", "Hand-assembled summaries")):
        A("## %s" % title)
        A("")
        A("| file | variant | producing script / source | supports |")
        A("|---|---|---|---|")
        for name, e in by_kind.get(kind, []):
            var, feeds = VARIANT.get(name, ("primary", ""))
            if name in SUPERSEDED:
                var = "legacy (superseded)"
            if kind == "direct":
                src = ", ".join("`%s`" % s for s in e.get("written_by", []))
            elif kind == "collected":
                src = ", ".join("`%s`" % s for s in e.get("sources", []))
            else:
                src = "_(none)_"
            A("| `%s` | %s | %s | %s |" % (name, var, src, feeds))
        A("")

    A("## Held locally only (never committed)")
    A("")
    A("These carry record-level identifiers or source-report text and are excluded by")
    A("`.gitignore`. They are listed here so that the gap between the local workspace and the")
    A("repository is explicit rather than implied.")
    A("")
    A("| file | producing script |")
    A("|---|---|")
    for name, e in by_kind.get("local-only", []):
        A("| `%s` | %s |" % (name, ", ".join("`%s`" % s for s in e.get("written_by", []))))
    A("")

    L.append(_variant_summary())
    return "\n".join(L) + "\n"


def _variant_summary() -> str:
    L = []
    A = L.append
    A("## How to cite a variant in the manuscript")
    A("")
    A("| manuscript statement | correct file |")
    A("|---|---|")
    A("| Table 1 / layer-3 transport and layer-4 benchmark summary | `layer34_dedup_transport_and_benchmark.csv` |")
    A("| Table 4 core / std / LASSO unweighted AUC | `layer3_external_validation_AUC.csv` |")
    A("| Table 4 MEC-weighted AUC and calibration slopes | `layer3_weighted_sensitivity_MEC.csv` |")
    A("| Section 3.4 three-row comparison and the 0.218 / 0.029 decomposition | `layer4_three_row_AUC.csv` |")
    A("| Section 3.4 weighted re-check, Brier 0.176 / 0.203, slope -1.52 | `layer4_weighted_full_vs_lipid.csv`, `layer4_weighted_calibration.csv` |")
    A("| Section 3.2 complete-case AUC 0.598 / 0.579 and calibration 1.000 / 0.000 | `layer2_complete_case_performance.csv` |")
    A("| Section 3.2 multiply-imputed AUC 0.582 and the CC-MI comparison | `layer2_dedup_discrimination_and_missingness.csv` |")
    A("| Methods 2.4 grade distribution 1,497 / 49 / 386 / 15 / 24 | `layer2_outcome_rebuild_patient_grades.csv` |")
    A("| Methods 2.4 exclusions 511 = 386 + 53 + 72 | `layer2_outcome_rebuild_exclusions.csv` |")
    A("| Methods 2.4 / Table S11 record-level cross-tabulation | `layer2_outcome_crosstab_full.csv` |")
    A("")
    A("`_dedup` files are the **repeat-screening sensitivity analysis** (one record per patient)")
    A("and must never be quoted as if they were the primary analysis. Note that the primary and")
    A("`_dedup` layer-3 files carry near-identical confidence intervals but different point")
    A("estimates and very different calibration slopes; always cite the file name alongside the")
    A("number.")
    return "\n".join(L)


def check(m: dict) -> int:
    have = set(committed())
    known = set(m)
    unknown = sorted(have - known)
    print("artefact map entries: %d | files present: %d" % (len(known), len(have)))
    bad = 0
    if unknown:
        bad += len(unknown)
        print("\nPRESENT BUT NOT IN THE MAP (add them to tools/artefact_map.json):")
        for u in unknown:
            print("   ", u)
    stale = sorted(known - have)
    if stale:
        print("\nIN THE MAP BUT NOT PRESENT LOCALLY (fine if untracked/aggregate-only):")
        for s in stale:
            print("   ", s)
    for name, e in m.items():
        if e["provenance"] == "hand-assembled" and name not in HAND:
            print("\nHAND-ASSEMBLED BUT UNDESCRIBED:", name)
            bad += 1
    return 1 if bad else 0


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--check", action="store_true", help="verify the map covers every committed table")
    a = ap.parse_args(argv)
    m = load_map()
    if a.check:
        return check(m)
    out = paths.out_dir() / "PROVENANCE.md"
    out.write_text(render(m), encoding="utf-8")
    print("wrote", out)
    return check(m)


if __name__ == "__main__":
    sys.exit(main())
