from __future__ import annotations

from pathlib import Path

import pandas as pd
import streamlit as st


ROOT = Path(__file__).resolve().parent

MUTATION_EVIDENCE_PATH = ROOT / "references" / "mutation_evidence_table.tsv"
VERIFIED_EVIDENCE_PATH = ROOT / "data" / "processed" / "verified_mutation_evidence.tsv"
SPECIES_SUMMARY_PATH = ROOT / "data" / "processed" / "mutation_summary_by_species.tsv"
HERBICIDE_SUMMARY_PATH = ROOT / "data" / "processed" / "mutation_summary_by_herbicide.tsv"
MODEL_REGISTRY_PATH = ROOT / "structures" / "model_registry.tsv"
MODEL_QUALITY_PATH = ROOT / "data" / "processed" / "model_quality_summary.tsv"
REFERENCE_SEQUENCES_PATH = ROOT / "sequences" / "final_reference_sequences.tsv"
MUTATION_MAPPING_PATH = ROOT / "sequences" / "mutation_reference_mapping.tsv"


@st.cache_data
def load_tsv(path: Path) -> pd.DataFrame:
    if not path.exists():
        return pd.DataFrame()
    return pd.read_csv(path, sep="\t", dtype=str).fillna("")


def multiselect_filter(
    frame: pd.DataFrame, label: str, column: str, key: str
) -> pd.DataFrame:
    if column not in frame.columns:
        return frame
    values = sorted(value for value in frame[column].dropna().unique() if str(value).strip())
    selected = st.multiselect(label, values, key=key)
    if selected:
        return frame[frame[column].isin(selected)]
    return frame


def metric_value(frame: pd.DataFrame, column: str) -> int:
    if frame.empty or column not in frame.columns:
        return 0
    return frame[column].replace("", pd.NA).dropna().nunique()


def overview_tab() -> None:
    verified = load_tsv(VERIFIED_EVIDENCE_PATH)
    registry = load_tsv(MODEL_REGISTRY_PATH)
    species_summary = load_tsv(SPECIES_SUMMARY_PATH)
    herbicide_summary = load_tsv(HERBICIDE_SUMMARY_PATH)
    references = load_tsv(REFERENCE_SEQUENCES_PATH)
    mappings = load_tsv(MUTATION_MAPPING_PATH)

    st.write(
        "Current prototype integrates verified literature evidence, PPO2/PPX2 "
        "mutation mapping, and ColabFold-predicted protein structures. Docking "
        "and herbicide-binding analysis are planned for the next phase."
    )

    completed_models = 0
    review_models = 0
    if not registry.empty and "model_status" in registry.columns:
        completed_models = int((registry["model_status"] == "completed").sum())
        review_models = int((registry["model_status"] == "needs_manual_check").sum())

    metric_columns = st.columns(5)
    metric_columns[0].metric("Verified mutation evidence rows", len(verified))
    metric_columns[1].metric("Species count", metric_value(verified, "species"))
    metric_columns[2].metric("Unique mutation count", metric_value(verified, "mutation"))
    metric_columns[3].metric("Completed ColabFold models", completed_models)
    metric_columns[4].metric("Models needing manual check", review_models)

    with st.expander("Curated data inputs"):
        st.write(
            {
                "species summary rows": len(species_summary),
                "herbicide summary rows": len(herbicide_summary),
                "reference sequence rows": len(references),
                "mutation mapping rows": len(mappings),
            }
        )


def mutation_evidence_tab() -> None:
    evidence = load_tsv(MUTATION_EVIDENCE_PATH)
    if evidence.empty:
        st.warning("Mutation evidence table is not available.")
        return

    filtered = evidence.copy()
    filter_columns = st.columns(4)
    with filter_columns[0]:
        filtered = multiselect_filter(filtered, "Species", "species", "species_filter")
    with filter_columns[1]:
        filtered = multiselect_filter(filtered, "Mutation", "mutation", "mutation_filter")
    with filter_columns[2]:
        filtered = multiselect_filter(
            filtered,
            "Specific herbicides tested",
            "specific_herbicides_tested",
            "herbicide_filter",
        )
    with filter_columns[3]:
        filtered = multiselect_filter(
            filtered,
            "Verification status",
            "verification_status",
            "verification_filter",
        )

    st.dataframe(filtered, use_container_width=True, hide_index=True)

    st.subheader("Summary Counts")
    count_columns = st.columns(2)
    with count_columns[0]:
        if "species" in filtered.columns:
            st.write("By species")
            st.dataframe(
                filtered.groupby("species", dropna=False)
                .size()
                .reset_index(name="evidence_rows"),
                use_container_width=True,
                hide_index=True,
            )
    with count_columns[1]:
        if "mutation" in filtered.columns:
            st.write("By mutation")
            st.dataframe(
                filtered.groupby("mutation", dropna=False)
                .size()
                .reset_index(name="evidence_rows"),
                use_container_width=True,
                hide_index=True,
            )


def protein_models_tab() -> None:
    registry = load_tsv(MODEL_REGISTRY_PATH)
    quality = load_tsv(MODEL_QUALITY_PATH)
    model_data = quality if not quality.empty else registry

    st.write(
        "Predicted protein structures are available, but docking and "
        "binding-affinity interpretation are not yet implemented in this version."
    )

    if model_data.empty:
        st.warning("Model registry data is not available.")
        return

    display_columns = [
        "species",
        "mutation",
        "reference_accession",
        "mean_plddt",
        "pae_file_present",
        "model_status",
        "best_model_file_name",
    ]
    available_columns = [column for column in display_columns if column in model_data.columns]
    st.dataframe(model_data[available_columns], use_container_width=True, hide_index=True)

    if "model_status" in model_data.columns:
        review_rows = model_data[model_data["model_status"] == "needs_manual_check"]
        if not review_rows.empty:
            st.warning("Models needing manual check")
            st.dataframe(review_rows[available_columns], use_container_width=True, hide_index=True)


def future_work_tab() -> None:
    planned_items = [
        ("A. Ligand preparation", "1. Herbicide ligand table with PubChem IDs"),
        ("A. Ligand preparation", "2. Ligand 3D structure preparation"),
        ("B. Receptor/docking setup", "3. Receptor preparation from PPO models"),
        ("B. Receptor/docking setup", "4. Binding pocket/grid definition"),
        ("B. Receptor/docking setup", "5. FAD/cofactor decision"),
        ("C. Docking analysis", "6. Docking validation"),
        ("C. Docking analysis", "7. Docking against WT and mutant proteins"),
        ("C. Docking analysis", "8. Binding affinity summary"),
        ("C. Docking analysis", "9. Interpretation of cross-resistance predictions"),
        ("D. Dashboard integration", "10. Dashboard integration of docking results"),
    ]

    future_work = pd.DataFrame(
        [{"section": section, "planned_item": item, "status": "planned"} for section, item in planned_items]
    )
    for section, section_rows in future_work.groupby("section", sort=False):
        st.subheader(section)
        st.dataframe(
            section_rows[["planned_item", "status"]],
            use_container_width=True,
            hide_index=True,
        )


def main() -> None:
    st.set_page_config(
        page_title="PPO-Inhibitor Resistance Dashboard",
        layout="wide",
    )
    st.title("PPO-Inhibitor Resistance Dashboard for Amaranthus spp.")

    tabs = st.tabs(
        [
            "Overview",
            "Mutation Evidence",
            "Protein Models",
            "Future Work / Next Phase",
        ]
    )
    with tabs[0]:
        overview_tab()
    with tabs[1]:
        mutation_evidence_tab()
    with tabs[2]:
        protein_models_tab()
    with tabs[3]:
        future_work_tab()


if __name__ == "__main__":
    main()
