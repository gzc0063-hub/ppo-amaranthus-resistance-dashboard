from __future__ import annotations

from html import escape
from pathlib import Path

import pandas as pd
import py3Dmol
import streamlit as st
import streamlit.components.v1 as components


ROOT = Path(__file__).resolve().parent

MUTATION_EVIDENCE_PATH = ROOT / "references" / "mutation_evidence_table.tsv"
VERIFIED_EVIDENCE_PATH = ROOT / "data" / "processed" / "verified_mutation_evidence.tsv"
SPECIES_SUMMARY_PATH = ROOT / "data" / "processed" / "mutation_summary_by_species.tsv"
HERBICIDE_SUMMARY_PATH = ROOT / "data" / "processed" / "mutation_summary_by_herbicide.tsv"
MODEL_REGISTRY_PATH = ROOT / "structures" / "model_registry.tsv"
MODEL_QUALITY_PATH = ROOT / "data" / "processed" / "model_quality_summary.tsv"
REFERENCE_SEQUENCES_PATH = ROOT / "sequences" / "final_reference_sequences.tsv"
MUTATION_MAPPING_PATH = ROOT / "sequences" / "mutation_reference_mapping.tsv"


FUTURE_WORK_GROUPS = {
    "A. Ligand preparation": [
        (
            "Herbicide ligand table with PubChem IDs",
            "Curate ligand identifiers so later docking inputs are traceable.",
        ),
        (
            "Ligand 3D structure preparation",
            "Prepare standardized ligand structures before any docking runs.",
        ),
    ],
    "B. Receptor/docking setup": [
        (
            "Receptor preparation",
            "Prepare selected PPO2/PPX2 model structures for a documented docking workflow.",
        ),
        (
            "Binding-pocket/grid definition",
            "Define the docking search region before comparing WT and mutant models.",
        ),
        (
            "FAD/cofactor decision",
            "Document whether cofactors should be included in receptor preparation.",
        ),
    ],
    "C. Docking analysis": [
        (
            "Docking validation",
            "Validate the workflow before interpreting computational outputs.",
        ),
        (
            "Docking against WT and mutant proteins",
            "Compare predicted binding behavior only after setup is validated.",
        ),
        (
            "Binding-affinity summary",
            "Summarize docking outputs without treating scores as resistance proof.",
        ),
        (
            "Cross-resistance prediction interpretation",
            "Frame computational patterns as hypotheses, separate from verified evidence.",
        ),
    ],
    "D. Decision-support integration": [
        (
            "Dashboard integration of docking results",
            "Add docking summaries later once validated results exist.",
        )
    ],
}


COMPLETED_ITEMS = [
    "verified mutation evidence curation",
    "sequence/reference mapping",
    "ColabFold model generation",
    "3D model visualization",
]
IN_PROGRESS_ITEMS = [
    "final model completeness review",
    "dashboard refinement",
    "model interpretation workflow",
]
NOT_IMPLEMENTED_ITEMS = [
    "herbicide ligand table with PubChem IDs",
    "ligand 3D preparation",
    "receptor preparation",
    "binding-pocket/grid definition",
    "FAD/cofactor decision",
    "docking validation",
    "docking against WT and mutant proteins",
    "binding-affinity summary",
    "cross-resistance prediction interpretation",
    "docking integration into dashboard",
]


@st.cache_data
def load_tsv(path: Path) -> pd.DataFrame:
    if not path.exists():
        return pd.DataFrame()
    return pd.read_csv(path, sep="\t", dtype=str).fillna("")


def inject_css() -> None:
    st.markdown(
        """
        <style>
        :root {
            --bg: #f5f7f2;
            --ink: #17231d;
            --muted: #5e6b63;
            --line: #d9e2d7;
            --panel: #ffffff;
            --green: #1f6f4a;
            --green-dark: #174f38;
            --gold: #b8872d;
            --red: #b84a42;
            --blue: #3b6f87;
            --shadow: 0 18px 42px rgba(31, 55, 41, 0.10);
        }

        .stApp {
            background:
                radial-gradient(circle at top left, rgba(85, 137, 94, 0.18), transparent 34rem),
                linear-gradient(180deg, #f6f8f3 0%, #edf3ed 42%, #f7f9f5 100%);
            color: var(--ink);
        }

        [data-testid="stSidebar"] {
            background: linear-gradient(180deg, #173523 0%, #244f35 100%);
        }

        [data-testid="stSidebar"] * {
            color: #f6fbf3 !important;
        }

        [data-testid="stSidebar"] .stMarkdown a {
            color: #d7efbf !important;
            text-decoration: none;
        }

        .block-container {
            padding-top: 2rem;
            padding-bottom: 3rem;
            max-width: 1260px;
        }

        h1, h2, h3 {
            letter-spacing: 0;
        }

        div[data-testid="stMetric"] {
            background: #ffffff;
            border: 1px solid var(--line);
            border-radius: 10px;
            padding: 1rem 1.1rem;
            box-shadow: 0 10px 26px rgba(31, 55, 41, 0.08);
        }

        div[data-testid="stMetric"] label {
            color: var(--muted) !important;
            font-weight: 700;
        }

        .hero {
            display: grid;
            grid-template-columns: minmax(0, 1.45fr) minmax(280px, 0.8fr);
            gap: 1.25rem;
            align-items: stretch;
            padding: 2rem;
            border: 1px solid rgba(31, 111, 74, 0.18);
            border-radius: 16px;
            background:
                linear-gradient(135deg, rgba(255,255,255,0.96), rgba(233,244,231,0.94)),
                linear-gradient(135deg, rgba(31,111,74,0.08), rgba(184,135,45,0.10));
            box-shadow: var(--shadow);
            margin-bottom: 1.4rem;
        }

        .hero h1 {
            margin: 0 0 0.7rem 0;
            font-size: clamp(2.1rem, 5vw, 4.1rem);
            line-height: 1.02;
            color: #112019;
        }

        .hero p {
            color: #3e4d43;
            font-size: 1.08rem;
            line-height: 1.62;
            max-width: 760px;
            margin: 0.55rem 0 0 0;
        }

        .badge-row {
            display: flex;
            flex-wrap: wrap;
            gap: 0.55rem;
            margin-bottom: 1rem;
        }

        .badge {
            display: inline-flex;
            align-items: center;
            border-radius: 999px;
            padding: 0.34rem 0.72rem;
            font-size: 0.78rem;
            font-weight: 800;
            letter-spacing: 0.02em;
            border: 1px solid transparent;
        }

        .badge.green { color: #174f38; background: #dff0df; border-color: #bbd9bd; }
        .badge.gold { color: #6e4d13; background: #f5e6be; border-color: #e0c371; }
        .badge.red { color: #7f2d28; background: #f6d8d4; border-color: #e7aaa2; }
        .badge.blue { color: #244f63; background: #dcecf2; border-color: #accbd7; }

        .feature-panel {
            min-height: 100%;
            border-radius: 14px;
            background:
                linear-gradient(150deg, rgba(31,111,74,0.96), rgba(23,58,39,0.96)),
                repeating-linear-gradient(45deg, rgba(255,255,255,0.08) 0 1px, transparent 1px 14px);
            color: #f4fbf1;
            padding: 1.35rem;
            display: flex;
            flex-direction: column;
            justify-content: space-between;
            overflow: hidden;
            position: relative;
        }

        .feature-panel::after {
            content: "";
            position: absolute;
            width: 220px;
            height: 220px;
            right: -80px;
            bottom: -95px;
            border-radius: 50%;
            border: 34px solid rgba(215, 239, 191, 0.14);
        }

        .feature-panel h3 {
            color: #ffffff;
            margin: 0;
            font-size: 1.25rem;
        }

        .feature-panel p {
            color: #dcebd7;
            font-size: 0.95rem;
        }

        .section-card {
            background: rgba(255,255,255,0.93);
            border: 1px solid var(--line);
            border-radius: 14px;
            padding: 1.25rem;
            box-shadow: 0 12px 28px rgba(31,55,41,0.08);
            margin: 0.8rem 0 1rem 0;
        }

        .section-card:hover {
            box-shadow: 0 18px 38px rgba(31,55,41,0.12);
            transform: translateY(-1px);
            transition: all 160ms ease;
        }

        .section-kicker {
            color: var(--green-dark);
            font-size: 0.82rem;
            text-transform: uppercase;
            font-weight: 850;
            letter-spacing: 0.08em;
            margin-bottom: 0.25rem;
        }

        .section-title {
            font-size: 1.55rem;
            font-weight: 850;
            margin: 0 0 0.45rem 0;
            color: #15231b;
        }

        .muted {
            color: var(--muted);
            line-height: 1.58;
        }

        .info-grid {
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(220px, 1fr));
            gap: 0.9rem;
        }

        .info-card {
            background: #ffffff;
            border: 1px solid var(--line);
            border-left: 5px solid var(--green);
            border-radius: 12px;
            padding: 1rem;
            min-height: 132px;
            box-shadow: 0 8px 24px rgba(31,55,41,0.07);
        }

        .info-card h4 {
            margin: 0 0 0.35rem 0;
            color: #1d3528;
            font-size: 1rem;
        }

        .info-card p, .info-card li {
            color: #526057;
            line-height: 1.5;
            margin-bottom: 0.25rem;
        }

        .status-card {
            background: #ffffff;
            border: 1px solid var(--line);
            border-radius: 14px;
            padding: 1.1rem;
            min-height: 320px;
            box-shadow: 0 10px 24px rgba(31,55,41,0.08);
        }

        .status-card.completed { border-top: 6px solid var(--green); }
        .status-card.progress { border-top: 6px solid var(--gold); }
        .status-card.pending { border-top: 6px solid var(--red); }

        .status-card h3 {
            margin-top: 0;
            font-size: 1.22rem;
        }

        .status-card ul {
            padding-left: 1.1rem;
            margin-bottom: 0;
        }

        .status-card li {
            margin-bottom: 0.42rem;
            color: #435047;
        }

        .model-meta {
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(170px, 1fr));
            gap: 0.7rem;
        }

        .meta-cell {
            background: #f8fbf7;
            border: 1px solid var(--line);
            border-radius: 10px;
            padding: 0.85rem;
        }

        .meta-label {
            color: var(--muted);
            font-size: 0.78rem;
            font-weight: 800;
            text-transform: uppercase;
            letter-spacing: 0.04em;
            margin-bottom: 0.2rem;
        }

        .meta-value {
            color: #17231d;
            font-weight: 800;
            word-break: break-word;
        }

        .divider {
            height: 1px;
            background: linear-gradient(90deg, transparent, #cbd9cd, transparent);
            margin: 1.4rem 0;
        }

        @media (max-width: 900px) {
            .hero {
                grid-template-columns: 1fr;
                padding: 1.35rem;
            }
        }
        </style>
        """,
        unsafe_allow_html=True,
    )


def html_list(items: list[str]) -> str:
    return "".join(f"<li>{escape(item)}</li>" for item in items)


def section_intro(kicker: str, title: str, body: str) -> None:
    st.markdown(
        f"""
        <div class="section-card">
            <div class="section-kicker">{escape(kicker)}</div>
            <div class="section-title">{escape(title)}</div>
            <div class="muted">{escape(body)}</div>
        </div>
        """,
        unsafe_allow_html=True,
    )


def status_badge(label: str, color: str) -> str:
    return f'<span class="badge {color}">{escape(label)}</span>'


def render_hero() -> None:
    st.markdown(
        f"""
        <div class="hero">
            <div>
                <div class="badge-row">
                    {status_badge("Prototype", "green")}
                    {status_badge("In Progress", "gold")}
                    {status_badge("Docking Not Yet Implemented", "red")}
                </div>
                <h1>PPO-Inhibitor Resistance Dashboard for Amaranthus spp.</h1>
                <p>
                    A presentation-ready research dashboard for organizing PPO
                    target-site resistance evidence, sequence mapping, and
                    ColabFold-predicted PPO2/PPX2 protein structures.
                </p>
                <p>
                    This dashboard integrates curated PPO target-site resistance
                    evidence, mutation mapping, and ColabFold-predicted PPO2/PPX2
                    protein structures for Amaranthus species.
                </p>
            </div>
            <div class="feature-panel">
                <div>
                    <h3>Palmer amaranth / waterhemp resistance research prototype</h3>
                    <p>
                        Built for classroom demonstration: evidence-first,
                        transparent about missing pieces, and ready for structural
                        exploration without claiming docking results.
                    </p>
                </div>
                <div class="badge-row">
                    {status_badge("Literature evidence", "blue")}
                    {status_badge("PPO2/PPX2 models", "green")}
                </div>
            </div>
        </div>
        """,
        unsafe_allow_html=True,
    )


def metric_value(frame: pd.DataFrame, column: str) -> int:
    if frame.empty or column not in frame.columns:
        return 0
    return frame[column].replace("", pd.NA).dropna().nunique()


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


def sidebar_navigation() -> None:
    st.sidebar.title("Dashboard Guide")
    st.sidebar.markdown(
        """
        **Quick jumps**

        - [Overview](#overview)
        - [Mutation Evidence](#mutation-evidence)
        - [Protein Models](#protein-models)
        - [Current Status](#current-status)
        - [Future Work / Next Phase](#future-work-next-phase)
        """
    )
    st.sidebar.markdown("---")
    st.sidebar.markdown(
        """
        **Presentation framing**

        Prototype research dashboard  
        Evidence curation implemented  
        3D model viewing implemented  
        Docking planned, not implemented
        """
    )


def render_metric_cards() -> None:
    verified = load_tsv(VERIFIED_EVIDENCE_PATH)
    registry = load_tsv(MODEL_REGISTRY_PATH)

    completed_models = 0
    review_models = 0
    if not registry.empty and "model_status" in registry.columns:
        completed_models = int((registry["model_status"] == "completed").sum())
        review_models = int((registry["model_status"] == "needs_manual_check").sum())

    metric_columns = st.columns(5)
    metric_columns[0].metric("Verified evidence rows", len(verified))
    metric_columns[1].metric("Species count", metric_value(verified, "species"))
    metric_columns[2].metric("Unique mutations", metric_value(verified, "mutation"))
    metric_columns[3].metric("Completed models", completed_models)
    metric_columns[4].metric("Manual-check models", review_models)


def overview_tab() -> None:
    st.markdown('<span id="overview"></span>', unsafe_allow_html=True)
    render_hero()
    render_metric_cards()

    st.markdown('<div class="divider"></div>', unsafe_allow_html=True)
    section_intro(
        "Overview",
        "What this prototype is built to show",
        "The dashboard connects literature-backed PPO resistance evidence with "
        "curated reference sequences, mutation-to-reference residue mapping, and "
        "predicted PPO2/PPX2 protein structures.",
    )

    cards = [
        (
            "Project objective",
            "Organize verified PPO-inhibitor resistance evidence for Amaranthus spp. "
            "and connect that evidence to reference sequences and structural models.",
        ),
        (
            "Why waterhemp and Palmer amaranth matter",
            "These Amaranthus weeds are major resistance-management concerns in row "
            "crop systems, and PPO-inhibitor resistance can complicate herbicide programs.",
        ),
        (
            "Biological question",
            "Which PPO2/PPX2 mutations have verified resistance evidence, and how can "
            "those mutations be placed onto curated reference sequences and predicted "
            "protein structures for follow-up analysis?",
        ),
    ]
    st.markdown(
        '<div class="info-grid">'
        + "".join(
            f'<div class="info-card"><h4>{escape(title)}</h4><p>{escape(body)}</p></div>'
            for title, body in cards
        )
        + "</div>",
        unsafe_allow_html=True,
    )

    st.markdown(
        """
        <div class="section-card">
            <div class="section-kicker">Current State Of The Project</div>
            <div class="muted">
                This is a prototype research dashboard. Evidence curation,
                mutation mapping, and protein structure prediction are implemented.
                Docking, ligand preparation, and herbicide-binding interpretation
                are planned for the next phase.
            </div>
        </div>
        """,
        unsafe_allow_html=True,
    )

    with st.expander("Curated data inputs used by this dashboard"):
        species_summary = load_tsv(SPECIES_SUMMARY_PATH)
        herbicide_summary = load_tsv(HERBICIDE_SUMMARY_PATH)
        references = load_tsv(REFERENCE_SEQUENCES_PATH)
        mappings = load_tsv(MUTATION_MAPPING_PATH)
        st.write(
            {
                "species summary rows": len(species_summary),
                "herbicide summary rows": len(herbicide_summary),
                "reference sequence rows": len(references),
                "mutation mapping rows": len(mappings),
            }
        )


def mutation_evidence_tab() -> None:
    st.markdown('<span id="mutation-evidence"></span>', unsafe_allow_html=True)
    section_intro(
        "Mutation Evidence",
        "Literature-backed resistance evidence",
        "Filters below use the curated mutation evidence table. Verified evidence "
        "is literature-based and remains separate from computational predictions.",
    )

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

    chart_columns = st.columns(2)
    with chart_columns[0]:
        st.markdown("**Evidence rows by species**")
        if "species" in filtered.columns and not filtered.empty:
            species_counts = filtered["species"].value_counts().rename_axis("species")
            st.bar_chart(species_counts)
    with chart_columns[1]:
        st.markdown("**Evidence rows by mutation**")
        if "mutation" in filtered.columns and not filtered.empty:
            mutation_counts = filtered["mutation"].value_counts().rename_axis("mutation")
            st.bar_chart(mutation_counts)

    st.markdown("**Filtered evidence table**")
    st.dataframe(filtered, use_container_width=True, hide_index=True, height=460)

    with st.expander("Summary counts for the filtered table"):
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


def render_model_metadata(selected_model: pd.Series) -> None:
    fields = [
        ("Species", selected_model.get("species", "")),
        ("Mutation", selected_model.get("mutation", "")),
        ("Reference accession", selected_model.get("reference_accession", "")),
        ("Mean pLDDT", selected_model.get("mean_plddt", "")),
        ("PAE present", selected_model.get("pae_file_present", "")),
        ("Model status", selected_model.get("model_status", "")),
    ]
    st.markdown(
        '<div class="model-meta">'
        + "".join(
            f"""
            <div class="meta-cell">
                <div class="meta-label">{escape(label)}</div>
                <div class="meta-value">{escape(str(value))}</div>
            </div>
            """
            for label, value in fields
        )
        + "</div>",
        unsafe_allow_html=True,
    )


def protein_models_tab() -> None:
    st.markdown('<span id="protein-models"></span>', unsafe_allow_html=True)
    section_intro(
        "Protein Models",
        "Interactive PPO2/PPX2 structure exploration",
        "These are predicted PPO2/PPX2 structures from ColabFold. They support "
        "structural exploration, but docking and binding-affinity interpretation "
        "are not implemented in this version.",
    )

    registry = load_tsv(MODEL_REGISTRY_PATH)
    quality = load_tsv(MODEL_QUALITY_PATH)

    if registry.empty:
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
        "dashboard_pdb_path",
    ]
    available_columns = [column for column in display_columns if column in registry.columns]

    filter_columns = st.columns([1, 1, 1.4])
    filtered = registry.copy()
    with filter_columns[0]:
        filtered = multiselect_filter(filtered, "Species", "species", "model_species_filter")
    with filter_columns[1]:
        filtered = multiselect_filter(filtered, "Mutation", "mutation", "model_mutation_filter")
    with filter_columns[2]:
        status_options = sorted(
            value for value in registry.get("model_status", pd.Series(dtype=str)).unique() if value
        )
        selected_status = st.multiselect("Model status", status_options, key="model_status_filter")
        if selected_status:
            filtered = filtered[filtered["model_status"].isin(selected_status)]

    st.dataframe(filtered[available_columns], use_container_width=True, hide_index=True, height=300)

    review_rows = registry[registry["model_status"] == "needs_manual_check"]
    if not review_rows.empty:
        st.warning("APAL_G399A remains needs_manual_check because no dashboard PDB is available.")

    viewable = filtered[filtered.get("dashboard_pdb_path", "").astype(str).str.strip() != ""]
    if viewable.empty:
        st.warning("No dashboard PDB file is available for the selected filters.")
        return

    viewable = viewable.reset_index(drop=True)
    labels = [
        f"{row['model_id']} - {row['species']} - {row['mutation']}"
        for _, row in viewable.iterrows()
    ]
    selected_label = st.selectbox("Choose a model for 3D viewing", labels)
    selected_model = viewable.iloc[labels.index(selected_label)]

    st.markdown("**Selected model metadata**")
    render_model_metadata(selected_model)
    st.caption(f"Best model file: {selected_model.get('best_model_file_name', '')}")

    pdb_path = ROOT / selected_model.get("dashboard_pdb_path", "")
    if not pdb_path.exists() or pdb_path.stat().st_size == 0:
        st.warning("No dashboard PDB file is available for the selected model.")
        return

    pdb_text = pdb_path.read_text(encoding="utf-8", errors="replace")
    viewer = py3Dmol.view(width=940, height=620)
    viewer.addModel(pdb_text, "pdb")
    viewer.setStyle({"cartoon": {"color": "spectrum"}})
    viewer.zoomTo()
    components.html(viewer._make_html(), height=650, scrolling=False)

    if not quality.empty:
        with st.expander("Model quality summary"):
            quality_columns = [
                column
                for column in [
                    "species",
                    "mutation",
                    "sequence_type",
                    "reference_accession",
                    "mean_plddt",
                    "pae_file_present",
                    "model_status",
                ]
                if column in quality.columns
            ]
            st.dataframe(quality[quality_columns], use_container_width=True, hide_index=True)


def current_status_tab() -> None:
    st.markdown('<span id="current-status"></span>', unsafe_allow_html=True)
    section_intro(
        "Current Status",
        "Prototype progress at a glance",
        "The dashboard is functional for evidence review, mutation mapping context, "
        "and 3D model exploration. Docking and ligand workflows remain planned.",
    )

    st.markdown(
        f"""
        <div class="info-grid">
            <div class="status-card completed">
                <h3>Completed</h3>
                <ul>{html_list(COMPLETED_ITEMS)}</ul>
            </div>
            <div class="status-card progress">
                <h3>In Progress</h3>
                <ul>{html_list(IN_PROGRESS_ITEMS)}</ul>
            </div>
            <div class="status-card pending">
                <h3>Not Yet Implemented</h3>
                <ul>{html_list(NOT_IMPLEMENTED_ITEMS)}</ul>
            </div>
        </div>
        """,
        unsafe_allow_html=True,
    )


def future_work_tab() -> None:
    st.markdown('<span id="future-work-next-phase"></span>', unsafe_allow_html=True)
    section_intro(
        "Future Work / Next Phase",
        "From structural viewing to decision-support hypotheses",
        "The next phase should add ligand preparation, receptor setup, validated "
        "docking, and careful interpretation without mixing computational predictions "
        "with experimentally confirmed resistance evidence.",
    )

    for group, items in FUTURE_WORK_GROUPS.items():
        st.subheader(group)
        cols = st.columns(2)
        for index, (title, explanation) in enumerate(items):
            with cols[index % 2]:
                st.markdown(
                    f"""
                    <div class="info-card">
                        <div class="badge-row">{status_badge("planned", "gold")}</div>
                        <h4>{escape(title)}</h4>
                        <p>{escape(explanation)}</p>
                    </div>
                    """,
                    unsafe_allow_html=True,
                )

    st.markdown(
        """
        <div class="section-card">
            <div class="section-kicker">Why This Matters</div>
            <div class="muted">
                The eventual goal is to help compare PPO mutations and explore how
                they may influence herbicide binding and possible cross-resistance
                risk. Any docking-based interpretation should remain clearly labeled
                as computational prediction unless supported by experimental evidence.
            </div>
        </div>
        """,
        unsafe_allow_html=True,
    )


def main() -> None:
    st.set_page_config(
        page_title="PPO-Inhibitor Resistance Dashboard",
        layout="wide",
    )
    inject_css()
    sidebar_navigation()

    tabs = st.tabs(
        [
            "Overview",
            "Mutation Evidence",
            "Protein Models",
            "Current Status",
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
        current_status_tab()
    with tabs[4]:
        future_work_tab()


if __name__ == "__main__":
    main()
