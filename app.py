from __future__ import annotations

from html import escape
from pathlib import Path
import re

import altair as alt
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
DASHBOARD_FIGURES_DIR = ROOT / "structures" / "dashboard_figures"


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
    "ColabFold model generation for all 12 registry models",
    "3D model visualization",
    "APAL_G399A dashboard model and quality figures",
]
IN_PROGRESS_ITEMS = [
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

        header[data-testid="stHeader"],
        [data-testid="stToolbar"],
        [data-testid="stDecoration"],
        [data-testid="stStatusWidget"],
        [data-testid="manage-app-button"],
        .stDeployButton,
        #MainMenu,
        footer {
            display: none !important;
            visibility: hidden !important;
            height: 0 !important;
        }

        [data-testid="collapsedControl"] {
            top: 0.65rem !important;
        }

        div[data-testid="stTabs"] [data-baseweb="tab-list"] {
            border-bottom: 1px solid rgba(23, 35, 29, 0.14);
        }

        [role="tablist"],
        [role="tablist"] * {
            opacity: 1 !important;
            text-shadow: none !important;
        }

        [role="tablist"] [role="tab"],
        [role="tablist"] [role="tab"] *,
        [role="tablist"] button,
        [role="tablist"] button *,
        button[role="tab"],
        button[role="tab"] *,
        [data-baseweb="tab"],
        [data-baseweb="tab"] * {
            color: #17231d !important;
            -webkit-text-fill-color: #17231d !important;
            opacity: 1 !important;
            font-weight: 650 !important;
            text-shadow: none !important;
        }

        [role="tablist"] [role="tab"]:hover,
        [role="tablist"] [role="tab"]:hover *,
        [role="tablist"] button:hover,
        [role="tablist"] button:hover *,
        button[role="tab"]:hover,
        button[role="tab"]:hover *,
        [data-baseweb="tab"]:hover,
        [data-baseweb="tab"]:hover * {
            color: var(--green-dark) !important;
            -webkit-text-fill-color: var(--green-dark) !important;
        }

        [role="tablist"] [role="tab"][aria-selected="true"],
        [role="tablist"] [role="tab"][aria-selected="true"] *,
        [role="tablist"] button[aria-selected="true"],
        [role="tablist"] button[aria-selected="true"] *,
        button[role="tab"][aria-selected="true"],
        button[role="tab"][aria-selected="true"] *,
        [data-baseweb="tab"][aria-selected="true"],
        [data-baseweb="tab"][aria-selected="true"] * {
            color: #ff4b4b !important;
            -webkit-text-fill-color: #ff4b4b !important;
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
            padding-top: 1rem;
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

        div[data-testid="stMetricLabel"],
        div[data-testid="stMetricLabel"] *,
        div[data-testid="stMetricValue"],
        div[data-testid="stMetricValue"] *,
        div[data-testid="stMetricDelta"],
        div[data-testid="stMetricDelta"] * {
            color: #17231d !important;
            -webkit-text-fill-color: #17231d !important;
            opacity: 1 !important;
            text-shadow: none !important;
        }

        div[data-testid="stMetric"] [data-testid],
        div[data-testid="stMetric"] [data-testid] *,
        div[data-testid="stMetric"] div,
        div[data-testid="stMetric"] span,
        div[data-testid="stMetric"] p {
            color: #17231d !important;
            -webkit-text-fill-color: #17231d !important;
            opacity: 1 !important;
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

        .collab-credit {
            display: inline-block;
            margin-top: 1rem;
            padding: 0.65rem 0.85rem;
            border-left: 4px solid var(--green);
            border-radius: 8px;
            background: rgba(255,255,255,0.76);
            color: #304239;
            font-size: 0.92rem;
            line-height: 1.45;
        }

        .collab-credit strong {
            color: #173523;
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

        .info-card.planned {
            border-left-color: var(--gold);
            background: #fffaf0;
        }

        .info-card.done {
            border-left-color: var(--green);
            background: #f6fbf6;
        }

        .method-step {
            background: #ffffff;
            border: 1px solid var(--line);
            border-left: 5px solid var(--blue);
            border-radius: 12px;
            padding: 1rem;
            box-shadow: 0 8px 22px rgba(31,55,41,0.07);
        }

        .method-step.done {
            border-left-color: var(--green);
            background: #f6fbf6;
        }

        .method-step.planned {
            border-left-color: var(--gold);
            background: #fffaf0;
        }

        .method-step h4 {
            margin: 0 0 0.4rem 0;
            color: #1d3528;
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

        .plddt-legend {
            background: #ffffff;
            border: 1px solid var(--line);
            border-radius: 12px;
            padding: 1rem 1.1rem;
            box-shadow: 0 10px 24px rgba(31,55,41,0.08);
            margin: 0.55rem 0 1rem 0;
        }

        .plddt-scale {
            display: flex;
            flex-wrap: wrap;
            align-items: center;
            gap: 0.7rem 1.1rem;
            margin: 0.35rem 0 0.75rem 0;
        }

        .plddt-title {
            font-weight: 900;
            color: #17231d;
            margin-right: 0.15rem;
        }

        .plddt-item {
            display: inline-flex;
            align-items: center;
            gap: 0.45rem;
            color: #2f3b33;
            font-weight: 700;
            white-space: nowrap;
        }

        .plddt-swatch {
            width: 36px;
            height: 13px;
            border-radius: 2px;
            border: 1px solid rgba(23,35,29,0.16);
            display: inline-block;
        }

        .compare-divider {
            height: 610px;
            margin: 3rem auto 0 auto;
            width: 1px;
            background: linear-gradient(180deg, transparent, rgba(31,111,74,0.35), transparent);
            position: relative;
        }

        .compare-divider::after {
            content: "Compare";
            position: absolute;
            top: 45%;
            left: 50%;
            transform: translate(-50%, -50%) rotate(-90deg);
            padding: 0.25rem 0.55rem;
            border-radius: 999px;
            background: #eef7ec;
            border: 1px solid rgba(31,111,74,0.22);
            color: #1f6f4a;
            font-size: 0.74rem;
            font-weight: 850;
            letter-spacing: 0.06em;
            text-transform: uppercase;
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
                <h1>PPO-Inhibitor Resistance Dashboard for <em>Amaranthus</em> spp.</h1>
                <p>
                    This project brings together PPO-inhibitor resistance evidence
                    for <em>Amaranthus</em> weeds and connects each verified mutation to
                    reference sequences, residue mapping, and predicted PPO2/PPX2
                    protein structures.
                </p>
                <p>
                    This dashboard integrates curated PPO target-site resistance
                    evidence, mutation mapping, and ColabFold-predicted PPO2/PPX2
                    protein structures for <em>Amaranthus</em> species.
                </p>
                <div class="collab-credit">
                    <strong>Collaborative idea:</strong> Gourav Chahal, PhD student,
                    Auburn University; and Rishabh Singh, PhD student,
                    University of Illinois Urbana-Champaign.
                </div>
            </div>
            <div class="feature-panel">
                <div>
                    <h3>Palmer amaranth / waterhemp resistance research prototype</h3>
                    <p>
                        The focus is target-site resistance in Palmer amaranth,
                        waterhemp, and related <em>Amaranthus</em> species: which PPO
                        mutations are supported by literature evidence, where they
                        map on the protein, and which predicted structures are
                        ready for careful follow-up.
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


def abbreviated_species_name(species: str) -> str:
    parts = species.split()
    if len(parts) >= 2:
        return f"{parts[0][0]}. {parts[1]}"
    return species


def evidence_count_chart(
    frame: pd.DataFrame,
    category_column: str,
    label_column: str,
    italic_labels: bool = False,
    label_angle: int = 0,
) -> alt.Chart:
    counts = (
        frame.groupby(category_column, dropna=False)
        .size()
        .reset_index(name="evidence_rows")
        .sort_values("evidence_rows", ascending=False)
    )
    if label_column != category_column:
        counts[label_column] = counts[category_column].map(abbreviated_species_name)

    axis = alt.Axis(
        labelAngle=label_angle,
        labelAlign="right" if label_angle else "center",
        labelBaseline="middle" if label_angle else "top",
        labelFontSize=12,
        labelFontStyle="italic" if italic_labels else "normal",
        labelLimit=140,
        labelOverlap=False,
        title=None,
    )

    return (
        alt.Chart(counts)
        .mark_bar(color="#1f6f4a", cornerRadiusTopLeft=4, cornerRadiusTopRight=4)
        .encode(
            x=alt.X(f"{label_column}:N", sort=None, axis=axis),
            y=alt.Y(
                "evidence_rows:Q",
                title="Evidence rows",
                axis=alt.Axis(grid=True, labelFontSize=12),
            ),
            tooltip=[
                alt.Tooltip(f"{category_column}:N", title=category_column.replace("_", " ").title()),
                alt.Tooltip("evidence_rows:Q", title="Evidence rows"),
            ],
        )
        .properties(height=330)
    )


def sidebar_navigation() -> None:
    st.sidebar.title("Dashboard Guide")
    st.sidebar.markdown(
        """
        **Quick jumps**

        - [Overview](#overview)
        - [Mutation Evidence](#mutation-evidence)
        - [Protein Models](#protein-models)
        - [Methodology](#methodology)
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
    metric_columns[4].metric("Models needing review", review_models)


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
            "Organize verified PPO-inhibitor resistance evidence for <em>Amaranthus</em> spp. "
            "and connect that evidence to reference sequences and structural models.",
        ),
        (
            "Why waterhemp and Palmer amaranth matter",
            "These <em>Amaranthus</em> weeds are major resistance-management concerns in row "
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
            f'<div class="info-card"><h4>{escape(title)}</h4><p>{body}</p></div>'
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
            st.altair_chart(
                evidence_count_chart(
                    filtered,
                    category_column="species",
                    label_column="species_label",
                    italic_labels=True,
                    label_angle=0,
                ),
                use_container_width=True,
            )
    with chart_columns[1]:
        st.markdown("**Evidence rows by mutation**")
        if "mutation" in filtered.columns and not filtered.empty:
            st.altair_chart(
                evidence_count_chart(
                    filtered,
                    category_column="mutation",
                    label_column="mutation",
                    italic_labels=False,
                    label_angle=-45,
                ),
                use_container_width=True,
            )

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
    metadata_cells = "".join(
        f'<div class="meta-cell">'
        f'<div class="meta-label">{escape(label)}</div>'
        f'<div class="meta-value">{escape(str(value))}</div>'
        f"</div>"
        for label, value in fields
    )
    st.markdown(
        f'<div class="model-meta">{metadata_cells}</div>',
        unsafe_allow_html=True,
    )


def viewer_style_controls() -> dict[str, str | bool]:
    st.markdown("**3D viewer display controls**")
    control_columns = st.columns(4)
    with control_columns[0]:
        mainchain_style = st.selectbox(
            "Main chain",
            ["Cartoon", "Backbone sticks", "Trace", "Hide"],
            index=0,
        )
    with control_columns[1]:
        sidechain_style = st.selectbox(
            "Side chains",
            ["Hide", "Sticks", "Lines"],
            index=0,
        )
    with control_columns[2]:
        color_scheme = st.selectbox(
            "Color",
            ["Spectrum", "Confidence-style", "Green", "Slate"],
            index=1,
        )
    with control_columns[3]:
        background = st.selectbox(
            "Background",
            ["White", "Light", "Dark"],
            index=0,
        )
    option_columns = st.columns(2)
    with option_columns[0]:
        show_surface = st.checkbox("Show translucent surface", value=False)
    with option_columns[1]:
        auto_rotate = st.checkbox("Auto-rotate model(s)", value=False)
    return {
        "mainchain_style": mainchain_style,
        "sidechain_style": sidechain_style,
        "color_scheme": color_scheme,
        "background": background,
        "show_surface": show_surface,
        "auto_rotate": auto_rotate,
    }


def color_for_scheme(color_scheme: str) -> str | dict[str, str]:
    if color_scheme == "Spectrum":
        return "spectrum"
    if color_scheme == "Confidence-style":
        return {"prop": "b", "gradient": "roygb", "min": 50, "max": 100}
    if color_scheme == "Green":
        return "#2f7d4f"
    return "#5d6a73"


def style_color_kwargs(color: str | dict[str, str]) -> dict[str, str | dict[str, str]]:
    if isinstance(color, dict):
        return {"colorscheme": color}
    return {"color": color}


def add_viewer_styles(viewer: py3Dmol.view, controls: dict[str, str | bool]) -> None:
    color = color_for_scheme(str(controls["color_scheme"]))
    color_kwargs = style_color_kwargs(color)
    mainchain_style = str(controls["mainchain_style"])
    sidechain_style = str(controls["sidechain_style"])
    backbone_atoms = ["N", "CA", "C", "O"]
    sidechain_selection = {"not": {"atom": backbone_atoms}}

    if mainchain_style == "Cartoon":
        viewer.setStyle({"cartoon": color_kwargs})
    elif mainchain_style == "Backbone sticks":
        viewer.setStyle(
            {"atom": backbone_atoms},
            {"stick": {**color_kwargs, "radius": 0.18}},
        )
    elif mainchain_style == "Trace":
        viewer.setStyle({"atom": "CA"}, {"sphere": {**color_kwargs, "radius": 0.45}})
    else:
        viewer.setStyle({})

    if sidechain_style == "Sticks":
        viewer.addStyle(
            sidechain_selection,
            {"stick": {**color_kwargs, "radius": 0.14}},
        )
    elif sidechain_style == "Lines":
        viewer.addStyle(
            sidechain_selection,
            {"line": color_kwargs},
        )

    if bool(controls["show_surface"]):
        viewer.addSurface(
            py3Dmol.VDW,
            {"opacity": 0.18, "color": "#78a883"},
        )

    background = str(controls["background"])
    if background == "Dark":
        viewer.setBackgroundColor("#111816")
    elif background == "Light":
        viewer.setBackgroundColor("#eef4ec")
    else:
        viewer.setBackgroundColor("#ffffff")


def render_3d_viewer_legend(controls: dict[str, str | bool]) -> None:
    color_mode = str(controls["color_scheme"])
    if color_mode == "Spectrum":
        color_explanation = (
            "Spectrum colors the protein as a rainbow along the model, which helps "
            "visually follow the fold from one region to another."
        )
    elif color_mode == "Confidence-style":
        color_explanation = (
            "Confidence-style uses the confidence values stored in the PDB B-factor "
            "field from the ColabFold output; warmer colors are lower confidence and "
            "cooler blue/green colors are higher confidence in this viewer."
        )
    elif color_mode == "Green":
        color_explanation = "Green is a single-color structural view for clean presentation."
    else:
        color_explanation = "Slate is a single-color structural view for clean presentation."

    legend_cells = "".join(
        [
            (
                f'<div class="meta-cell"><div class="meta-label">Main chain</div>'
                f'<div class="meta-value">{escape(str(controls["mainchain_style"]))}</div></div>'
            ),
            (
                f'<div class="meta-cell"><div class="meta-label">Side chains</div>'
                f'<div class="meta-value">{escape(str(controls["sidechain_style"]))}</div></div>'
            ),
            (
                f'<div class="meta-cell"><div class="meta-label">Color mode</div>'
                f'<div class="meta-value">{escape(str(controls["color_scheme"]))}</div></div>'
            ),
            (
                f'<div class="meta-cell"><div class="meta-label">Surface</div>'
                f'<div class="meta-value">{"Shown" if controls["show_surface"] else "Hidden"}</div></div>'
            ),
            (
                f'<div class="meta-cell"><div class="meta-label">Auto-rotate</div>'
                f'<div class="meta-value">{"On" if controls.get("auto_rotate", False) else "Off"}</div></div>'
            ),
        ]
    )
    st.markdown(
        f'<div class="section-card">'
        f'<div class="section-kicker">3D Viewer Legend</div>'
        f'<div class="info-grid">{legend_cells}</div>'
        f'<p class="muted" style="margin-top:0.85rem;">'
        f"Cartoon, trace, and backbone styles are display choices for exploring "
        f"the predicted PPO2/PPX2 model. Side-chain sticks or lines can be added "
        f"for closer inspection. These views do not represent docking or binding-affinity results."
        f"</p>"
        f'<p class="muted" style="margin-top:0.45rem;">'
        f"<strong>Corner axis marker:</strong> The X/Y/Z triad stays fixed in the viewer "
        f"corner and follows the model view while ignoring zoom. It is an orientation guide only, "
        f"not a biological coordinate system or residue measurement."
        f"</p>"
        f'<p class="muted" style="margin-top:0.45rem;">'
        f"<strong>Color meaning:</strong> {escape(color_explanation)}"
        f"</p>"
        f"</div>",
        unsafe_allow_html=True,
    )


def render_plddt_confidence_legend(controls: dict[str, str | bool]) -> None:
    active_note = (
        "This pLDDT scale matches the current Confidence-style model coloring."
        if controls["color_scheme"] == "Confidence-style"
        else "This pLDDT scale applies when the model color control is set to Confidence-style."
    )
    legend_items = [
        ("#ff2400", "Very low (&lt;50)"),
        ("#ffff33", "Low (60)"),
        ("#37f43b", "OK (70)"),
        ("#36e6e6", "Confident (80)"),
        ("#0615e8", "Very high (&gt;90)"),
    ]
    legend_html = "".join(
        f'<span class="plddt-item">'
        f'<span class="plddt-swatch" style="background:{color};"></span>'
        f"<span>{label}</span>"
        f"</span>"
        for color, label in legend_items
    )
    st.markdown(
        f'<div class="plddt-legend">'
        f'<div class="plddt-scale"><span class="plddt-title">pLDDT:</span>{legend_html}</div>'
        f'<div class="muted"><strong>Caption:</strong> {active_note} pLDDT is the '
        f"ColabFold/AlphaFold per-residue confidence score used for structure model "
        f"confidence. This legend describes model confidence only; it is not docking, "
        f"herbicide-binding, or field-resistance evidence.</div>"
        f"</div>",
        unsafe_allow_html=True,
    )


def add_corner_axis_overlay(viewer_html: str, auto_rotate: bool = False) -> str:
    container_match = re.search(r'id="(3dmolviewer_[^"]+)"', viewer_html)
    viewer_match = re.search(r"var (viewer_[A-Za-z0-9_]+) = null;", viewer_html)
    if not container_match or not viewer_match:
        return viewer_html

    container_id = container_match.group(1)
    viewer_var = viewer_match.group(1)
    auto_rotate_js = "true" if auto_rotate else "false"
    overlay_script = f"""
<script>
(function() {{
  const container = document.getElementById("{container_id}");
  if (!container) return;
  const canvas = document.createElement("canvas");
  canvas.width = 150;
  canvas.height = 150;
  canvas.style.position = "absolute";
  canvas.style.left = "18px";
  canvas.style.bottom = "18px";
  canvas.style.width = "112px";
  canvas.style.height = "112px";
  canvas.style.pointerEvents = "none";
  canvas.style.background = "transparent";
  canvas.style.border = "0";
  canvas.style.boxShadow = "none";
  canvas.style.zIndex = "20";
  container.appendChild(canvas);

  const ctx = canvas.getContext("2d");
  const center = {{ x: 62, y: 95 }};
  const axisLength = 42;
  const autoRotate = {auto_rotate_js};
  let lastSpin = 0;
  const axes = [
    {{ label: "X", color: "#d43f3a", vector: [1, 0, 0] }},
    {{ label: "Y", color: "#2f7d4f", vector: [0, 1, 0] }},
    {{ label: "Z", color: "#244f9e", vector: [0, 0, 1] }}
  ];

  function viewerInstance() {{
    return window["{viewer_var}"] || null;
  }}

  function normalizeQuaternion(q) {{
    const length = Math.sqrt(q.x * q.x + q.y * q.y + q.z * q.z + q.w * q.w) || 1;
    return {{ x: q.x / length, y: q.y / length, z: q.z / length, w: q.w / length }};
  }}

  function rotateVector(v, q) {{
    const x = v[0], y = v[1], z = v[2];
    const qx = q.x, qy = q.y, qz = q.z, qw = q.w;
    const ix =  qw * x + qy * z - qz * y;
    const iy =  qw * y + qz * x - qx * z;
    const iz =  qw * z + qx * y - qy * x;
    const iw = -qx * x - qy * y - qz * z;
    return [
      ix * qw + iw * -qx + iy * -qz - iz * -qy,
      iy * qw + iw * -qy + iz * -qx - ix * -qz,
      iz * qw + iw * -qz + ix * -qy - iy * -qx
    ];
  }}

  function currentQuaternion() {{
    const viewer = viewerInstance();
    if (viewer && typeof viewer.getView === "function") {{
      const view = viewer.getView();
      if (view && view.length >= 8) {{
        return normalizeQuaternion({{ x: view[4] || 0, y: view[5] || 0, z: view[6] || 0, w: view[7] || 1 }});
      }}
      if (view && view.length >= 7) {{
        return normalizeQuaternion({{ x: view[3] || 0, y: view[4] || 0, z: view[5] || 0, w: view[6] || 1 }});
      }}
    }}
    return {{ x: 0, y: 0, z: 0, w: 1 }};
  }}

  function projectedAxis(axis, q) {{
    const rotated = rotateVector(axis.vector, q);
    return {{
      ...axis,
      rotated,
      end: {{
        x: center.x + rotated[0] * axisLength,
        y: center.y - rotated[1] * axisLength - rotated[2] * 10
      }},
      start: {{
        x: center.x - rotated[0] * axisLength * 0.45,
        y: center.y + rotated[1] * axisLength * 0.45 + rotated[2] * 4
      }}
    }};
  }}

  function drawArrow(axis) {{
    const rotated = axis.rotated;
    const depth = Math.max(0.28, Math.min(1, (rotated[2] + 1.45) / 2.45));
    const start = axis.start;
    const end = axis.end;
    ctx.strokeStyle = axis.color;
    ctx.fillStyle = axis.color;
    ctx.globalAlpha = depth;
    ctx.lineWidth = 5.4;
    ctx.lineCap = "round";
    ctx.shadowColor = "rgba(255,255,255,0.92)";
    ctx.shadowBlur = 4;
    ctx.shadowOffsetY = 0;
    ctx.beginPath();
    ctx.moveTo(start.x, start.y);
    ctx.lineTo(end.x, end.y);
    ctx.stroke();
    ctx.shadowColor = "transparent";

    const angle = Math.atan2(end.y - center.y, end.x - center.x);
    ctx.beginPath();
    ctx.moveTo(end.x, end.y);
    ctx.lineTo(end.x - 11 * Math.cos(angle - Math.PI / 6), end.y - 11 * Math.sin(angle - Math.PI / 6));
    ctx.lineTo(end.x - 11 * Math.cos(angle + Math.PI / 6), end.y - 11 * Math.sin(angle + Math.PI / 6));
    ctx.closePath();
    ctx.fill();

    ctx.globalAlpha = 1;
    ctx.font = "bold 17px Arial, sans-serif";
    ctx.textAlign = "center";
    ctx.textBaseline = "middle";
    ctx.shadowColor = "rgba(255,255,255,0.95)";
    ctx.shadowBlur = 5;
    ctx.fillStyle = axis.color;
    ctx.fillText(axis.label, end.x + Math.sign(rotated[0] || 1) * 12, end.y - Math.sign(rotated[1] || 1) * 11);
    ctx.shadowColor = "transparent";
  }}

  function maybeRotate(viewer, timestamp) {{
    if (!autoRotate || !viewer) return;
    if (timestamp - lastSpin < 38) return;
    lastSpin = timestamp;
    try {{
      if (typeof viewer.rotate === "function") {{
        viewer.rotate(1.0, "y");
        viewer.render();
      }}
    }} catch (error) {{}}
  }}

  function draw(timestamp) {{
    const viewer = viewerInstance();
    maybeRotate(viewer, timestamp || 0);

    ctx.clearRect(0, 0, canvas.width, canvas.height);
    const q = currentQuaternion();
    ctx.save();
    ctx.strokeStyle = "rgba(255,255,255,0.82)";
    ctx.lineWidth = 8;
    ctx.lineCap = "round";
    axes
      .map(axis => projectedAxis(axis, q))
      .sort((a, b) => a.rotated[2] - b.rotated[2])
      .forEach(axis => drawArrow(axis));
    ctx.restore();
    ctx.fillStyle = "#17231d";
    ctx.shadowColor = "rgba(255,255,255,0.95)";
    ctx.shadowBlur = 4;
    ctx.beginPath();
    ctx.arc(center.x, center.y, 5, 0, Math.PI * 2);
    ctx.fill();
    ctx.shadowColor = "transparent";
    window.requestAnimationFrame(draw);
  }}
  draw();
}})();
</script>
"""
    return viewer_html + overlay_script


def render_structure_viewer(
    selected_model: pd.Series,
    controls: dict[str, str | bool],
    height: int = 620,
    viewer_width: int | str = "100%",
) -> bool:
    pdb_path = ROOT / selected_model.get("dashboard_pdb_path", "")
    if not pdb_path.exists() or pdb_path.stat().st_size == 0:
        st.warning("No dashboard PDB file is available for the selected model.")
        return False

    pdb_text = pdb_path.read_text(encoding="utf-8", errors="replace")
    viewer = py3Dmol.view(width=viewer_width, height=height)
    viewer.addModel(pdb_text, "pdb")
    add_viewer_styles(viewer, controls)
    viewer.zoomTo()
    viewer_html = add_corner_axis_overlay(
        viewer._make_html(),
        auto_rotate=bool(controls.get("auto_rotate", False)),
    )
    components.html(viewer_html, height=height + 30, scrolling=False)
    return True


def model_figure_paths(model_id: str) -> dict[str, Path]:
    return {
        "coverage": DASHBOARD_FIGURES_DIR / f"{model_id}_coverage.png",
        "plddt": DASHBOARD_FIGURES_DIR / f"{model_id}_plddt.png",
        "pae": DASHBOARD_FIGURES_DIR / f"{model_id}_pae.png",
    }


def render_model_figures(selected_model: pd.Series) -> None:
    model_id = selected_model.get("model_id", "")
    paths = model_figure_paths(model_id)
    available = {key: path for key, path in paths.items() if path.exists() and path.stat().st_size}
    if not available:
        st.info("No ColabFold summary figures are available for this selected model.")
        return

    st.markdown("**ColabFold model summary figures**")
    captions = {
        "plddt": (
            "pLDDT confidence plot. Higher pLDDT values indicate higher model "
            "confidence for those residue positions."
        ),
        "pae": (
            "Predicted aligned error (PAE) plot. Lower values indicate higher "
            "confidence in relative domain/residue placement."
        ),
        "coverage": (
            "Sequence coverage plot from the ColabFold run. This summarizes MSA "
            "coverage supporting the prediction."
        ),
    }
    ordered_keys = ["plddt", "pae", "coverage"]
    columns = st.columns(3)
    for index, key in enumerate(ordered_keys):
        if key not in available:
            continue
        with columns[index]:
            st.image(str(available[key]), use_container_width=True)
            st.caption(captions[key])


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
        review_ids = ", ".join(review_rows["model_id"].astype(str).tolist())
        st.warning(
            f"{len(review_rows)} model row(s) still need manual review before dashboard viewing: "
            f"{review_ids}."
        )

    viewable = filtered[filtered.get("dashboard_pdb_path", "").astype(str).str.strip() != ""]
    if viewable.empty:
        st.warning("No dashboard PDB file is available for the selected filters.")
        return

    viewable = viewable.reset_index(drop=True)
    labels = [
        f"{row['model_id']} - {row['species']} - {row['mutation']}"
        for _, row in viewable.iterrows()
    ]
    viewer_mode = st.selectbox(
        "Viewer mode",
        ["Single model", "Compare two models"],
        index=0,
    )
    viewer_controls = viewer_style_controls()
    render_3d_viewer_legend(viewer_controls)

    if viewer_mode == "Compare two models":
        st.markdown("**Compare two predicted PPO2/PPX2 models**")
        st.caption(
            "Side-by-side comparison uses the same viewer settings for both structures. "
            "The models are not structurally aligned or superposed, and this is not a docking comparison."
        )
        compare_columns = st.columns([1, 0.06, 1])
        with compare_columns[0]:
            left_label = st.selectbox(
                "Left model",
                labels,
                index=0,
                key="left_compare_model",
            )
            left_model = viewable.iloc[labels.index(left_label)]
            render_model_metadata(left_model)
            render_structure_viewer(left_model, viewer_controls, height=520)
        with compare_columns[1]:
            st.markdown('<div class="compare-divider"></div>', unsafe_allow_html=True)
        with compare_columns[2]:
            right_default = 1 if len(labels) > 1 else 0
            right_label = st.selectbox(
                "Right model",
                labels,
                index=right_default,
                key="right_compare_model",
            )
            right_model = viewable.iloc[labels.index(right_label)]
            render_model_metadata(right_model)
            render_structure_viewer(right_model, viewer_controls, height=520)
        render_plddt_confidence_legend(viewer_controls)
        return

    selected_label = st.selectbox("Choose a model for 3D viewing", labels)
    selected_model = viewable.iloc[labels.index(selected_label)]

    st.markdown("**Selected model metadata**")
    render_model_metadata(selected_model)
    st.caption(f"Best model file: {selected_model.get('best_model_file_name', '')}")
    if render_structure_viewer(selected_model, viewer_controls, height=620):
        render_plddt_confidence_legend(viewer_controls)
        render_model_figures(selected_model)

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
        f'<div class="info-grid">'
        f'<div class="status-card completed"><h3>Completed</h3><ul>{html_list(COMPLETED_ITEMS)}</ul></div>'
        f'<div class="status-card progress"><h3>In Progress</h3><ul>{html_list(IN_PROGRESS_ITEMS)}</ul></div>'
        f'<div class="status-card pending"><h3>Not Yet Implemented</h3><ul>{html_list(NOT_IMPLEMENTED_ITEMS)}</ul></div>'
        f"</div>",
        unsafe_allow_html=True,
    )


def methodology_tab() -> None:
    st.markdown('<span id="methodology"></span>', unsafe_allow_html=True)
    section_intro(
        "Methodology",
        "How the project was built and where it goes next",
        "This workflow explains how literature evidence, NCBI reference accessions, "
        "mutation mapping, and ColabFold protein models were connected. Planned "
        "ligand and docking work is shown separately so computational next steps "
        "are not mixed with verified resistance evidence.",
    )

    st.markdown(
        """
        <div class="section-card">
            <div class="section-kicker">Research Workflow</div>
            <div class="muted">
                The dashboard starts with experimentally reported PPO-inhibitor
                resistance evidence, links each usable mutation record to curated
                PPO2/PPX2 reference sequences, then displays ColabFold-predicted
                structures for structural exploration. Docking has not been run in
                this version.
            </div>
        </div>
        """,
        unsafe_allow_html=True,
    )

    completed_steps = [
        (
            "1. Literature and paper screening",
            "Papers were used when they directly supported PPO-inhibitor resistance "
            "evidence in <em>Amaranthus</em> spp. The useful records had to give "
            "enough information to connect a species, PPO gene or isoform, mutation, "
            "resistance evidence, and citation metadata.",
        ),
        (
            "2. Mutation evidence decisions",
            "A mutation was marked as verified only for a specific mutation-paper "
            "combination when the evidence was clear. Background mentions, uncertain "
            "records, or computational-only statements stayed separate and could remain "
            "<code>needs_manual_check</code>.",
        ),
        (
            "3. NCBI accession and reference sequence curation",
            "NCBI protein and nucleotide accession records were used to choose the "
            "working PPO2/PPX2 references. This was needed because mutation positions "
            "must be tied to a specific accession instead of mixing species, PPO "
            "isoforms, or residue numbering systems.",
        ),
        (
            "4. Mutation-to-reference residue mapping",
            "The curated mutation records were connected to the selected reference "
            "accessions using the project mapping tables. This keeps Palmer amaranth, "
            "waterhemp, and other <em>Amaranthus</em> records traceable to the exact "
            "reference sequence used for downstream modeling.",
        ),
        (
            "5. ColabFold / AlphaFold-style protein modeling",
            "Wild-type and mutant PPO2/PPX2 FASTA files were submitted to ColabFold. "
            "ColabFold uses an AlphaFold-style structure-prediction workflow with "
            "sequence-search support, producing ranked model files and confidence "
            "outputs for each submitted sequence.",
        ),
        (
            "6. Model registry and dashboard visualization",
            "The full ColabFold ZIP outputs stay outside the repository. The dashboard "
            "uses a registry with model status, mean pLDDT, PAE availability, external "
            "ZIP path, and the selected rank_001 PDB path for local 3D viewing. Small "
            "summary figures for pLDDT, PAE, and sequence coverage are shown as model "
            "quality context.",
        ),
    ]

    planned_steps = [
        (
            "7. Ligand records and 3D ligand preparation",
            "The next phase should curate herbicide ligand records with PubChem IDs, "
            "then prepare ligand structures in a documented way. Candidate tools to "
            "evaluate include PubChem plus RDKit or Open Babel for structure cleanup "
            "and conversion.",
        ),
        (
            "8. Receptor preparation and binding-site setup",
            "Predicted PPO2/PPX2 receptor structures will need documented preparation "
            "before docking. That includes receptor inspection, cleanup choices, "
            "binding-pocket or grid definition, and a clear FAD/cofactor decision. "
            "PyMOL or ChimeraX are candidate tools to evaluate for inspection.",
        ),
        (
            "9. Docking validation and WT-versus-mutant comparison",
            "Only after receptor and ligand preparation are documented should a docking "
            "engine such as AutoDock Vina or a similar tool be evaluated. Any docking "
            "scores would be treated as computational predictions, not proof of field "
            "resistance or cross-resistance.",
        ),
        (
            "10. Dashboard interpretation layer",
            "Validated docking outputs can later be added to the dashboard as a separate "
            "computational layer beside the literature evidence. The goal is to help "
            "compare mutations and generate hypotheses without changing the verified "
            "evidence calls.",
        ),
    ]

    methodology_groups = [
        ("Completed So Far", "done", completed_steps),
        ("Planned Next Phase", "planned", planned_steps),
    ]

    for group_title, status, steps in methodology_groups:
        badge_label = "completed" if status == "done" else "planned"
        badge_color = "green" if status == "done" else "gold"
        method_cards = "".join(
            f'<div class="method-step {status}">'
            f'<div class="badge-row">{status_badge(badge_label, badge_color)}</div>'
            f"<h4>{escape(title)}</h4>"
            f"<p>{body}</p>"
            f"</div>"
            for title, body in steps
        )
        st.markdown(
            f'<div class="badge-row">{status_badge(group_title, badge_color)}</div>',
            unsafe_allow_html=True,
        )
        st.markdown(
            f'<div class="info-grid">{method_cards}</div>',
            unsafe_allow_html=True,
        )

    st.markdown(
        """
        <div class="section-card">
            <div class="section-kicker">Important Interpretation Rule</div>
            <div class="muted">
                Literature evidence, sequence mapping, and predicted structures are
                implemented. Ligand preparation, docking, binding-affinity summaries,
                and cross-resistance interpretation are future work and should remain
                labeled as computational analysis once added.
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
                future_card = (
                    f'<div class="info-card planned">'
                    f'<div class="badge-row">{status_badge("planned", "gold")}</div>'
                    f"<h4>{escape(title)}</h4>"
                    f"<p>{escape(explanation)}</p>"
                    f"</div>"
                )
                st.markdown(
                    future_card,
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
            "Methodology",
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
        methodology_tab()
    with tabs[4]:
        current_status_tab()
    with tabs[5]:
        future_work_tab()


if __name__ == "__main__":
    main()
