import streamlit as st
import pandas as pd
import numpy as np
import os
import glob


# ============================================================
# PAGE CONFIGURATION
# ============================================================

st.set_page_config(
    page_title="OncoGeneX",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded"
)


# ============================================================
# CONSTANTS
# ============================================================

PROJECT_NAME = "OncoGeneX"
DATASET_ID = "GSE40791"
CANCER_TYPE = "Lung Adenocarcinoma"

P_VALUE_THRESHOLD = 0.05
LOG2FC_THRESHOLD = 1.0


# ============================================================
# CUSTOM CSS
# ============================================================

st.markdown("""
<style>

    /* ================================
       THEME 2 — PURPLE / MAGENTA
    ================================= */

    .stApp {
        background:
            radial-gradient(circle at 82% 12%, rgba(119, 63, 190, 0.18), transparent 28%),
            radial-gradient(circle at 18% 82%, rgba(208, 62, 180, 0.08), transparent 30%),
            #100a1d;
        color: #f2edf8;
    }

    .main .block-container {
        max-width: 1450px;
        padding-top: 2.5rem;
        padding-bottom: 3rem;
    }

    h1, h2, h3 {
        color: #f6f0ff !important;
    }

    p, label, span, div {
        color: #e4dceb;
    }

    /* ================================
       SIDEBAR
    ================================= */

    [data-testid="stSidebar"] {
        background: linear-gradient(180deg, #170b2b 0%, #23113b 100%);
        border-right: 1px solid #513273;
    }

    [data-testid="stSidebar"] * {
        color: #eee6f7;
    }

    /* ================================
       TITLE
    ================================= */

    .main-title {
        font-size: 3.2rem;
        font-weight: 800;
        margin-bottom: 0.3rem;
        background: linear-gradient(90deg, #ffffff, #e6d8ff, #c07cff);
        -webkit-background-clip: text;
        -webkit-text-fill-color: transparent;
    }

    .subtitle {
        font-size: 1.2rem;
        color: #cdb8e3 !important;
        margin-bottom: 2rem;
    }

    /* ================================
       DASHBOARD CARDS
    ================================= */

    .dashboard-card {
        background: linear-gradient(145deg, rgba(42, 24, 66, 0.96), rgba(24, 13, 43, 0.96));
        border: 1px solid #68458d;
        border-radius: 18px;
        padding: 1.5rem;
        box-shadow: 0 8px 30px rgba(0, 0, 0, 0.30), inset 0 1px 0 rgba(255,255,255,0.05);
        min-height: 130px;
    }

    .dashboard-card:hover {
        border: 1px solid #c262e5;
        transform: translateY(-2px);
        transition: 0.2s;
    }

    /* ================================
       GENE RESULT CARDS
    ================================= */

    .gene-card {
        background: linear-gradient(145deg, #2a1842, #180d2b);
        border: 1px solid #76539a;
        border-radius: 20px;
        padding: 1.5rem;
        box-shadow: 0 10px 30px rgba(0,0,0,0.30);
        min-height: 145px;
    }

    .gene-card-label {
        font-size: 0.95rem;
        color: #cdb8e3 !important;
        margin-bottom: 0.8rem;
    }

    .gene-card-value {
        font-size: 2.1rem;
        font-weight: 700;
        color: #ffffff !important;
    }

    /* ================================
       STATUS BOXES
    ================================= */

    .status-significant {
        background: linear-gradient(90deg, rgba(46, 137, 102, 0.30), rgba(21, 79, 62, 0.18));
        border: 1px solid #43c59a;
        border-radius: 16px;
        padding: 1.2rem;
        font-size: 1.1rem;
        font-weight: 700;
        color: #8cf0cc !important;
    }

    .status-up {
        background: linear-gradient(90deg, rgba(119, 78, 210, 0.32), rgba(63, 34, 126, 0.18));
        border: 1px solid #a979ff;
        border-radius: 16px;
        padding: 1.2rem;
        color: #d4bcff !important;
        font-weight: 700;
    }

    .status-down {
        background: linear-gradient(90deg, rgba(194, 45, 130, 0.34), rgba(106, 18, 76, 0.20));
        border: 1px solid #ef68ba;
        border-radius: 16px;
        padding: 1.2rem;
        color: #ffc0e3 !important;
        font-weight: 700;
    }

    .status-not-significant {
        background: linear-gradient(90deg, rgba(154, 105, 31, 0.30), rgba(93, 59, 10, 0.18));
        border: 1px solid #c99842;
        border-radius: 16px;
        padding: 1.2rem;
        color: #ffe0a0 !important;
        font-weight: 700;
    }

    .info-box {
        background: rgba(71, 38, 103, 0.48);
        border-left: 4px solid #bd75ef;
        padding: 1rem 1.2rem;
        border-radius: 10px;
        color: #eadcff !important;
    }

    /* ================================
       INPUTS
    ================================= */

    .stTextInput input,
    .stSelectbox div[data-baseweb="select"] > div {
        background-color: #24143b !important;
        color: #ffffff !important;
        border-color: #76519a !important;
    }

    .stButton > button {
        background: linear-gradient(90deg, #7a3db2, #b353d3);
        color: white !important;
        border: none;
        border-radius: 10px;
        min-height: 46px;
        font-weight: 600;
    }

    .stButton > button:hover {
        background: linear-gradient(90deg, #9a52ca, #d267e0);
        color: white !important;
    }

    /* ================================
       DATAFRAME
    ================================= */

    [data-testid="stDataFrame"] {
        border: 1px solid #60427f;
        border-radius: 12px;
        overflow: hidden;
    }

    /* ================================
       DIVIDER
    ================================= */

    hr {
        border-color: #4b3465 !important;
    }

    /* ================================
       HERO BANNER
    ================================= */

    .oncogenex-hero {
        position: relative;
        overflow: hidden;
        min-height: 420px;
        border-radius: 28px;
        padding: 3.2rem 3.4rem;
        margin: 0.2rem 0 2.2rem 0;
        background:
            radial-gradient(circle at 72% 35%, rgba(179, 74, 222, 0.30), transparent 18%),
            radial-gradient(circle at 90% 75%, rgba(224, 54, 161, 0.18), transparent 25%),
            linear-gradient(135deg, #160b29 0%, #2a1245 52%, #1b0d31 100%);
        border: 1px solid rgba(192, 104, 235, 0.78);
        box-shadow: 0 25px 70px rgba(0, 0, 0, 0.52), inset 0 1px 0 rgba(255,255,255,0.16), inset 0 0 55px rgba(189, 83, 227, 0.09);
    }

    .oncogenex-hero:before {
        content: "";
        position: absolute;
        inset: 10px;
        border-radius: 22px;
        border: 1px solid rgba(232, 175, 255, 0.18);
        pointer-events: none;
    }

    .hero-grid {
        position: relative;
        z-index: 2;
        display: grid;
        grid-template-columns: minmax(0, 1.05fr) minmax(280px, 0.95fr);
        gap: 2rem;
        align-items: center;
    }

    .hero-brand {
        display: flex;
        align-items: center;
        gap: 1rem;
        margin-bottom: 0.6rem;
    }

    .hero-dna {
        font-size: 5.4rem;
        filter: drop-shadow(0 8px 14px rgba(205, 82, 235, 0.45));
        transform: rotate(-6deg);
    }

    .hero-name {
        font-size: clamp(3rem, 6vw, 6rem);
        font-weight: 900;
        letter-spacing: -0.06em;
        line-height: 0.95;
        margin: 0;
        color: #fbf8ff !important;
        text-shadow: 0 4px 0 #442263, 0 13px 25px rgba(0,0,0,0.48);
    }

    .hero-name .x {
        color: #c86cff !important;
        text-shadow: 0 4px 0 #74359a, 0 13px 25px rgba(184, 78, 225, 0.38);
    }

    .hero-tagline {
        font-size: 1.65rem;
        font-weight: 700;
        color: #eadff3 !important;
        margin: 0.7rem 0 1.7rem 6.6rem;
        padding-bottom: 1rem;
        border-bottom: 2px solid rgba(197, 101, 236, 0.78);
    }

    .hero-message, .hero-build {
        background: rgba(20, 9, 36, 0.68);
        border: 1px solid rgba(190, 99, 230, 0.54);
        border-radius: 18px;
        box-shadow: 0 14px 30px rgba(0,0,0,0.32), inset 0 1px 0 rgba(255,255,255,0.07);
    }

    .hero-message {
        display: flex;
        gap: 1rem;
        align-items: center;
        padding: 1.2rem 1.4rem;
        max-width: 760px;
        font-size: 1.22rem;
        line-height: 1.55;
        color: #ead9ff !important;
    }

    .hero-message-icon { font-size: 2.8rem; }

    .hero-build {
        margin-top: 1rem;
        padding: 1rem 1.4rem;
        font-size: 1.05rem;
        color: #d8c9e7 !important;
    }

    .hero-build strong {
        color: #d78cff !important;
        font-size: 1.15rem;
    }

    .hero-build .sep {
        color: #ff71bd !important;
        margin: 0 0.8rem;
    }

    .hero-visual {
        position: relative;
        min-height: 350px;
        display: flex;
        align-items: center;
        justify-content: center;
    }

    .dna-3d {
        width: min(100%, 520px);
        height: 350px;
        display: flex;
        align-items: center;
        justify-content: center;
        transform: rotate(12deg) scaleX(0.92);
        filter: drop-shadow(0 0 16px rgba(197, 92, 255, 0.50)) drop-shadow(18px 20px 20px rgba(0,0,0,0.52));
        animation: floatDNA 5s ease-in-out infinite;
        z-index: 2;
    }

    .dna-3d svg {
        width: 100%;
        height: 100%;
        overflow: visible;
    }

    .dna-3d .dna-rail-a,
    .dna-3d .dna-rail-b {
        fill: none;
        stroke-linecap: round;
        stroke-width: 18;
    }

    .dna-3d .dna-rung {
        stroke: rgba(241, 224, 255, 0.78);
        stroke-width: 7;
        stroke-linecap: round;
    }

    .molecule {
        position: absolute;
        width: 20px;
        height: 20px;
        border-radius: 50%;
        background: #d374ff;
        box-shadow: 0 0 24px #a54be3;
        animation: pulseMolecule 2.8s ease-in-out infinite;
    }

    .m1 { top: 8%; left: 15%; }
    .m2 { top: 22%; right: 10%; animation-delay: .5s; }
    .m3 { bottom: 20%; left: 7%; animation-delay: 1s; }
    .m4 { bottom: 8%; right: 22%; animation-delay: 1.5s; }

    @keyframes floatDNA {
        0%, 100% { transform: rotate(12deg) scaleX(0.82) translateY(0); }
        50% { transform: rotate(8deg) scaleX(0.82) translateY(-10px); }
    }

    @keyframes pulseMolecule {
        0%, 100% { transform: scale(0.8); opacity: 0.55; }
        50% { transform: scale(1.35); opacity: 1; }
    }

    @media (max-width: 900px) {
        .oncogenex-hero { padding: 2rem 1.5rem; min-height: auto; }
        .hero-grid { grid-template-columns: 1fr; }
        .hero-tagline { margin-left: 0; font-size: 1.25rem; }
        .hero-visual { min-height: 250px; }
        .dna-3d { font-size: 11rem; }
    }



    /* ================================
       MOBILE RESPONSIVE FIX
    ================================= */
    html, body {
        max-width: 100%;
        overflow-x: hidden;
    }

    .stApp,
    [data-testid="stAppViewContainer"],
    [data-testid="stMain"] {
        max-width: 100%;
        overflow-x: hidden;
    }

    *, *::before, *::after {
        box-sizing: border-box;
    }

    @media (max-width: 768px) {
        .main .block-container {
            width: 100%;
            max-width: 100%;
            padding: 1rem 0.75rem 2rem 0.75rem;
        }

        .oncogenex-hero {
            width: 100%;
            max-width: 100%;
            min-height: auto;
            padding: 1.25rem;
            margin: 0 0 1.25rem 0;
            border-radius: 20px;
        }

        .oncogenex-hero:before {
            inset: 7px;
            border-radius: 14px;
        }

        .hero-grid {
            width: 100%;
            max-width: 100%;
            grid-template-columns: minmax(0, 1fr);
            gap: 1.25rem;
        }

        .hero-brand {
            flex-wrap: wrap;
            gap: 0.75rem;
        }

        .hero-dna {
            font-size: 3.5rem;
        }

        .hero-name {
            font-size: clamp(2.2rem, 13vw, 3.5rem);
            letter-spacing: -0.04em;
            overflow-wrap: anywhere;
        }

        .hero-tagline {
            margin: 0.75rem 0 1rem 0;
            font-size: 1.05rem;
            line-height: 1.4;
            padding-bottom: 0.75rem;
            overflow-wrap: anywhere;
        }

        .hero-message {
            width: 100%;
            max-width: 100%;
            flex-direction: column;
            align-items: flex-start;
            gap: 0.6rem;
            padding: 1rem;
            font-size: 1rem;
            line-height: 1.5;
            overflow-wrap: anywhere;
        }

        .hero-message-icon {
            font-size: 2.2rem;
        }

        .hero-build {
            width: 100%;
            max-width: 100%;
            padding: 0.9rem 1rem;
            font-size: 0.95rem;
            line-height: 1.7;
            overflow-wrap: anywhere;
        }

        .hero-build strong {
            font-size: 1rem;
        }

        .hero-build .sep {
            display: inline-block;
            margin: 0 0.35rem;
        }

        .hero-visual {
            width: 100%;
            max-width: 100%;
            min-height: 210px;
        }

        .dna-3d {
            width: min(100%, 300px);
            height: 210px;
            transform: rotate(8deg) scaleX(0.9);
        }

        /* Stack Streamlit column layouts vertically on phones.
           This prevents 4 desktop columns from becoming extremely narrow. */
        [data-testid="stHorizontalBlock"] {
            flex-direction: column !important;
            flex-wrap: nowrap !important;
            width: 100% !important;
            gap: 1rem !important;
        }

        [data-testid="stHorizontalBlock"] > div,
        [data-testid="column"] {
            width: 100% !important;
            min-width: 0 !important;
            max-width: 100% !important;
            flex: 1 1 100% !important;
        }

        .dashboard-card,
        .gene-card,
        .info-box,
        .status-significant,
        .status-up,
        .status-down,
        .status-not-significant {
            width: 100%;
            max-width: 100%;
            overflow-wrap: anywhere;
        }
    }

    /* ================================
       FOOTER
    ================================= */

    .footer {
        text-align: center;
        color: #ad91c7 !important;
        padding: 2rem;
        font-size: 0.9rem;
    }

</style>
""", unsafe_allow_html=True)


# ============================================================
# HELPER FUNCTIONS
# ============================================================

@st.cache_data
def load_csv(path):
    try:
        return pd.read_csv(path)
    except Exception:
        return None


def find_first_existing(paths):
    for path in paths:
        if os.path.exists(path):
            return path
    return None


def find_file(patterns):
    for pattern in patterns:
        files = glob.glob(pattern, recursive=True)
        if files:
            return files[0]
    return None


def find_column(df, possible_names):
    if df is None:
        return None

    normalized_columns = {
        str(col).strip().lower(): col
        for col in df.columns
    }

    for name in possible_names:
        name_lower = name.strip().lower()

        if name_lower in normalized_columns:
            return normalized_columns[name_lower]

    return None


def get_gene_symbol_column(df):

    return find_column(df, [
        "Gene Symbol",
        "Gene_Symbol",
        "gene_symbol",
        "GeneSymbol",
        "Symbol",
        "SYMBOL",
        "symbol",
        "Gene"
    ])


def get_probe_column(df):

    return find_column(df, [
        "Gene",
        "Probe",
        "Probe_ID",
        "Probe ID",
        "ID",
        "ID_REF"
    ])


def get_logfc_column(df):

    return find_column(df, [
        "Log2FC",
        "log2fc",
        "logFC",
        "LogFC",
        "log_fold_change",
        "Log2 Fold Change"
    ])


def get_pvalue_column(df):

    return find_column(df, [
        "Adjusted_P_Value",
        "Adjusted P Value",
        "Adjusted_Pvalue",
        "adj.P.Val",
        "adj_p_value",
        "padj",
        "FDR",
        "P_adj"
    ])


def get_mean_cancer_column(df):

    return find_column(df, [
        "Mean_Cancer",
        "Mean Cancer",
        "Cancer_Mean",
        "Cancer Expression"
    ])


def get_mean_normal_column(df):

    return find_column(df, [
        "Mean_Normal",
        "Mean Normal",
        "Normal_Mean",
        "Normal Expression"
    ])


def clean_gene_query(gene):

    if gene is None:
        return ""

    return str(gene).strip().upper()


def search_gene_records(df, query):
    """Search gene symbol, probe ID, display name and gene title safely."""
    if df is None or df.empty:
        return pd.DataFrame()

    query = clean_gene_query(query)
    if not query:
        return pd.DataFrame()

    candidate_columns = [
        "Gene Symbol", "Gene_Symbol", "gene_symbol", "GeneSymbol",
        "Symbol", "SYMBOL", "symbol", "Display_Name", "Gene",
        "Probe", "Probe_ID", "Probe ID", "Gene Title", "Title"
    ]

    matches = []
    seen_columns = set()

    for name in candidate_columns:
        col = find_column(df, [name])
        if col is None or col in seen_columns:
            continue
        seen_columns.add(col)
        series = df[col].fillna("").astype(str).str.strip().str.upper()
        exact = df[series == query]
        if not exact.empty:
            matches.append(exact)

    if not matches:
        # Fallback: partial search across likely annotation columns.
        for col in df.columns:
            col_name = str(col).lower()
            if any(key in col_name for key in ["gene", "symbol", "probe", "title", "display"]):
                series = df[col].fillna("").astype(str).str.upper()
                hit = df[series.str.contains(query, regex=False, na=False)]
                if not hit.empty:
                    matches.append(hit)

    if not matches:
        return pd.DataFrame(columns=df.columns)

    return pd.concat(matches, ignore_index=True).drop_duplicates()


# ============================================================
# LOAD ANALYSIS DATA
# ============================================================

# Prefer the annotated results because gene-symbol searches (e.g. BRCA1)
# require the annotation columns.
ALL_RESULTS_PATH = find_first_existing([
    "results/all_differential_expression_results_annotated.csv",
    "results/significant_genes_with_symbols.csv",
    "results/all_differential_expression_results.csv"
])

if ALL_RESULTS_PATH is None:
    ALL_RESULTS_PATH = find_file([
        "results/**/*differential*expression*.csv",
        "results/**/*results*.csv"
    ])


UP_PATH = find_first_existing([
    "results/upregulated_genes.csv"
])

DOWN_PATH = find_first_existing([
    "results/downregulated_genes.csv"
])

ALL_DF = load_csv(ALL_RESULTS_PATH) if ALL_RESULTS_PATH else None
UP_DF = load_csv(UP_PATH) if UP_PATH else None
DOWN_DF = load_csv(DOWN_PATH) if DOWN_PATH else None


# ============================================================
# SESSION STATE
# ============================================================

if "page" not in st.session_state:
    st.session_state.page = "Dashboard"

if "gene_query" not in st.session_state:
    st.session_state.gene_query = ""

if "search_results" not in st.session_state:
    st.session_state.search_results = pd.DataFrame()


# ============================================================
# SIDEBAR
# ============================================================

with st.sidebar:

    st.html("""
    <div style="padding-top: 1rem;">
        <h1 style="
            font-size: 2rem;
            margin-bottom: 0.3rem;
        ">
            🧬 OncoGeneX
        </h1>

        <p style="
            color: #9fc5dd;
            font-size: 0.95rem;
        ">
            Cancer Gene Expression Analysis Platform
        </p>
    </div>
    """)

    st.divider()

    st.markdown("### Navigation")

    pages = [
        "🏠 Dashboard",
        "🔎 Gene Explorer",
        "📊 Differential Expression",
        "📈 Visualizations",
        "🧪 Pathway Analysis",
        "🔗 PPI & Hub Genes",
        "🤖 Machine Learning"
    ]

    selected_page = st.radio(
        "Navigation",
        pages,
        label_visibility="collapsed"
    )

    st.divider()

    st.markdown("### Dataset")

    st.write(DATASET_ID)
    st.write(CANCER_TYPE)
    st.write("Cancer vs Normal")


# Remove emoji for internal logic
PAGE = selected_page.split(" ", 1)[1]


# ============================================================
# HEADER
# ============================================================

if PAGE != "Dashboard":
    st.markdown("""
    <div>
        <div class="main-title">🧬 OncoGeneX</div>
        <div class="subtitle">Interactive Cancer Gene Expression Analysis Platform</div>
    </div>
    """, unsafe_allow_html=True)


# ============================================================
# DASHBOARD
# ============================================================

if PAGE == "Dashboard":

    # HERO BANNER — rendered first on the dashboard
    st.html("""
    <section class="oncogenex-hero">
        <div class="hero-grid">
            <div>
                <div class="hero-brand">
                    <div class="hero-dna">🧬</div>
                    <h1 class="hero-name">OncoGene<span class="x">X</span></h1>
                </div>
                <div class="hero-tagline">Cancer Gene Expression Analysis Platform</div>
                <div class="hero-message">
                    <div class="hero-message-icon">🔬</div>
                    <div>Exploring differential gene expression, cancer versus normal samples, biological pathways and molecular networks.</div>
                </div>
                <div class="hero-build">
                    Built by <strong>Mansoor</strong><span class="sep">|</span><strong>Kalindhy</strong><span class="sep">|</span><strong>Heama shri</strong><span class="sep">|</span><strong>Kireeti</strong>
                </div>
            </div>
            <div class="hero-visual" aria-label="DNA analysis visual">
                <span class="molecule m1"></span><span class="molecule m2"></span>
                <span class="molecule m3"></span><span class="molecule m4"></span>
                <div class="dna-3d" aria-hidden="true">
                    <svg viewBox="0 0 520 360" role="img" aria-label="Professional DNA double helix">
                        <defs>
                            <linearGradient id="railPink" x1="0" y1="0" x2="1" y2="1">
                                <stop offset="0%" stop-color="#ff4b8b"/>
                                <stop offset="50%" stop-color="#e51d68"/>
                                <stop offset="100%" stop-color="#9d174d"/>
                            </linearGradient>
                            <linearGradient id="railBlue" x1="1" y1="0" x2="0" y2="1">
                                <stop offset="0%" stop-color="#73e8ff"/>
                                <stop offset="48%" stop-color="#24b9ed"/>
                                <stop offset="100%" stop-color="#2563eb"/>
                            </linearGradient>
                            <linearGradient id="rungGlow" x1="0" y1="0" x2="1" y2="0">
                                <stop offset="0%" stop-color="#ff7ab1"/>
                                <stop offset="50%" stop-color="#d8f7ff"/>
                                <stop offset="100%" stop-color="#42d8ff"/>
                            </linearGradient>
                            <filter id="dnaGlow" x="-30%" y="-30%" width="160%" height="160%">
                                <feGaussianBlur stdDeviation="7" result="blur"/>
                                <feMerge>
                                    <feMergeNode in="blur"/>
                                    <feMergeNode in="SourceGraphic"/>
                                </feMerge>
                            </filter>
                        </defs>

                        <g filter="url(#dnaGlow)">
                            <path class="dna-rail-a" stroke="url(#railPink)"
                                d="M135,40 C315,85 320,125 180,170 C40,215 55,265 275,320"/>
                            <path class="dna-rail-b" stroke="url(#railBlue)"
                                d="M385,40 C205,85 200,125 340,170 C480,215 465,265 245,320"/>

                            <g stroke="url(#rungGlow)">
                                <line class="dna-rung" x1="154" y1="52" x2="366" y2="52"/>
                                <line class="dna-rung" x1="211" y1="76" x2="309" y2="76"/>
                                <line class="dna-rung" x1="266" y1="103" x2="254" y2="103"/>
                                <line class="dna-rung" x1="292" y1="130" x2="228" y2="130"/>
                                <line class="dna-rung" x1="257" y1="151" x2="263" y2="151"/>
                                <line class="dna-rung" x1="182" y1="178" x2="338" y2="178"/>
                                <line class="dna-rung" x1="115" y1="205" x2="405" y2="205"/>
                                <line class="dna-rung" x1="103" y1="232" x2="417" y2="232"/>
                                <line class="dna-rung" x1="149" y1="259" x2="371" y2="259"/>
                                <line class="dna-rung" x1="220" y1="287" x2="300" y2="287"/>
                                <line class="dna-rung" x1="266" y1="310" x2="254" y2="310"/>
                            </g>
                        </g>
                    </svg>
                </div>
            </div>
        </div>
    </section>
    """)

    # ABOUT THE WEBSITE — directly after the banner
    st.header("🔬 About OncoGeneX")
    st.markdown("""
    <div class="info-box">
    OncoGeneX is an interactive bioinformatics platform for comparing gene-expression patterns between cancer and normal samples.
    It supports differential-expression analysis, biological pathway exploration, protein-protein interaction and hub-gene analysis,
    visualisation, and research-oriented machine-learning results.
    </div>
    """, unsafe_allow_html=True)

    st.divider()

    # SEARCH OPTION — after About section
    st.header("🔎 Search a Gene")
    st.write("Enter a gene symbol to quickly search the analysed dataset.")

    dash_query = st.text_input(
        "Gene symbol",
        value=st.session_state.gene_query,
        placeholder="Examples: TP53, BRCA1, EGFR, KRAS",
        key="dashboard_gene_search"
    )

    search_col1, search_col2, search_col3 = st.columns([1, 2, 1])
    with search_col2:
        dash_search = st.button("🔍 Search Gene", width="stretch", key="dashboard_search_button")

    if dash_search:
        query = clean_gene_query(dash_query)
        st.session_state.gene_query = query
        st.session_state.search_results = search_gene_records(ALL_DF, query)
        if not query:
            st.warning("Please enter a gene symbol.")
        elif st.session_state.search_results.empty:
            st.warning(f"{query} was not found in the analysed dataset.")
        else:
            st.success(f"Found {len(st.session_state.search_results)} matching analysed record(s) for {query}.")
            st.dataframe(st.session_state.search_results, width="stretch")
            st.info("Open Gene Explorer from the sidebar for detailed analysis of the matching record(s).")

    st.divider()

    # SUMMARY CARDS
    total_records = len(ALL_DF) if ALL_DF is not None else 0
    up_count = len(UP_DF) if UP_DF is not None else 0
    down_count = len(DOWN_DF) if DOWN_DF is not None else 0

    st.subheader("📊 Analysis Summary")
    col1, col2, col3, col4 = st.columns(4)

    for col, label, value in [
        (col1, "Analysed Records", f"{total_records:,}"),
        (col2, "Upregulated DEGs", f"{up_count:,}"),
        (col3, "Downregulated DEGs", f"{down_count:,}"),
        (col4, "Dataset", DATASET_ID),
    ]:
        with col:
            st.markdown(f'<div class="dashboard-card"><div class="gene-card-label">{label}</div><div class="gene-card-value">{value}</div></div>', unsafe_allow_html=True)

    st.divider()
    st.subheader("⚙️ Analysis Workflow")
    workflow = [
        ("1️⃣ Data Processing", "Gene-expression data is prepared and analysed."),
        ("2️⃣ Differential Expression", "Cancer and normal samples are statistically compared."),
        ("3️⃣ Biological Analysis", "Pathways and molecular networks are investigated."),
        ("4️⃣ Prediction", "Machine-learning models evaluate cancer classification."),
    ]
    for col, (title, text) in zip(st.columns(4), workflow):
        with col:
            st.markdown(f"""<div class="dashboard-card"><b>{title}</b><br><br>{text}</div>""", unsafe_allow_html=True)


# ============================================================
# GENE EXPLORER
# ============================================================

elif PAGE == "Gene Explorer":

    st.header("🔎 Gene Explorer")

    st.write(
        "Search for a gene symbol to examine its expression pattern "
        "in the cancer versus normal analysis."
    )

    gene_input = st.text_input(
        "Enter a gene symbol",
        value=st.session_state.gene_query,
        placeholder="Examples: TP53, BRCA1, EGFR, KRAS"
    )

    search_col1, search_col2, search_col3 = st.columns([1, 2, 1])

    with search_col2:
        analyze_button = st.button(
            "🔍 Analyze Gene",
            width="stretch"
        )

    if analyze_button:

        query = clean_gene_query(gene_input)
        st.session_state.gene_query = query

        if query == "":
            st.warning("Please enter a gene symbol.")

        elif ALL_DF is None:
            st.error(
                "The complete differential-expression dataset could not be found."
            )

        else:

            gene_column = get_gene_symbol_column(ALL_DF)

            if gene_column is None:

                st.error(
                    "A gene-symbol column could not be identified in the dataset."
                )

                st.write("Available columns:")
                st.write(list(ALL_DF.columns))

            else:

                matches = search_gene_records(ALL_DF, query)
                st.session_state.search_results = matches

    # --------------------------------------------------------
    # DISPLAY RESULTS
    # --------------------------------------------------------

    results = st.session_state.search_results

    if (
        st.session_state.gene_query != ""
        and results is not None
        and not results.empty
    ):

        st.divider()

        gene = st.session_state.gene_query

        st.success(
            f"Found {len(results)} matching analysed record(s)."
        )

        if len(results) > 1:

            st.markdown("""
            <div class="info-box">
            Multiple microarray probes were found for this gene.
            Select a probe below for detailed analysis.
            </div>
            """, unsafe_allow_html=True)

        probe_column = get_probe_column(results)

        if probe_column is not None:

            probe_options = results[probe_column].astype(str).tolist()

            selected_probe = st.selectbox(
                "Select probe",
                probe_options
            )

            selected_record = results[
                results[probe_column].astype(str) == selected_probe
            ].iloc[0]

        else:
            selected_record = results.iloc[0]

        # ----------------------------------------------------
        # EXTRACT VALUES
        # ----------------------------------------------------

        logfc_col = get_logfc_column(results)
        pvalue_col = get_pvalue_column(results)
        cancer_col = get_mean_cancer_column(results)
        normal_col = get_mean_normal_column(results)

        logfc = None
        pvalue = None
        cancer_mean = None
        normal_mean = None

        if logfc_col:
            logfc = pd.to_numeric(
                selected_record[logfc_col],
                errors="coerce"
            )

        if pvalue_col:
            pvalue = pd.to_numeric(
                selected_record[pvalue_col],
                errors="coerce"
            )

        if cancer_col:
            cancer_mean = pd.to_numeric(
                selected_record[cancer_col],
                errors="coerce"
            )

        if normal_col:
            normal_mean = pd.to_numeric(
                selected_record[normal_col],
                errors="coerce"
            )

        # ----------------------------------------------------
        # GENE HEADER
        # ----------------------------------------------------

        st.markdown(
            f"## 🧬 {gene}"
        )

        if probe_column:
            st.write(
                f"**Probe ID:** {selected_record[probe_column]}"
            )

        # ----------------------------------------------------
        # SIGNIFICANCE CALCULATION
        # ----------------------------------------------------

        significant = False
        regulation = "Not significant"

        if (
            logfc is not None
            and pvalue is not None
            and not pd.isna(logfc)
            and not pd.isna(pvalue)
        ):

            significant = (
                pvalue < P_VALUE_THRESHOLD
                and abs(logfc) >= LOG2FC_THRESHOLD
            )

            if significant:

                if logfc > 0:
                    regulation = "Upregulated"
                else:
                    regulation = "Downregulated"

        # ----------------------------------------------------
        # STATUS
        # ----------------------------------------------------

        if significant and regulation == "Upregulated":

            st.markdown("""
            <div class="status-up">
            🟢 SIGNIFICANTLY UPREGULATED IN CANCER
            </div>
            """, unsafe_allow_html=True)

        elif significant and regulation == "Downregulated":

            st.markdown("""
            <div class="status-down">
            🔴 SIGNIFICANTLY DOWNREGULATED IN CANCER
            </div>
            """, unsafe_allow_html=True)

        else:

            st.markdown("""
            <div class="status-not-significant">
            🟡 NOT SIGNIFICANTLY DIFFERENTIALLY EXPRESSED
            </div>
            """, unsafe_allow_html=True)

        st.markdown("## 📊 Expression Comparison")

        # ----------------------------------------------------
        # RESPONSIVE METRIC CARDS
        # ----------------------------------------------------

        metric_col1, metric_col2 = st.columns(2)

        with metric_col1:

            cancer_display = (
                f"{cancer_mean:.4f}"
                if cancer_mean is not None
                and not pd.isna(cancer_mean)
                else "N/A"
            )

            st.markdown(f'''<div class="gene-card">
<div class="gene-card-label">Mean Cancer Expression</div>
<div class="gene-card-value">{cancer_display}</div>
</div>''', unsafe_allow_html=True)

        with metric_col2:

            normal_display = (
                f"{normal_mean:.4f}"
                if normal_mean is not None
                and not pd.isna(normal_mean)
                else "N/A"
            )

            st.markdown(f'''<div class="gene-card">
<div class="gene-card-label">Mean Normal Expression</div>
<div class="gene-card-value">{normal_display}</div>
</div>''', unsafe_allow_html=True)

        st.write("")

        metric_col3, metric_col4 = st.columns(2)

        with metric_col3:

            logfc_display = (
                f"{logfc:.4f}"
                if logfc is not None
                and not pd.isna(logfc)
                else "N/A"
            )

            st.markdown(f'''<div class="gene-card">
<div class="gene-card-label">Log2 Fold Change</div>
<div class="gene-card-value">{logfc_display}</div>
</div>''', unsafe_allow_html=True)

        with metric_col4:

            pvalue_display = (
                f"{pvalue:.2e}"
                if pvalue is not None
                and not pd.isna(pvalue)
                else "N/A"
            )

            st.markdown(f'''<div class="gene-card">
<div class="gene-card-label">Adjusted P-value</div>
<div class="gene-card-value">{pvalue_display}</div>
</div>''', unsafe_allow_html=True)

        # ----------------------------------------------------
        # INTERPRETATION
        # ----------------------------------------------------

        st.markdown("## 🔬 Analysis Interpretation")

        if significant:

            if regulation == "Upregulated":

                message = (
                    f"{gene} showed statistically significant higher "
                    f"expression in cancer compared with normal samples."
                )

            else:

                message = (
                    f"{gene} showed statistically significant lower "
                    f"expression in cancer compared with normal samples."
                )

            st.markdown(
                f"""
                <div class="info-box">
                {message}
                </div>
                """,
                unsafe_allow_html=True
            )

        else:

            st.markdown(
                f"""
                <div class="info-box">
                {gene} was detected in the dataset but did not meet
                the project's defined differential-expression
                significance criteria.
                </div>
                """,
                unsafe_allow_html=True
            )

        # ----------------------------------------------------
        # COMPLETE RECORD
        # ----------------------------------------------------

        with st.expander("View Complete Analysis Record"):

            complete_record = (
                selected_record
                .to_frame()
                .T
            )

            st.dataframe(
                complete_record,
                width="stretch"
            )

    # --------------------------------------------------------
    # NO RESULTS
    # --------------------------------------------------------

    elif (
        st.session_state.gene_query != ""
        and results is not None
        and results.empty
    ):

        st.divider()

        st.warning(
            f"{st.session_state.gene_query} was not found "
            f"in the analysed dataset."
        )

        st.markdown("""
        <div class="info-box">
        This does not necessarily mean the gene is biologically absent.
        It may not be represented by a matching probe or annotation in
        the analysed microarray dataset.
        </div>
        """, unsafe_allow_html=True)


# ============================================================
# DIFFERENTIAL EXPRESSION
# ============================================================

elif PAGE == "Differential Expression":

    st.header("📊 Differential Expression Analysis")

    st.write(
        "Overview of genes identified through the cancer versus "
        "normal differential-expression analysis."
    )

    if ALL_DF is not None:

        col1, col2, col3 = st.columns(3)

        with col1:
            st.metric(
                "Total Analysed Records",
                f"{len(ALL_DF):,}"
            )

        with col2:
            st.metric(
                "Upregulated DEGs",
                f"{len(UP_DF):,}" if UP_DF is not None else "N/A"
            )

        with col3:
            st.metric(
                "Downregulated DEGs",
                f"{len(DOWN_DF):,}" if DOWN_DF is not None else "N/A"
            )

        st.divider()

        view_option = st.selectbox(
            "Select result dataset",
            [
                "Complete Differential Expression Results",
                "Upregulated Genes",
                "Downregulated Genes"
            ]
        )

        if view_option == "Complete Differential Expression Results":

            st.dataframe(
                ALL_DF,
                width="stretch",
                height=500
            )

        elif view_option == "Upregulated Genes":

            if UP_DF is not None:

                st.dataframe(
                    UP_DF,
                    width="stretch",
                    height=500
                )

            else:
                st.warning("Upregulated gene results were not found.")

        elif view_option == "Downregulated Genes":

            if DOWN_DF is not None:

                st.dataframe(
                    DOWN_DF,
                    width="stretch",
                    height=500
                )

            else:
                st.warning("Downregulated gene results were not found.")

    else:

        st.error(
            "Differential-expression results could not be loaded."
        )


# ============================================================
# VISUALIZATIONS
# ============================================================

elif PAGE == "Visualizations":

    st.header("📈 Analysis Visualizations")

    visualization_files = {
        "PCA Visualization":
            "results/visualization/pca_visualization.png",

        "Top 20 Gene Heatmap":
            "results/visualization/top_20_gene_heatmap.png",

        "Clustered Top 20 Gene Heatmap":
            "results/clustered_top_20_gene_heatmap.png",

        "Feature Importance":
            "results/top_20_feature_importance.png",

        "ROC Curve":
            "results/roc_curve.png",

        "Confusion Matrix":
            "results/confusion_matrix.png"
    }

    available = {
        name: path
        for name, path in visualization_files.items()
        if os.path.exists(path)
    }

    if available:

        selected_visualization = st.selectbox(
            "Select visualization",
            list(available.keys())
        )

        st.image(
            available[selected_visualization],
            width="stretch"
        )

    else:

        st.warning(
            "No visualization images were found."
        )


# ============================================================
# PATHWAY ANALYSIS
# ============================================================

elif PAGE == "Pathway Analysis":

    st.header("🧪 Pathway Enrichment Analysis")

    enrichment_files = {
        "GO Biological Process — Upregulated":
            "results/enrichment/go_biological_process_upregulated_barplot.png",

        "GO Biological Process — Downregulated":
            "results/enrichment/go_biological_process_downregulated_barplot.png",

        "KEGG Pathway — Upregulated":
            "results/enrichment/kegg_pathway_upregulated_barplot.png",

        "KEGG Pathway — Downregulated":
            "results/enrichment/kegg_pathway_downregulated_barplot.png",

        "Upregulated Pathways":
            "results/enrichment/upregulated_pathway_barplot.png",

        "Downregulated Pathways":
            "results/enrichment/downregulated_pathway_barplot.png"
    }

    available = {
        name: path
        for name, path in enrichment_files.items()
        if os.path.exists(path)
    }

    if available:

        selected = st.selectbox(
            "Select enrichment analysis",
            list(available.keys())
        )

        st.image(
            available[selected],
            width="stretch"
        )

    else:

        st.warning(
            "Pathway enrichment visualizations were not found."
        )


# ============================================================
# PPI AND HUB GENES
# ============================================================

elif PAGE == "PPI & Hub Genes":

    st.header("🔗 PPI Network & Hub Gene Analysis")

    analysis_type = st.selectbox(
        "Select analysis",
        [
            "Upregulated PPI Network",
            "Downregulated PPI Network",
            "Upregulated Hub Genes",
            "Downregulated Hub Genes"
        ]
    )

    if analysis_type == "Upregulated PPI Network":

        possible_paths = [
            "results/ppi/upregulated_ppi_network.png",
            "results/ppi/upregulated_network.png"
        ]

        path = find_first_existing(possible_paths)

        if path:
            st.image(
                path,
                width="stretch"
            )
        else:
            st.warning("Upregulated PPI network not found.")

    elif analysis_type == "Downregulated PPI Network":

        possible_paths = [
            "results/ppi/downregulated_ppi_network.png",
            "results/ppi/downregulated_network.png"
        ]

        path = find_first_existing(possible_paths)

        if path:
            st.image(
                path,
                width="stretch"
            )
        else:
            st.warning("Downregulated PPI network not found.")

    elif analysis_type == "Upregulated Hub Genes":

        path = find_first_existing([
            "results/ppi/upregulated_hub_genes.csv"
        ])

        if path:

            hub_df = load_csv(path)

            st.dataframe(
                hub_df,
                width="stretch"
            )

        else:
            st.warning("Upregulated hub genes not found.")

    elif analysis_type == "Downregulated Hub Genes":

        path = find_first_existing([
            "results/ppi/downregulated_hub_genes.csv"
        ])

        if path:

            hub_df = load_csv(path)

            st.dataframe(
                hub_df,
                width="stretch"
            )

        else:
            st.warning("Downregulated hub genes not found.")


# ============================================================
# MACHINE LEARNING
# ============================================================

elif PAGE == "Machine Learning":

    st.header("🤖 Machine Learning Analysis")

    st.markdown("""
    <div class="info-box">
    The machine-learning component uses selected gene-expression
    features to distinguish cancer samples from normal samples.
    </div>
    """, unsafe_allow_html=True)

    ml_images = {
        "ROC Curve":
            "results/roc_curve.png",

        "Confusion Matrix":
            "results/confusion_matrix.png",

        "Feature Importance":
            "results/top_20_feature_importance.png"
    }

    available_ml = {
        name: path
        for name, path in ml_images.items()
        if os.path.exists(path)
    }

    if available_ml:

        selected_ml = st.selectbox(
            "Select machine-learning result",
            list(available_ml.keys())
        )

        st.image(
            available_ml[selected_ml],
            width="stretch"
        )

    else:

        st.warning(
            "Machine-learning visualization files were not found."
        )


# ============================================================
# FOOTER
# ============================================================

st.divider()

st.markdown("""
<div class="footer">
🧬 <b>OncoGeneX</b> |
Cancer Gene Expression Analysis Platform |
Cancer vs Normal Differential Expression Analysis
</div>
""", unsafe_allow_html=True)