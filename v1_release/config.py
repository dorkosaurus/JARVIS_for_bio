"""Project configuration: paths, constants, hyperparameters."""

from pathlib import Path

# --- Paths ---
PROJECT_ROOT = Path(__file__).resolve().parent
DATA_DIR = PROJECT_ROOT / "data"
SEQUENCES_DIR = DATA_DIR / "sequences"
EMBEDDINGS_DIR = DATA_DIR / "embeddings" / "esm3"
EMBEDDINGS_INDEX = EMBEDDINGS_DIR / "index.parquet"
FEATURES_DIR = DATA_DIR / "features"
PUBLIC_DATA_DIR = DATA_DIR / "public"
RESULTS_DB = DATA_DIR / "results.db"
PRETRAINED_POLICY = DATA_DIR / "pretrained_policy.pt"
OUTPUTS_DIR = PROJECT_ROOT / "outputs"

# --- AAV2 wildtype (NCBI Protein) ---
AAV2_ACCESSION = "YP_680426.1"
AAV2_FASTA = SEQUENCES_DIR / "aav2_vp1_wildtype.fasta"

# --- AAV.7m8 reference (Dalkara 2013) ---
# AAV2 backbone with LALGETTRP inserted between residues 587 and 588.
AAV7M8_FASTA = SEQUENCES_DIR / "aav7m8_vp1_reference.fasta"
AAV7M8_PEPTIDE = "LALGETTRP"
AAV7M8_INSERTION_AFTER = 587  # 1-indexed residue position

# --- Capsid variants library ---
CAPSID_VARIANTS_FASTA = SEQUENCES_DIR / "capsid_variants.fasta"

# --- VR loops (AAV2 VP1, 1-indexed inclusive ranges) ---
# Govindasamy 2006 / DiMattia 2012 assignment.
VR_LOOPS = {
    "VR-I":    (263, 268),
    "VR-IV":   (449, 468),
    "VR-V":    (488, 505),
    "VR-VIII": (581, 593),
    "VR-IX":   (704, 714),
}

# Loops eligible for substitution variants
SUBSTITUTION_LOOPS = ["VR-IV", "VR-V", "VR-VIII", "VR-IX"]

# 7-mer insertion site (AAV.7m8 trick)
INSERTION_AFTER_RESIDUE = 587  # 1-indexed; insert between 587 and 588

# --- Variant generation ---
N_SUBSTITUTION_VARIANTS = 50
N_INSERTION_VARIANTS = 30
N_STACKED_VARIANTS = 40        # 7-mer insertion + 2-4 substitutions at non-VIII VR loops
MIN_SUBSTITUTIONS_PER_VARIANT = 1
MAX_SUBSTITUTIONS_PER_VARIANT = 6
MIN_STACKED_SUBSTITUTIONS = 2
MAX_STACKED_SUBSTITUTIONS = 4
RNG_SEED = 0

# Non-VIII VR loops eligible for stacking with a 7-mer insertion
# (VR-VIII is the insertion site itself; mutating it would interfere)
STACKED_SUBSTITUTION_LOOPS = ["VR-IV", "VR-V", "VR-IX"]

# --- Closed-loop hyperparameters (used in later phases) ---
INFLAMMATION_THRESHOLD = 0.40
N_CYCLES = 10
K_PER_CYCLE = 5
N_PRETRAIN_CAMPAIGNS = 1000
SIMULATOR_NOISE_STD = 0.08
SIMULATOR_VERSION = "simulator_v1.0"

# --- ESM3 ---
ESM3_MODEL = "esmc-6b-2024-12"
ESM3_FORGE_URL = "https://biohub.ai"  # forge.evolutionaryscale.ai sunsets 2026-11-28
ESM3_EMBEDDING_DIM = 2560

# --- NCBI Entrez (compliance) ---
NCBI_EMAIL = "vramaswamy@gmail.com"
