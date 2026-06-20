# CLAUDE.md — AMD intravitreal AAV capsid optimization

## Project overview

Build a realistic implementation of a closed-loop world-model pipeline for
**AAV capsid optimization** targeting **dry AMD (geographic atrophy)** in the
retinal pigment epithelium. This is the *Experiment* leg of the JARVIS-for-bio
trinity, following the v0 *Explore* demo (`../v0_release/`).

The goal: demonstrate that a Gaussian-process world model + a small RL
policy (pretrained on a biologically grounded simulator, then deployed in
a 10-cycle closed loop) explores the Pareto frontier of
**RPE transduction × neutralizing-antibody escape** faster than random
selection — using real AAV2 VP1 sequences, real ESM3 embeddings, and a
wet-lab simulator calibrated to published intravitreal AAV literature.

### Threading from v0

v0 ran an AMD GWAS (Fritsche LG et al. 2016, GCST003219) through the
inference-oriented architecture and surfaced six top targets. **C9** (rank 3,
L2G 0.961, intronic regulatory) was assigned the hypothesis:

> An intronic variant cis-regulating C9 expression raises local complement
> C9 dosage in the RPE-choroid, biasing the terminal complement pathway
> toward more efficient C5b-9 (MAC) assembly on Bruch's membrane.

v1 acts on that hypothesis. The therapeutic angle: deliver **soluble CD59**
(sCD59) to RPE cells. sCD59 binds C8 and prevents C9 polymerization into the
MAC pore, neutralizing the dysregulated terminal cascade locally. Real-world
precedent: **HMR59 / JNJ-1887** (Hemera Biosciences → J&J), AAV2-sCD59
intravitreal, Phase 2 in dry and wet AMD.

The cassette is fixed. The capsid is the design variable.

### The load-bearing mock

The **wet-lab simulator is the only synthetic component**. It is calibrated
against published intravitreal AAV data (Dalkara 2013, Byrne 2020,
Kotterman 2015) but it is not a wet lab. Everything else — the AAV2
reference sequence, the ESM3 embeddings, the variant generation, the GP
world model, the RL policy, the SQLite store — is real. The README and
final report must label the simulator visibly. The honesty is load-bearing
for the architecture pitch.

The simulator does double duty: training ground for the policy during
pretraining, and the v1 evidence source during the live cycles. Both roles
use the same calibrated function. In v2+ the live evidence source swaps
to the actual wet lab; pretraining stays on the simulator.

---

## Biological context

### Disease

- Dry AMD with geographic atrophy: leading cause of irreversible blindness
  in adults >60 in developed economies
- ~5M US patients with intermediate or advanced AMD; ~1M with GA
- Approved drugs (Syfovre / pegcetacoplan, Izervay / avacincaptad pegol)
  hit complement but require **monthly intravitreal injections**
- One-time AAV gene therapy is the durability play

### Therapeutic objective

Deliver engineered AAV capsids carrying **secreted sCD59** to RPE cells via
**intravitreal injection**. The intravitreal route is what makes the capsid
optimization story matter: subretinal injection bypasses the vitreal
antibody compartment (immune-privileged), but intravitreal puts the capsid
directly in front of patient neutralizing antibodies. Capsid engineering
earns its keep on this route.

Therapeutic threshold (qualitative): durable RPE transduction sufficient to
hold sCD59 above the local MAC-neutralization Kd for >12 months, in patients
with anti-AAV2 NAb titers up to 1:160 (representative of the seroprevalent
adult population, Calcedo 2009).

### Target cell type

**Retinal pigment epithelium (RPE)** — monolayer between photoreceptors and
choroid. Site of drusen deposition and MAC accumulation in dry AMD. Also the
biosynthetic factory in the eye for secreted proteins, making it the right
cell for sCD59 expression.

Secondary targets: photoreceptors and Müller glia (engineered intravitreal
AAVs hit all three with different efficiencies — Dalkara 2013).

### Payload (fixed)

Cassette (NOT a design variable in this pipeline):

```
5'-ITR — CMV promoter — IgK signal peptide — mature CD59 (residues 26-102) — WPRE — bGH polyA — 3'-ITR
```

- Mature CD59 is ~77 residues / ~9 kDa (extracellular domain only,
  GPI-anchor stripped)
- IgK signal peptide replaces native GPI-anchor signal → secretion instead
  of membrane attachment
- Total cassette: ~2.3 kb, well within AAV's ~4.7 kb packaging limit
- Reference for the secreted-sCD59 design: Cashman et al. 2011, Mol Ther 19:1640

### Capsid (the design variable)

**Baseline: AAV2 VP1 wildtype.** NCBI accession `YP_680426.1` (or AAV2 capsid
PDB 1LP3 sequence). 735 amino acids. Human seroprevalence ~30–70%
(Calcedo 2009), which is exactly why capsid engineering matters for
intravitreal delivery.

**Reference engineered capsid: AAV.7m8** (Dalkara 2013). AAV2 backbone with
a 7-mer peptide insertion `LALGETTRP` between residues 587 and 588 (VR-VIII
loop), plus a small number of point mutations from the directed evolution
selection. 30–50× higher intravitreal RPE transduction than AAV2 in
murine + primate retina. Used as the anchor for the simulator's transduction
axis.

**Variable regions** (residue ranges in AAV2 VP1 numbering, from
Govindasamy 2006 / DiMattia 2012):

| Loop | Residues | Role |
|---|---|---|
| VR-I | 263–268 | Surface |
| VR-IV | 449–468 | **Major antibody epitope + tropism determinant** |
| VR-V | 488–505 | Surface |
| VR-VIII | 581–593 | **AAV.7m8 insertion site, tropism + epitope** |
| VR-IX | 704–714 | Surface, antibody epitope |

Variant generation focuses on VR-IV, VR-V, VR-VIII, VR-IX — the same loops
Tseng & Agbandje-McKenna 2014 mapped as antibody hotspots. Two move classes:

1. **Point substitutions** (1–6 amino acid swaps drawn from the loop's
   natural substitution profile across the AAV serotype family)
2. **7-mer peptide insertions** at position 587–588 (the AAV.7m8 trick),
   drawn from a curated library of retinal-tropic 7-mer sequences

Target: generate ~80 capsid variants for the initial candidate pool.

---

## Pipeline architecture

```
Layer 1 — Generative
  └── capsid_generator.py
        AAV2 VP1 wildtype as base
        Generate variants by:
          (a) substituting 1-6 residues at VR-IV, V, VIII, IX
          (b) inserting 7-mer peptides between residues 587-588
        Output: library of capsid variant sequences (FASTA)

        v2+ swap-in: replace this with a learned generator
        (RFdiffusion + ProteinMPNN over the VR loops, conditioned on the
        same Pareto reward signal). Layer 1's interface is just "yield
        VP1 sequences"; the rest of the pipeline doesn't change.

Layer 2 — Feature engineering
  └── esm3_embedder.py
        ESM3 Forge API (esmc-6b-2024-12)
        Input: capsid variant VP1 protein sequences (FASTA)
        Output: HDF5 directory, one file per sequence
        Embedding: mean-pooled across residues -> 2560-dim vector
        Also compute pseudo-likelihood score per variant
        Cache wildtype AAV2 and AAV.7m8 reference embeddings separately

Layer 3 — World model
  └── world_model.py
        Multi-output Gaussian process (GPyTorch + BoTorch)
        Input: ESM3 embedding (2560-dim)
          + 4 engineered features (peptide insertion present, insertion
            length, edit distance to AAV2, edit distance to AAV.7m8)
        Multi-output: predict
          (1) rpe_transduction   (0-1, intravitreal RPE transduction)
          (2) neut_escape        (0-1, fraction of NAb panel escaped)
          (3) inflammation_score (0-1, relative IFN-beta + cytokine response)
                                 -> CONSTRAINT, not Pareto axis
        Calibrated uncertainty estimates (posterior variance)
        Pre-train on public data (data/public/)

Layer 4 — RL policy (closed loop)
  ├── policy.py
  │     Small MLP (~10k params). Inputs per candidate:
  │       - Candidate features: ESM3 embedding (2560) + 4 engineered features
  │       - World-model output: (mean, var) for transduction, escape, inflammation
  │     Output: scalar selection_score
  │     Selection: top-k by score, filtered by predicted
  │       inflammation_score < INFLAMMATION_THRESHOLD (default 0.40)
  │     Trained by REINFORCE on Pareto hypervolume improvement reward
  │     Pretrained on ~1000 simulated 10-cycle campaigns before the live run
  │
  ├── wet_lab_simulator.py  *** THE LOAD-BEARING MOCK ***
  │     Evidence source. v1 stand-in for the real wet lab.
  │     Same interface as a real wet lab; different implementation.
  │     Calibrated to: Dalkara 2013 (AAV.7m8), Byrne 2020 (4D-R100),
  │     Kotterman 2015 (vitreal NAb), Reichel 2017 / Bucher 2021 (inflammation),
  │     HMR59 / ADVM-022 trial readouts (where public).
  │     Returns (rpe_transduction, neut_escape, inflammation_score) per variant.
  │     Add Gaussian noise (sigma=0.08) to simulate biological + technical variation.
  │
  └── closed_loop.py
        Run N=10 cycles, select k=5 candidates per cycle.
        Each cycle: GP world model refits on accumulated data; RL policy
        scores untested candidates; top-k go to evidence source (simulator
        in v1, wet lab in v2+).
        Log Pareto frontier evolution + constraint violations.
        Save all results to SQLite.
        Compare against random-baseline selection (run in parallel).
```

Note what's gone from the original DMD CLAUDE.md: no CRISPOR, no
ViennaRNA, no guide_features.py, no encode_features.py. No CRISPR layer
at all. The cassette is fixed; the capsid is the design variable; that's
the entire story. The Pareto axes are cleaner as a result.

---

## Data storage conventions

```
data/
├── sequences/
│   ├── aav2_vp1_wildtype.fasta         # NCBI YP_680426.1
│   ├── aav7m8_vp1_reference.fasta      # AAV2 + LALGETTRP @ 587-588
│   └── capsid_variants.fasta            # generated variants (~80)
├── embeddings/
│   └── esm3/
│       ├── {sequence_hash}.h5           # shape: (2560,) mean-pooled embedding
│       └── index.parquet                # sequence_id -> file path + pseudo-LL score
├── features/
│   └── capsid_features.parquet          # engineered features per variant
├── public/                              # pre-training data from literature
│   ├── dalkara_2013.csv                # AAV.7m8 + directed-evolution retinal hits
│   ├── byrne_2020.csv                  # engineered intravitreal capsid library
│   └── kotterman_2015.csv              # vitreal NAb seroprevalence + escape
├── pretrained_policy.pt                 # RL policy weights from scripts/pretrain_policy.py
└── results.db                           # SQLite: all experimental results
```

---

## ESM3 API usage

```python
from esm.sdk.forge import ESM3ForgeInferenceClient
from esm.sdk.api import ESMProtein, LogitsConfig

forge_client = ESM3ForgeInferenceClient(
    model="esmc-6b-2024-12",
    url="https://forge.evolutionaryscale.ai",
    token=os.environ["ESM3_API_TOKEN"],
)

def embed_sequence(sequence: str) -> np.ndarray:
    """Embed AAV2 VP1 capsid variant. Returns 2560-dim mean-pooled vector."""
    protein = ESMProtein(sequence=sequence)
    protein_tensor = forge_client.encode(protein)
    logits_output = forge_client.logits(
        protein_tensor,
        LogitsConfig(sequence=True, return_embeddings=True),
    )
    embedding = logits_output.embeddings.mean(dim=1).squeeze().numpy()
    return embedding
```

Cache to HDF5 immediately. Check cache before any API call. Compute and
cache embeddings for AAV2 wildtype and AAV.7m8 reference *first* — every
simulator call depends on distances to these anchors.

---

## Capsid variant generation

```python
def generate_variants(
    wildtype_seq: str,
    n_substitution_variants: int = 50,
    n_insertion_variants: int = 30,
    rng: np.random.Generator = np.random.default_rng(0),
) -> list[dict]:
    """
    Generate AAV2 VP1 capsid variants.

    Substitution variants: 1-6 residue swaps drawn from VR-IV, V, VIII, IX.
    Substitution choices weighted by natural amino-acid frequency in
    homologous AAV serotype VR loops (AAV1, 3, 5, 6, 8, 9, rh10, rh74).

    Insertion variants: 7-mer peptide inserted between residues 587 and 588.
    7-mer library seeded with AAV.7m8 (LALGETTRP) + 29 engineered variants
    drawn from the Dalkara 2013 / Byrne 2020 retinal-tropic peptide pool.
    """
    pass
```

Resulting variants get assigned a stable `capsid_id` (8-char hash of
sequence) and written to `data/sequences/capsid_variants.fasta`.

---

## World model implementation

```python
import gpytorch, botorch, torch
from botorch.models import MultiTaskGP
from botorch.fit import fit_gpytorch_mll

class AMDCapsidWorldModel:
    """
    Multi-output Gaussian process world model predicting:
      1. rpe_transduction       (0-1, intravitreal RPE transduction)  [Pareto axis]
      2. neut_escape            (0-1, fraction of serum panel escaped) [Pareto axis]
      3. inflammation_score     (0-1, relative innate immune activation) [CONSTRAINT]

    Input features: ESM3 embedding (2560-dim) + 4 engineered features:
      - has_7mer_insertion      (binary)
      - insertion_length        (int, 0 if none)
      - hamming_to_aav2         (int, normalized)
      - cosine_to_aav7m8_embed  (float, ESM3-space)

    Total input dim: 2564
    """

    def __init__(self):
        self.model = None
        self.train_X = None
        self.train_Y = None
        self.feature_dim = 2564

    def fit(self, X: torch.Tensor, Y: torch.Tensor):
        """Fit MultiTaskGP on accumulated experimental data."""
        # normalize inputs; standardize Y
        # MultiTaskGP for correlated outputs
        # fit_gpytorch_mll
        pass

    def predict(self, X: torch.Tensor) -> tuple[torch.Tensor, torch.Tensor]:
        """Return (mean, variance) for both outputs over candidate set."""
        pass
```

---

## Wet lab simulator (THE MOCK)

Biologically grounded, not random. Calibrated against published intravitreal
AAV literature. **Label visibly as a simulator in every artifact downstream.**

```python
def simulate_wet_lab(
    capsid_embedding: np.ndarray,
    aav2_embedding: np.ndarray,
    aav7m8_embedding: np.ndarray,
    has_7mer_insertion: bool,
    insertion_length: int,
    hamming_to_aav2: int,
    noise_std: float = 0.08,
    rng: np.random.Generator = np.random.default_rng(42),
) -> dict:
    """
    Simulate intravitreal AAV outcomes from a single capsid variant.

    rpe_transduction:
      - Anchored on AAV.7m8 (~0.45 mean intravitreal RPE transduction, ~30-50x
        baseline AAV2 per Dalkara 2013)
      - Baseline AAV2: ~0.015 (intravitreal AAV2 reaches RPE poorly)
      - Boosted by 7-mer insertion at 587-588 (the AAV.7m8 trick)
      - Boosted by cosine similarity to AAV.7m8 in ESM3 space
      - Degraded if hamming_to_aav2 is too large (capsid folding / packaging fails)

    neut_escape:
      - Anchored on AAV2 NAb prevalence (~50% population GMT >= 1:20, Calcedo 2009)
      - AAV2 escape baseline: ~0.10 of population escapes wildtype-targeting NAbs
      - AAV.7m8 partial escape (~0.30) per Kotterman 2015 vitreal NAb data
      - Escape rises with distance from AAV2 in ESM3 space (epitope erosion)
      - Falls off sharply when capsid structure breaks (hamming too high)

    Pareto tension: VR-VIII drives BOTH axes. AAV.7m8-like insertions raise
    transduction AND escape, but past some divergence threshold capsid fitness
    collapses (packaging yield, virion stability). The GP has to learn this.
    """

    cos_to_aav7m8 = float(
        np.dot(capsid_embedding, aav7m8_embedding)
        / (np.linalg.norm(capsid_embedding) * np.linalg.norm(aav7m8_embedding))
    )
    cos_to_aav2 = float(
        np.dot(capsid_embedding, aav2_embedding)
        / (np.linalg.norm(capsid_embedding) * np.linalg.norm(aav2_embedding))
    )
    embed_distance_to_aav2 = 1.0 - cos_to_aav2  # 0 = identical, larger = farther

    # Fitness cliff: capsid folding/packaging fails when too divergent
    fitness_penalty = np.tanh(max(0, hamming_to_aav2 - 10) / 6.0)  # 0-1

    # --- Axis 1: RPE transduction (intravitreal) ---
    insertion_bonus = 0.35 if has_7mer_insertion else 0.0
    similarity_bonus = 0.40 * max(0, cos_to_aav7m8)
    baseline_aav2 = 0.015
    rpe_transduction = (
        baseline_aav2 + insertion_bonus + similarity_bonus
    ) * (1.0 - 0.6 * fitness_penalty)

    # --- Axis 2: NAb escape ---
    escape_baseline = 0.10
    escape_from_divergence = 0.65 * np.tanh(embed_distance_to_aav2 / 0.10)
    neut_escape = (escape_baseline + escape_from_divergence) * (
        1.0 - 0.4 * fitness_penalty
    )

    # --- Constraint: inflammation_score ---
    # AAV2 baseline ~0.15 (mild but present). AAV.7m8 ~0.35 (the dose-limiting
    # signal that triggered ADVM-022 dose reduction in clinic). Rises with
    # capsid divergence (more non-self -> more innate sensing); rises sharply
    # past fitness cliff (misfolded capsid -> aggregation -> TLR2).
    inflammation_baseline = 0.15
    inflammation_from_divergence = 0.30 * np.tanh(embed_distance_to_aav2 / 0.12)
    inflammation_from_cliff = 0.50 * fitness_penalty
    inflammation_score = (
        inflammation_baseline + inflammation_from_divergence + inflammation_from_cliff
    )

    # noise
    rpe_transduction += rng.normal(0, noise_std)
    neut_escape += rng.normal(0, noise_std)
    inflammation_score += rng.normal(0, noise_std)

    return {
        "rpe_transduction": float(np.clip(rpe_transduction, 0, 1)),
        "neut_escape": float(np.clip(neut_escape, 0, 1)),
        "inflammation_score": float(np.clip(inflammation_score, 0, 1)),
        "_simulator_version": "v1.0",  # bump if calibration changes
        "_is_simulated": True,
    }
```

The Pareto tension is real because both axes are pulled by VR-VIII
divergence, but the fitness cliff penalizes runaway divergence. The GP
learns the curved frontier.

---

## Public pre-training data

Synthetic-but-grounded CSVs in `data/public/`, calibrated to published
intravitreal AAV capsid datasets. These seed the GP before the simulated
campaign begins.

**dalkara_2013.csv** — AAV.7m8 directed-evolution hits
- Columns: `capsid_id, vp1_sequence, rpe_transduction, neut_escape, inflammation_score, source`
- ~15 rows: AAV2 wildtype, AAV.7m8, plus partial-selection intermediates
- Source: Dalkara D et al. 2013, Sci Transl Med 5:189ra76

**byrne_2020.csv** — engineered intravitreal library readouts
- Columns: `capsid_id, vp1_sequence, rpe_transduction, neut_escape, inflammation_score, source`
- ~12 rows: 4D-R100 family + control capsids
- Source: Byrne LC et al. 2020, J Clin Invest 130:4214

**kotterman_2015.csv** — vitreal NAb seroprevalence
- Columns: `capsid_id, vp1_sequence, neut_escape, n_serum_samples, source`
- ~8 rows: NAb prevalence vs panel of natural + engineered capsids
- Source: Kotterman MA et al. 2015, Gene Ther 22:116

**reichel_bucher_inflammation.csv** — intravitreal AAV ocular inflammation
- Columns: `capsid_id, vp1_sequence, inflammation_score, ifn_beta_fold, source`
- ~10 rows: AAV2, AAV.7m8 (ADVM-022 backbone, clinical inflammation
  data), 4D-R100 family, plus murine-only control capsids
- Source: Reichel FF et al. 2017, Mol Ther 25:2648 (intravitreal AAV
  inflammation in NHP); Bucher K et al. 2021, Mol Ther 29:1985
  (CpG-driven innate immune sensing of AAV genomes)

Total pre-training rows: ~45.

---

## RL policy + closed loop

Layer 4 has three pieces: a policy network, an evidence source, and the loop
that wires them together. The policy is small. The loop is short. The
"intelligence" lives in pretraining against the simulator.

### Policy network

Small MLP (~10k params). Inputs per candidate:
- Candidate features: ESM3 embedding (2560) + 4 engineered features
- World-model output: `(mean_t, var_t, mean_n, var_n, mean_i, var_i)` for
  rpe_transduction / neut_escape / inflammation_score

Output: scalar `selection_score`.

Selection: top-k by score, filtered by predicted
`inflammation_score < INFLAMMATION_THRESHOLD`.

Training: REINFORCE. Reward = Pareto hypervolume improvement after the
cycle's k candidates are observed. Policy weights updated each cycle.

### Pretraining (sim2real)

Before the figure-generating live run, the policy is pretrained against the
simulator:

```python
for campaign in range(1000):
    candidates = generate_random_variant_pool(n=80)
    observed = pretraining_anchors + public_data
    for cycle in range(10):
        world_model = fit_gp(observed)
        scores = policy(candidates, world_model.predict(candidates))
        picks = top_k_constrained(scores, k=5)
        outcomes = simulator(picks)
        reward = pareto_hypervolume_improvement(observed, outcomes)
        policy.reinforce(reward)
        observed = observed + outcomes
    # next campaign: regenerate candidate pool, reset observed
```

~30 seconds wall time on CPU. Outputs `data/pretrained_policy.pt`.

The 80 candidates of the figure-generating run are sampled fresh; the
policy never sees them during pretraining. The policy learns the *shape*
of the search space, not specific answers.

### Live cycle (the figure-generating run)

Same loop, k=5 per cycle, 10 cycles. Pretrained policy is the starting
point. Each cycle:

1. GP world model refits on accumulated data (public + observed in this run)
2. Policy scores all untested candidates
3. Top-k by score (constraint-filtered) go to the evidence source
4. Outcomes append to training set; policy REINFORCE step
5. Pareto frontier + hypervolume logged to SQLite

### Random baseline

Same budget, same constraint filter, random selection. Plotted alongside.
Headline plot: `outputs/rl_vs_random.png` — hypervolume convergence over
cycles.

---

## Closed-loop architecture: evidence sources

The closed loop is defined by what feeds it. v1 uses the wet-lab simulator
as its evidence source. v2+ swaps in the actual wet lab. The rest of the
architecture is unchanged.

Same interface:

```python
class EvidenceSource(Protocol):
    def measure(self, candidates: list[CapsidVariant]) -> list[Outcome]:
        """(rpe_transduction, neut_escape, inflammation_score) per variant."""
        ...
```

v1 implementations:
- `WetLabSimulator` — calibrated model in `wet_lab_simulator.py`. Cheap,
  fast, the load-bearing mock.

v2+ implementations (out of scope here, but the interface is the contract):
- `WetLabReal` — Opentrons-driven AAV production + HTRF/qPCR readouts.
  Slow, expensive, ground truth.
- `HybridSource` — sim for pretraining, real for live cycles. The sim2real
  setup we recommend in production.

The policy doesn't know or care which source it's reading from. The
`source` column in `experimental_results` carries the provenance so every
data point is traceable.

---

## SQLite schema

```sql
CREATE TABLE candidates (
    candidate_id TEXT PRIMARY KEY,
    capsid_id TEXT NOT NULL,
    vp1_sequence TEXT NOT NULL,
    has_7mer_insertion BOOLEAN,
    insertion_peptide TEXT,
    hamming_to_aav2 INTEGER,
    capsid_embedding_path TEXT,
    is_anchor BOOLEAN DEFAULT FALSE,    -- AAV2 wildtype + AAV.7m8 reference
    anchor_label TEXT,                  -- 'AAV2', 'AAV.7m8', NULL otherwise
    created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
);

CREATE TABLE experimental_results (
    result_id INTEGER PRIMARY KEY AUTOINCREMENT,
    candidate_id TEXT REFERENCES candidates(candidate_id),
    cycle INTEGER NOT NULL,             -- 0 for pretraining + anchors
    selection_strategy TEXT,            -- 'rl_policy' | 'random_baseline' | 'pretraining' | 'anchor'
    rpe_transduction REAL,
    neut_escape REAL,
    inflammation_score REAL,
    meets_constraint BOOLEAN,           -- inflammation_score < threshold
    is_on_pareto_frontier BOOLEAN,      -- recomputed at end of run
    source TEXT NOT NULL,               -- 'public_literature' | 'wet_lab_simulator' | 'wet_lab_real'
    source_version TEXT,                -- 'simulator_v1.0' | 'opentrons_2026-Q2' | DOI of public dataset
    is_simulated BOOLEAN DEFAULT TRUE,  -- redundant honesty flag (derived from source)
    assay_date TIMESTAMP DEFAULT CURRENT_TIMESTAMP
);

CREATE TABLE cycle_summary (
    cycle INTEGER,
    selection_strategy TEXT,
    n_candidates_tested INTEGER,
    n_constraint_violations INTEGER,    -- inflammation_score >= threshold
    best_rpe_transduction REAL,
    best_neut_escape REAL,
    pareto_frontier_size INTEGER,
    pareto_hypervolume REAL,
    timestamp TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
    PRIMARY KEY (cycle, selection_strategy)
);
```

The figure-generating export is a single parquet built by joining these tables:

```sql
-- outputs/pareto_data.parquet
SELECT
    c.candidate_id, c.capsid_id, c.vp1_sequence,
    c.has_7mer_insertion, c.insertion_peptide, c.hamming_to_aav2,
    c.is_anchor, c.anchor_label,
    r.cycle, r.selection_strategy,
    r.rpe_transduction, r.neut_escape, r.inflammation_score,
    r.meets_constraint, r.is_on_pareto_frontier,
    r.source, r.source_version, r.is_simulated
FROM candidates c
JOIN experimental_results r ON c.candidate_id = r.candidate_id;
```

Every column the Pareto figure needs comes from this single file.

---

## Project structure

```
v1_release/
├── CLAUDE.md
├── README.md
├── requirements.txt
├── .env.example                         # ESM3_API_TOKEN=your_token_here
├── config.py                            # constants, paths, hyperparameters
├── data/
│   ├── sequences/
│   ├── embeddings/esm3/
│   ├── features/
│   ├── public/
│   └── results.db
├── pipeline/
│   ├── __init__.py
│   ├── layer1_generative/
│   │   └── capsid_generator.py          # AAV2 VR-loop variant + 7-mer insertion
│   ├── layer2_features/
│   │   └── esm3_embedder.py             # ESM3 Forge API + HDF5 cache
│   ├── layer3_world_model/
│   │   └── world_model.py               # MultiTaskGP
│   └── layer4_closed_loop/
│       ├── policy.py                    # small MLP, REINFORCE
│       ├── wet_lab_simulator.py         # *** THE MOCK (evidence source v1) ***
│       └── closed_loop.py               # cycle driver
├── visualization/
│   ├── pareto.py                        # frontier evolution plot
│   └── rl_vs_random.py                  # headline comparison plot
├── scripts/
│   ├── fetch_sequences.py               # AAV2 VP1 from NCBI YP_680426.1
│   ├── build_aav7m8_reference.py        # AAV2 + LALGETTRP @ 587-588
│   ├── seed_public_data.py              # build pre-training CSVs
│   └── pretrain_policy.py               # 1000 simulated campaigns -> pretrained_policy.pt
└── main.py                              # end-to-end pipeline runner
```

---

## Requirements

```
# Core
numpy>=1.24
pandas>=2.0
scipy>=1.11
torch>=2.0
h5py>=3.9
pyarrow>=13.0
sqlalchemy>=2.0

# ESM3
esm>=3.0

# ML
gpytorch>=1.11
botorch>=0.9

# Bio
biopython>=1.81
requests>=2.31

# Viz
matplotlib>=3.7
seaborn>=0.13

# Utilities
python-dotenv>=1.0
tqdm>=4.66
loguru>=0.7
```

No ViennaRNA. No CRISPOR. The DMD-era guide-design dependencies are gone.

---

## Setup instructions (for README.md)

```bash
# 1. Dependencies (shared venv per project convention)
~/venv/bin/pip install -r requirements.txt

# 2. ESM3 API token
cp .env.example .env
# Edit .env and add your ESM3 Forge token

# 3. Fetch real sequences
~/venv/bin/python scripts/fetch_sequences.py        # AAV2 VP1 from NCBI
~/venv/bin/python scripts/build_aav7m8_reference.py # AAV.7m8 reference

# 4. Seed public pre-training data
~/venv/bin/python scripts/seed_public_data.py

# 5. Pretrain the RL policy on the simulator (~30 seconds)
~/venv/bin/python scripts/pretrain_policy.py    # -> data/pretrained_policy.pt

# 6. Run pipeline
~/venv/bin/python main.py
```

---

## Expected output

The pipeline's primary deliverable is **the variants + the data table that draws
the Pareto figure** (`pareto_frontier.png` design reference in repo root). Every
column the figure needs is in `pareto_data.parquet`; the renderer is a thin
matplotlib script over it.

```
=== AMD AAV Capsid Optimization Pipeline ===

[Layer 1] Generating capsid variants...
  AAV2 wildtype loaded (735 aa, YP_680426.1)
  AAV.7m8 reference built (LALGETTRP inserted at 587-588)
  Substitution variants: 50
  7-mer insertion variants: 30
  Total candidate pool: 80

[Layer 2] Computing ESM3 embeddings...
  Forge API: esmc-6b-2024-12
  82 sequences embedded (80 variants + AAV2 + AAV.7m8 anchors)
  Cache: data/embeddings/esm3/  (HDF5, 2560-dim mean-pooled)

[Layer 3] Pre-training world model on public data (n=45)...
  MultiTaskGP fitted:
    rpe_transduction    R^2 = 0.58
    neut_escape         R^2 = 0.61
    inflammation_score  R^2 = 0.54  (constraint)
  Initial Pareto frontier: 4 candidates (constraint-satisfying)

[Layer 4] Loading pretrained RL policy (data/pretrained_policy.pt)...
  Policy: 2-layer MLP, 10,272 parameters
  Trained on 1,000 simulated 10-cycle campaigns (50,000 simulator calls)
  Pretraining wall time: 28 seconds

[Layer 4] Live closed loop (10 cycles, k=5/cycle, RL policy vs random)...
  Cycle  1: rl     pareto_hv=0.043  random pareto_hv=0.038  (0 violations)
  Cycle  2: rl     pareto_hv=0.061  random pareto_hv=0.041  (1 violation, filtered)
  Cycle  3: rl     pareto_hv=0.083  random pareto_hv=0.048
  ...
  Cycle 10: rl     pareto_hv=0.151  random pareto_hv=0.072

  Final Pareto frontier (RL policy, constraint-satisfying):
    rpe_transduction range: 0.31 - 0.52
    neut_escape range:      0.40 - 0.78
    pareto_size:            14 capsids
```

### Deliverables (figure-generating dataset)

**`outputs/variants.fasta`** — every capsid variant the pipeline considered,
including anchors. FASTA headers carry `capsid_id` + metadata. This is the
"someone could synthesize and order these" file.

**`outputs/pareto_data.parquet`** — one row per (candidate, cycle) observation.
This is the table that draws the figure. Columns:

| Column | Type | Purpose in figure |
|---|---|---|
| `capsid_id` | str | identity |
| `vp1_sequence` | str | re-derivation |
| `has_7mer_insertion` | bool | variant class |
| `insertion_peptide` | str (nullable) | variant class |
| `hamming_to_aav2` | int | divergence |
| `is_anchor` | bool | special markers (AAV2 ★, AAV.7m8 ★) |
| `anchor_label` | str (nullable) | legend |
| `cycle` | int | when it was tested |
| `selection_strategy` | str | grey dots vs orange RL picks vs pretraining |
| `rpe_transduction` | float | x-axis |
| `neut_escape` | float | y-axis |
| `inflammation_score` | float | constraint shading |
| `meets_constraint` | bool | inside vs outside safe zone |
| `is_on_pareto_frontier` | bool | dashed-line points |
| `source` | str | evidence-source provenance (public / sim / wet lab) |
| `source_version` | str | honesty flag |
| `is_simulated` | bool | honesty flag |

**`outputs/pareto_frontier.png`** — rendered from `pareto_data.parquet` by
`visualization/pareto.py`. Layout matches the conceptual reference figure at
repo root: scatter of all candidates, RL selections highlighted, AAV2 and
AAV.7m8 anchors marked, threshold lines for min transduction + min escape,
inflammation-constraint zone shaded red, target zone boxed.

**`outputs/closed_loop_summary.csv`** — cycle_summary table dump (one row
per cycle × strategy). Used for the `rl_vs_random.png` companion plot
showing hypervolume convergence.

**`data/results.db`** — SQLite, the source of truth. The parquet + csv are
exports.

---

## Visualization spec

`visualization/pareto.py` reads `outputs/pareto_data.parquet` and produces
`outputs/pareto_frontier.png`. The figure reproduces the conceptual reference
at the repo root (`pareto_frontier.png` — design target).

Layout:

```
X axis: RPE Transduction Efficiency (intravitreal route)  0 -> 1
Y axis: Neutralizing Antibody Escape (fraction of panel)  0 -> 1

Scatter layers (back to front):
  1. Grey dots       — all candidates ever evaluated (pretraining + live)
  2. Orange circles  — RL-policy selections (one per cycle, colored by cycle)
  3. Red diamond     — AAV2 wildtype anchor                  (labeled)
  4. Green diamond   — AAV.7m8 anchor (Dalkara 2013)         (labeled)
  5. Purple star     — target zone centroid                  (labeled)
  6. Blue dashed     — Pareto frontier line (non-dominated, constraint-passing)

Reference lines:
  - Vertical dashed: min RPE transduction threshold (~0.10)
  - Horizontal dashed: min NAb escape threshold (~0.50)

Shaded regions:
  - Green: safe therapeutic window (above both thresholds, below tox cliff)
  - Red:   inflammation-constraint violation zone
           (capsids with inflammation_score >= INFLAMMATION_THRESHOLD,
            plotted with reduced opacity + red background hatch)

Title:     "World model is exploring a Pareto frontier"
Subtitle:  "RL policy navigates two engineering objectives
            under an inflammation-safety constraint"

Honesty footer (small, bottom-right): "Simulated data — wet-lab outcomes
generated by biologically grounded simulator (v1.0). See README."
```

The script is intentionally thin: ~150 lines of matplotlib over the
parquet. No business logic, no recomputation. If a column in
`pareto_data.parquet` is wrong, the figure is wrong — there's nowhere else
to hide.

A companion script `visualization/rl_vs_random.py` renders the cycle-vs-
hypervolume convergence curve from `closed_loop_summary.csv`. This is the
second headline plot: shows the RL policy beating random over the 10-cycle
budget.

---

## Success criteria

1. ESM3 embeddings are computed from real AAV2 capsid variant sequences and
   cached to HDF5
2. AAV2 wildtype and AAV.7m8 references are real sequences from NCBI / the
   Dalkara 2013 published construct
3. GP world model R² improves across cycles as data accumulates (all three
   outputs: rpe_transduction, neut_escape, inflammation_score)
4. `scripts/pretrain_policy.py` produces `data/pretrained_policy.pt` in
   under 60 seconds of CPU compute; main.py loads it before the live cycles
5. RL policy outperforms random baseline on Pareto hypervolume (visible
   separation by cycle 4–5) WITHOUT exceeding the inflammation constraint
   threshold on selected candidates
6. Pareto frontier advances toward `rpe_transduction > 0.35 AND
   neut_escape > 0.60 AND inflammation_score < 0.40` (the AAV.7m8 +
   intravitreal-engineered envelope the field has published)
7. `outputs/pareto_data.parquet` contains every column the figure needs;
   `outputs/variants.fasta` lists every capsid that was generated
8. `visualization/pareto.py` regenerates `outputs/pareto_frontier.png`
   matching the conceptual reference figure (anchors marked, constraint
   zone shaded, target zone visible)
9. Every row in `experimental_results` carries a `source` value; the
   simulator-sourced rows additionally carry `is_simulated=TRUE` and a
   `source_version` of `simulator_v1.0`. Figures footer the simulator
   visibly.
10. `python main.py` runs end to end

---

## Key references

- **Dalkara D et al. 2013, Sci Transl Med 5:189ra76** — AAV.7m8 directed
  evolution, intravitreal RPE transduction
- **Byrne LC et al. 2020, J Clin Invest 130:4214** — engineered intravitreal
  AAV capsids for retina
- **Kotterman MA et al. 2015, Gene Ther 22:116** — vitreal NAb escape in
  primate
- **Calcedo R et al. 2009, J Infect Dis 199:381** — global AAV2
  seroprevalence
- **Tseng YS & Agbandje-McKenna M 2014, Front Immunol 5:9** — AAV antibody
  epitope mapping
- **Govindasamy L et al. 2006, J Virol 80:11556** — AAV2 capsid structure +
  VR loop assignments
- **Cashman SM et al. 2011, Mol Ther 19:1640** — secreted-sCD59 AAV cassette
  design
- **Heier JS et al. 2020 (ARVO)** — HMR59 AAV2-sCD59 Phase 1 dry AMD
- **Fritsche LG et al. 2016, Nat Genet 48:134** — AMD GWAS (GCST003219, the
  v0 input study)
