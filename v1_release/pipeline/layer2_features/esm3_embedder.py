"""ESM3 Forge embeddings for AAV2 VP1 capsid variants.

Per-sequence cache: one HDF5 file per sequence (named by 8-char SHA-1).
An index parquet maps capsid_id -> embedding path + pseudo-LL score.

Run as a script to embed AAV2 + AAV.7m8 anchors + all 80 candidate variants.
"""

import hashlib
import os
import sys
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
import torch
from Bio import SeqIO
from dotenv import load_dotenv
from esm.sdk.api import ESMProtein, LogitsConfig
from esm.sdk.forge import ESM3ForgeInferenceClient
from tqdm import tqdm

sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent))
import config


def sequence_hash(seq: str) -> str:
    return hashlib.sha1(seq.encode()).hexdigest()[:8]


def make_client() -> ESM3ForgeInferenceClient:
    load_dotenv(config.PROJECT_ROOT / ".env")
    token = os.environ.get("ESM3_API_KEY")
    if not token:
        raise RuntimeError(
            "ESM3_API_KEY not set. Copy from v0_release/.env or set in v1_release/.env."
        )
    return ESM3ForgeInferenceClient(
        model=config.ESM3_MODEL,
        url=config.ESM3_FORGE_URL,
        token=token,
    )


def embed_one(
    client: ESM3ForgeInferenceClient, sequence: str
) -> tuple[np.ndarray, float]:
    """Return (mean-pooled embedding, pseudo log-likelihood per residue)."""
    protein = ESMProtein(sequence=sequence)
    protein_tensor = client.encode(protein)
    out = client.logits(
        protein_tensor,
        LogitsConfig(sequence=True, return_embeddings=True),
    )

    # Mean-pooled embedding across residues
    embedding = out.embeddings.mean(dim=1).squeeze().float().cpu().numpy()

    # Pseudo log-likelihood: mean log-prob per residue.
    # logits.sequence shape: [L+2, vocab]; protein_tensor.sequence shape: [L+2]
    # Positions 0 and -1 are special tokens (BOS/EOS), skip them.
    logits = out.logits.sequence.float()
    token_ids = protein_tensor.sequence
    log_probs = torch.log_softmax(logits, dim=-1)
    residue_positions = torch.arange(1, logits.shape[0] - 1)
    residue_tokens = token_ids[1:-1]
    per_residue_ll = log_probs[residue_positions, residue_tokens]
    pseudo_ll = float(per_residue_ll.mean().item())

    return embedding, pseudo_ll


def cache_path(hash_str: str) -> Path:
    return config.EMBEDDINGS_DIR / f"{hash_str}.h5"


def save_embedding(
    path: Path, embedding: np.ndarray, sequence: str, pseudo_ll: float
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with h5py.File(path, "w") as f:
        f.create_dataset("embedding", data=embedding)
        f.attrs["sequence"] = sequence
        f.attrs["sequence_hash"] = sequence_hash(sequence)
        f.attrs["pseudo_ll"] = pseudo_ll
        f.attrs["model"] = config.ESM3_MODEL
        f.attrs["embedding_dim"] = int(embedding.shape[0])


def load_embedding(path: Path) -> tuple[np.ndarray, float]:
    with h5py.File(path, "r") as f:
        return np.array(f["embedding"]), float(f.attrs["pseudo_ll"])


def parse_anchor_fasta(path: Path, anchor_label: str) -> dict:
    rec = SeqIO.read(path, "fasta")
    return {
        "capsid_id": anchor_label,
        "sequence": str(rec.seq),
        "is_anchor": True,
        "anchor_label": anchor_label,
    }


def parse_variants_fasta(path: Path) -> list[dict]:
    records = []
    for rec in SeqIO.parse(path, "fasta"):
        records.append(
            {
                "capsid_id": rec.id.split("|")[0],
                "sequence": str(rec.seq),
                "is_anchor": False,
                "anchor_label": None,
            }
        )
    return records


def parse_public_fasta(path: Path) -> list[dict]:
    """Parse the public-data FASTA. Skips records whose capsid_id collides
    with the anchor labels we already loaded; anchors get a single canonical
    entry from their dedicated FASTAs."""
    records = []
    for rec in SeqIO.parse(path, "fasta"):
        cid = rec.id.split("|")[0]
        if cid in ("AAV2", "AAV7m8"):
            continue
        records.append(
            {
                "capsid_id": cid,
                "sequence": str(rec.seq),
                "is_anchor": False,
                "anchor_label": None,
            }
        )
    return records


def embed_all(
    client: ESM3ForgeInferenceClient, records: list[dict]
) -> list[dict]:
    config.EMBEDDINGS_DIR.mkdir(parents=True, exist_ok=True)
    results = []
    for rec in tqdm(records, desc="Embedding"):
        seq = rec["sequence"]
        h = sequence_hash(seq)
        p = cache_path(h)
        if p.exists():
            _, pseudo_ll = load_embedding(p)
            cached = True
        else:
            emb, pseudo_ll = embed_one(client, seq)
            save_embedding(p, emb, seq, pseudo_ll)
            cached = False
        results.append(
            {
                **rec,
                "sequence_hash": h,
                "embedding_path": str(p.relative_to(config.PROJECT_ROOT)),
                "pseudo_ll": pseudo_ll,
                "sequence_length": len(seq),
                "cached": cached,
            }
        )
    return results


def build_index(records: list[dict]) -> pd.DataFrame:
    return pd.DataFrame(
        [
            {
                "capsid_id": r["capsid_id"],
                "sequence_hash": r["sequence_hash"],
                "embedding_path": r["embedding_path"],
                "pseudo_ll": r["pseudo_ll"],
                "is_anchor": r["is_anchor"],
                "anchor_label": r["anchor_label"],
                "sequence_length": r["sequence_length"],
            }
            for r in records
        ]
    )


def main() -> None:
    records = [
        parse_anchor_fasta(config.AAV2_FASTA, "AAV2"),
        parse_anchor_fasta(config.AAV7M8_FASTA, "AAV7m8"),
    ]
    records.extend(parse_variants_fasta(config.CAPSID_VARIANTS_FASTA))

    public_fasta = config.PUBLIC_DATA_DIR / "public_sequences.fasta"
    if public_fasta.exists():
        records.extend(parse_public_fasta(public_fasta))

    print(f"Embedding {len(records)} sequences with {config.ESM3_MODEL}")
    print(f"  endpoint: {config.ESM3_FORGE_URL}")
    print(f"  cache:    {config.EMBEDDINGS_DIR}")

    client = make_client()
    results = embed_all(client, records)

    df = build_index(results)
    df.to_parquet(config.EMBEDDINGS_INDEX, index=False)

    n_new = sum(not r["cached"] for r in results)
    n_cached = sum(r["cached"] for r in results)
    print("\nDone.")
    print(f"  newly embedded:   {n_new}")
    print(f"  loaded from cache: {n_cached}")
    print(f"  index parquet:    {config.EMBEDDINGS_INDEX}")

    aav2 = df[df.anchor_label == "AAV2"].iloc[0]
    aav7m8 = df[df.anchor_label == "AAV7m8"].iloc[0]
    print("\nAnchor embeddings:")
    print(
        f"  AAV2:    hash={aav2.sequence_hash}  len={aav2.sequence_length} "
        f"pseudo_ll={aav2.pseudo_ll:.3f}"
    )
    print(
        f"  AAV7m8:  hash={aav7m8.sequence_hash}  len={aav7m8.sequence_length} "
        f"pseudo_ll={aav7m8.pseudo_ll:.3f}"
    )

    var = df[~df.is_anchor]
    print(
        f"\nVariant pseudo_ll: min={var.pseudo_ll.min():.3f} "
        f"max={var.pseudo_ll.max():.3f} mean={var.pseudo_ll.mean():.3f}"
    )


if __name__ == "__main__":
    main()
