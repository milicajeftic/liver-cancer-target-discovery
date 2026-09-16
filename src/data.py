"""Download and load TCGA-LIHC data from the UCSC Xena GDC hub.

Usage:
    python -m src.data            # download all files into data/raw/

If a download fails (URLs on the hub occasionally change), open
https://xenabrowser.net/datapages/?cohort=GDC%20TCGA%20Liver%20Cancer%20(LIHC)
download the files by hand and put them in data/raw/ with the names below.
"""

from pathlib import Path

import numpy as np
import pandas as pd
import requests

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data" / "raw"
PROCESSED = ROOT / "data" / "processed"

HUB = "https://gdc-hub.s3.us-east-1.amazonaws.com/download"
FILES = {
    "counts": "TCGA-LIHC.star_counts.tsv.gz",
    "clinical": "TCGA-LIHC.clinical.tsv.gz",
    "survival": "TCGA-LIHC.survival.tsv.gz",
    "probemap": "gencode.v36.annotation.gtf.gene.probemap",
}


def download_all(force: bool = False) -> None:
    """Download every file in FILES to data/raw/ (skips files that exist)."""
    RAW.mkdir(parents=True, exist_ok=True)
    for name in FILES.values():
        out = RAW / name
        if out.exists() and not force:
            print(f"✔ {name} already exists")
            continue
        url = f"{HUB}/{name}"
        print(f"↓ {url}")
        try:
            with requests.get(url, stream=True, timeout=120) as r:
                r.raise_for_status()
                with open(out, "wb") as f:
                    for chunk in r.iter_content(chunk_size=1 << 20):
                        f.write(chunk)
        except requests.RequestException as e:
            out.unlink(missing_ok=True)
            print(f"✘ failed: {e}\n  Download it manually (see module docstring).")


# ---------------------------------------------------------------- loaders

def load_counts() -> pd.DataFrame:
    """Genes x samples matrix of raw integer counts.

    Xena stores STAR counts as log2(count + 1), so we undo that transform.
    DESeq2 needs raw integer counts, not log values.
    """
    df = pd.read_csv(RAW / FILES["counts"], sep="\t", index_col=0)
    df.index.name = "gene_id"
    counts = np.round(np.power(2.0, df) - 1).clip(lower=0).astype(np.int64)
    return counts


def load_gene_map() -> pd.Series:
    """Ensembl gene ID (with version) -> gene symbol."""
    pm = pd.read_csv(RAW / FILES["probemap"], sep="\t")
    return pm.set_index("id")["gene"]


def load_survival() -> pd.DataFrame:
    return pd.read_csv(RAW / FILES["survival"], sep="\t")


def load_clinical() -> pd.DataFrame:
    return pd.read_csv(RAW / FILES["clinical"], sep="\t", low_memory=False)


# ---------------------------------------------------------------- samples

SAMPLE_TYPES = {"01": "Tumor", "02": "Recurrent", "06": "Metastatic", "11": "Normal"}


def sample_table(sample_ids) -> pd.DataFrame:
    """Parse TCGA barcodes like TCGA-2V-A95S-01A into patient / type columns."""
    s = pd.Series(list(sample_ids), name="sample")
    meta = pd.DataFrame({
        "sample": s,
        "patient": s.str[:12],
        "type_code": s.str[13:15],
        "vial": s.str[15:16],
    })
    meta["condition"] = meta["type_code"].map(SAMPLE_TYPES).fillna("Other")
    return meta.set_index("sample")


if __name__ == "__main__":
    download_all()
