#!/usr/bin/env python3
"""
Standalone script: plot xcorr vs q-value dotplots per peptide.

Reads config from a YAML file (default: config/config_lumos_1.yaml),
extracts xcorr (COMET:xcorr) and Percolator q-value from idXML files,
and produces one dotplot per peptide with a linear regression line and
Pearson correlation annotation.

Only PSMs from Percolator-scored files (score type == "q-value") are used,
so raw Comet (psm.idXML) and MS2Rescore output (ms2rescore.idXML) are
automatically skipped.  PSMs are deduplicated by (spectrum_id, sequence)
so pout.idXML and pout_filtered.idXML do not double-count the same hit.

Usage:
    conda run -n mspecviz python scripts/plot_xcorr_vs_qvalue.py
    conda run -n mspecviz python scripts/plot_xcorr_vs_qvalue.py \\
        --config config/config_lumos_1.yaml \\
        --out /path/to/plots \\
        --workers 15
"""
import argparse
import glob
import logging
import sys
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import yaml

# --- Project root on sys.path so we can import from src/ ---
_SCRIPT_DIR = Path(__file__).resolve().parent
_PROJECT_ROOT = _SCRIPT_DIR.parent
if str(_PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(_PROJECT_ROOT))

from src.export import sanitize_filename  # noqa: E402

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%H:%M:%S",
)
logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Config loading
# ---------------------------------------------------------------------------

def load_config(config_path: str) -> dict:
    """Load and return the YAML config dict."""
    with open(config_path) as f:
        return yaml.safe_load(f)


def resolve_globs(patterns: List[str]) -> List[str]:
    """Expand a list of glob patterns into sorted unique file paths."""
    paths: Set[str] = set()
    for pattern in patterns:
        paths.update(glob.glob(pattern, recursive=True))
    return sorted(paths)


# ---------------------------------------------------------------------------
# idXML extraction — module-level worker (required for ProcessPoolExecutor)
# ---------------------------------------------------------------------------

def _worker_extract_one(
    idxml_path: str,
    target_set: Set[str],
) -> Tuple[str, Dict[str, List[Tuple[float, float, str]]], str]:
    """
    Worker: extract (xcorr, q-value, dedup_key) triples from one idXML.

    Only processes PeptideIdentifications whose score type is "q-value"
    (i.e. Percolator-rescored output).  Skips decoy hits.

    Must be a module-level function for ProcessPoolExecutor pickling.

    Returns:
        (idxml_path,
         {base_seq: [(xcorr, q_val, dedup_key), ...]},
         error_msg)
    """
    import pyopenms as oms

    local: Dict[str, List[Tuple[float, float, str]]] = {}
    error = ""
    try:
        protein_ids: list = []
        peptide_ids: list = []
        oms.IdXMLFile().load(str(idxml_path), protein_ids, peptide_ids)

        for pep_id in peptide_ids:
            # Only use PSMs where Percolator assigned a q-value
            score_type = pep_id.getScoreType().lower()
            if score_type != "q-value":
                continue

            # Spectrum native ID for deduplication across pout/pout_filtered
            spectrum_ref = ""
            if pep_id.metaValueExists("spectrum_reference"):
                rv = pep_id.getMetaValue("spectrum_reference")
                spectrum_ref = rv.decode() if isinstance(rv, bytes) else str(rv)

            for hit in pep_id.getHits():
                try:
                    # Skip decoys
                    if hit.metaValueExists("target_decoy"):
                        td = hit.getMetaValue("target_decoy")
                        if isinstance(td, bytes):
                            td = td.decode()
                        if "decoy" in str(td).lower():
                            continue

                    base_seq = hit.getSequence().toUnmodifiedString().upper()
                    if base_seq not in target_set:
                        continue

                    q_value = float(hit.getScore())

                    # xcorr: prefer COMET:xcorr, fall back to PSI-MS accession
                    xcorr: Optional[float] = None
                    for key in ("COMET:xcorr", "MS:1002252"):
                        if hit.metaValueExists(key):
                            try:
                                xcorr = float(hit.getMetaValue(key))
                                break
                            except (TypeError, ValueError):
                                pass

                    if xcorr is None:
                        continue

                    # Dedup key: file stem (without _pout/_pout_filtered) + spectrum
                    stem = Path(idxml_path).stem
                    for suffix in ("_pout_filtered", "_pout"):
                        if stem.endswith(suffix):
                            stem = stem[: -len(suffix)]
                            break
                    dedup_key = f"{stem}|{spectrum_ref}|{base_seq}"

                    if base_seq not in local:
                        local[base_seq] = []
                    local[base_seq].append((xcorr, q_value, dedup_key))

                except Exception:
                    continue

    except Exception as exc:
        error = f"{Path(idxml_path).name}: {exc}"

    return idxml_path, local, error


def extract_xcorr_qvalue(
    idxml_paths: List[str],
    target_peptides: List[str],
    workers: int = 1,
    progress_callback=None,
) -> Dict[str, List[Tuple[float, float]]]:
    """
    Extract (xcorr, q-value) pairs for each target peptide from idXML files.

    PSMs are deduplicated so the same spectrum from pout.idXML and
    pout_filtered.idXML is counted only once.

    Args:
        idxml_paths:       List of idXML file paths to scan.
        target_peptides:   Base sequences to match (case-insensitive).
        workers:           Parallel workers (1 = sequential).
        progress_callback: Optional callable(current, total, filename).

    Returns:
        Dict mapping base_sequence -> [(xcorr, q_value), ...].
    """
    from concurrent.futures import ProcessPoolExecutor, as_completed

    target_set: Set[str] = {p.strip().upper() for p in target_peptides if p.strip()}
    # Use sets per sequence to deduplicate by dedup_key
    seen: Dict[str, Set[str]] = {seq: set() for seq in target_set}
    results: Dict[str, List[Tuple[float, float]]] = {seq: [] for seq in target_set}

    def _merge(local: Dict[str, List[Tuple[float, float, str]]], error: str) -> None:
        if error:
            logger.warning(error)
        for seq, triples in local.items():
            if seq not in results:
                continue
            for xcorr, q_val, dedup_key in triples:
                if dedup_key not in seen[seq]:
                    seen[seq].add(dedup_key)
                    results[seq].append((xcorr, q_val))

    total = len(idxml_paths)
    workers = max(1, min(workers, total))

    if workers <= 1:
        for idx, path in enumerate(idxml_paths):
            if progress_callback:
                progress_callback(idx, total, Path(path).name)
            _, local, error = _worker_extract_one(path, target_set)
            _merge(local, error)
    else:
        with ProcessPoolExecutor(max_workers=workers) as pool:
            futures = {
                pool.submit(_worker_extract_one, p, target_set): p
                for p in idxml_paths
            }
            done = 0
            for future in as_completed(futures):
                done += 1
                path = futures[future]
                try:
                    _, local, error = future.result()
                    _merge(local, error)
                except Exception as exc:
                    logger.error(f"Worker crashed on {Path(path).name}: {exc}")
                if progress_callback:
                    progress_callback(done, total, Path(path).name)

    if progress_callback:
        progress_callback(total, total, "Done")

    return results


# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------

def plot_peptide(
    peptide: str,
    pairs: List[Tuple[float, float]],
    out_dir: Path,
) -> Optional[Path]:
    """
    Create a xcorr vs -log10(q-value) dotplot with regression line.

    Using -log10(q-value) on the y-axis spreads the low-q-value region
    (the scientifically interesting PSMs) and compresses the high-q noise.

    Args:
        peptide:  Base sequence (title and filename stem).
        pairs:    List of (xcorr, q_value) tuples.
        out_dir:  Output directory.

    Returns:
        Path to saved PNG, or None if fewer than 2 distinct xcorr values.
    """
    from scipy import stats

    if len(pairs) < 2:
        logger.info(f"Skipping {peptide}: only {len(pairs)} PSM(s), need >= 2")
        return None

    xcorr = np.array([p[0] for p in pairs])
    qvals = np.array([p[1] for p in pairs])

    # Guard against q-value == 0 before log transform (clip to small epsilon)
    eps = 1e-10
    neg_log_q = -np.log10(np.clip(qvals, eps, 1.0))

    # Guard against zero variance in xcorr (all identical values)
    if np.std(xcorr) == 0:
        logger.info(f"Skipping regression for {peptide}: all xcorr values identical")
        reg_label = "xcorr has no variance — regression not applicable"
        show_regression = False
    else:
        slope, intercept, r_value, p_value, _ = stats.linregress(xcorr, neg_log_q)
        pearson_r, pearson_p = stats.pearsonr(xcorr, neg_log_q)
        r_sq = r_value ** 2
        reg_label = f"r = {pearson_r:.3f}  (p = {pearson_p:.2e})\nR² = {r_sq:.3f}"
        show_regression = True

    fig, ax = plt.subplots(figsize=(7, 5))

    ax.scatter(xcorr, neg_log_q, alpha=0.4, s=15, color="#2171b5", linewidths=0)

    if show_regression:
        x_line = np.linspace(xcorr.min(), xcorr.max(), 200)
        ax.plot(
            x_line, slope * x_line + intercept,
            color="#d62728", linewidth=1.5, label=reg_label,
        )
        ax.legend(fontsize=9)
    else:
        ax.text(
            0.05, 0.95, reg_label, transform=ax.transAxes,
            fontsize=8, va="top", color="gray",
        )

    # Reference line at q = 0.01 (common significance threshold)
    q_thresh = 0.01
    ax.axhline(-np.log10(q_thresh), color="grey", linestyle="--",
               linewidth=0.8, label=f"q = {q_thresh}")

    ax.set_xlabel("Comet xcorr", fontsize=12)
    ax.set_ylabel("−log₁₀(q-value)", fontsize=12)
    ax.set_title(f"{peptide}  (n = {len(pairs)})", fontsize=13)
    ax.grid(True, alpha=0.3)
    fig.tight_layout()

    stem = sanitize_filename(peptide)
    out_path = out_dir / f"{stem}.png"
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    return out_path


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    parser = argparse.ArgumentParser(
        description="Plot Comet xcorr vs Percolator q-value per peptide"
    )
    parser.add_argument(
        "--config",
        default=str(_PROJECT_ROOT / "config" / "config_lumos_1.yaml"),
        help="Path to YAML config file",
    )
    parser.add_argument(
        "--out",
        default=None,
        help="Output directory for PNG plots (default: exports/xcorr_qvalue/)",
    )
    parser.add_argument(
        "--workers",
        type=int,
        default=None,
        help="Parallel workers for idXML scanning (overrides config value)",
    )
    args = parser.parse_args()

    # --- Load config ---
    logger.info(f"Loading config from {args.config}")
    cfg = load_config(args.config)

    idxml_patterns: List[str] = cfg.get("input", {}).get("idxml_files", [])
    peptide_file: str = cfg.get("peptides", "")
    workers: int = args.workers or cfg.get("performance", {}).get("workers", 1)
    out_dir_str: str = (
        args.out or str(_PROJECT_ROOT / "exports" / "xcorr_qvalue")
    )

    # --- Resolve idXML globs ---
    logger.info("Resolving idXML file globs...")
    idxml_paths = resolve_globs(idxml_patterns)
    if not idxml_paths:
        logger.error("No idXML files found — check glob patterns in config")
        sys.exit(1)
    logger.info(f"Found {len(idxml_paths)} idXML file(s)")

    # --- Load target peptides ---
    if not peptide_file or not Path(peptide_file).exists():
        logger.error(f"Peptide file not found: {peptide_file}")
        sys.exit(1)
    with open(peptide_file) as fh:
        peptides = [line.strip() for line in fh if line.strip()]
    logger.info(f"Loaded {len(peptides)} peptide(s) from {peptide_file}")

    # --- Prepare output directory ---
    out_dir = Path(out_dir_str)
    out_dir.mkdir(parents=True, exist_ok=True)
    logger.info(f"Output directory: {out_dir}")

    # --- Extract xcorr + q-value pairs ---
    logger.info(
        f"Scanning idXML files with {workers} worker(s) "
        "(Percolator q-value files only, deduplicated)..."
    )
    log_every = max(1, len(idxml_paths) // 20)

    def _progress(current: int, total: int, msg: str) -> None:
        if current % log_every == 0 or current == total:
            logger.info(f"  [{current}/{total}] {msg}")

    results = extract_xcorr_qvalue(
        idxml_paths, peptides, workers=workers, progress_callback=_progress
    )

    total_psms = sum(len(v) for v in results.values())
    logger.info(
        f"Extracted {total_psms} unique PSM pair(s) across {len(peptides)} peptide(s)"
    )

    not_found = [p for p in peptides if not results.get(p.upper())]
    if not_found:
        logger.warning(
            f"No PSMs found for {len(not_found)} peptide(s): "
            + ", ".join(not_found[:10])
            + ("..." if len(not_found) > 10 else "")
        )

    # --- Plot one file per peptide ---
    saved = 0
    for pep in peptides:
        pairs = results.get(pep.upper(), [])
        if not pairs:
            continue
        out_path = plot_peptide(pep.upper(), pairs, out_dir)
        if out_path:
            logger.info(f"  Saved {out_path.name}  ({len(pairs)} PSMs)")
            saved += 1

    logger.info(f"Done — {saved} plot(s) saved to {out_dir}")


if __name__ == "__main__":
    main()
