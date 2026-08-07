#!/usr/bin/env python3
"""
Assign reads to BLAZE cell barcodes (BLAZE step 3) for downstream BAM tagging.

Reads putative_bc.csv and whitelist.csv produced by BLAZE, applies the same
read-to-whitelist matching as blaze.read_assignment (subsequence edit distance,
max ED 2 by default), and writes a TSV of corrected assignments:

    read_id<TAB>CB<TAB>UB

Only successfully assigned reads are written. Unassigned and ambiguous reads
are omitted (matching BLAZE demultiplex behaviour).

This is the expensive step for large libraries: non-exact barcodes scan the
full whitelist with edit distance. Use --threads and watch --progress-log
(or stderr) for throughput.
"""

from __future__ import annotations

import argparse
import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from typing import Any

import pandas as pd
from blaze.config import DEFAULT_ASSIGNMENT_ED
from blaze.read_assignment import _match_bc_row

# Columns required by blaze.read_assignment._match_bc_row
_REQUIRED_COLS = (
    "read_id",
    "putative_bc",
    "putative_umi",
    "putative_bc_qscore",
    "polyT_end",
    "pre_bc_flanking",
    "post_umi_flanking",
)


def load_whitelist(path: str) -> set[str]:
    barcodes: list[str] = []
    with open(path, encoding="utf-8") as handle:
        for line in handle:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            # BLAZE whitelist.csv and Cell Ranger exports may include GEM suffix.
            token = line.split(",")[0].strip()
            barcodes.append(token.split("-")[0])
    return set(barcodes)


def _assign_chunk(
    frame: pd.DataFrame,
    whitelist: set[str],
    max_ed: int,
    min_q: int,
) -> list[tuple[str, str, str]]:
    """Worker: assign one chunk of putative_bc rows."""
    out: list[tuple[str, str, str]] = []
    frame = frame.fillna("")
    for row in frame.itertuples(index=False):
        bc, umi, _strand = _match_bc_row(row, whitelist, max_ed, min_q)
        if bc:
            out.append((str(row.read_id), str(bc), str(umi)))
    return out


def _log(msg: str, progress_log: str | None) -> None:
    line = msg if msg.endswith("\n") else msg + "\n"
    sys.stderr.write(line)
    sys.stderr.flush()
    if progress_log:
        with open(progress_log, "a", encoding="utf-8") as handle:
            handle.write(line)


def assign_reads(
    putative_bc_csv: str,
    whitelist_csv: str,
    output_tsv: str,
    max_ed: int = DEFAULT_ASSIGNMENT_ED,
    min_q: int = 0,
    threads: int = 1,
    chunksize: int = 50000,
    progress_log: str | None = None,
) -> tuple[int, int]:
    whitelist = load_whitelist(whitelist_csv)
    _log(
        "Loaded whitelist: {} barcodes from {}".format(len(whitelist), whitelist_csv),
        progress_log,
    )

    # Probe columns once
    header = pd.read_csv(putative_bc_csv, nrows=0)
    missing = [c for c in _REQUIRED_COLS if c not in header.columns]
    if missing:
        raise SystemExit(
            "putative_bc.csv missing columns required for BLAZE assignment: "
            + ", ".join(missing)
        )

    assigned = 0
    total = 0
    t0 = time.time()
    threads = max(1, int(threads))
    chunksize = max(1000, int(chunksize))

    _log(
        "Assigning reads from {} (threads={}, chunksize={})".format(
            putative_bc_csv, threads, chunksize
        ),
        progress_log,
    )

    reader = pd.read_csv(putative_bc_csv, chunksize=chunksize)

    with open(output_tsv, "w", encoding="utf-8") as out:
        if threads == 1:
            for chunk in reader:
                hits = _assign_chunk(chunk, whitelist, max_ed, min_q)
                for read_id, bc, umi in hits:
                    out.write(f"{read_id}\t{bc}\t{umi}\n")
                assigned += len(hits)
                total += len(chunk)
                rate = total / max(time.time() - t0, 1e-6)
                _log(
                    "Progress: {:,} reads ({:,} assigned, {:.0f} reads/s)".format(
                        total, assigned, rate
                    ),
                    progress_log,
                )
        else:
            # Submit chunks as they are read; bound in-flight work to ~2x threads
            in_flight: dict[Any, int] = {}
            with ProcessPoolExecutor(max_workers=threads) as pool:
                for chunk in reader:
                    fut = pool.submit(
                        _assign_chunk, chunk.copy(), whitelist, max_ed, min_q
                    )
                    in_flight[fut] = len(chunk)
                    if len(in_flight) >= threads * 2:
                        done = next(as_completed(in_flight))
                        n_rows = in_flight.pop(done)
                        hits = done.result()
                        for read_id, bc, umi in hits:
                            out.write(f"{read_id}\t{bc}\t{umi}\n")
                        assigned += len(hits)
                        total += n_rows
                        rate = total / max(time.time() - t0, 1e-6)
                        _log(
                            "Progress: {:,} reads ({:,} assigned, "
                            "{:.0f} reads/s)".format(total, assigned, rate),
                            progress_log,
                        )

                for fut in as_completed(in_flight):
                    n_rows = in_flight[fut]
                    hits = fut.result()
                    for read_id, bc, umi in hits:
                        out.write(f"{read_id}\t{bc}\t{umi}\n")
                    assigned += len(hits)
                    total += n_rows
                    rate = total / max(time.time() - t0, 1e-6)
                    _log(
                        "Progress: {:,} reads ({:,} assigned, "
                        "{:.0f} reads/s)".format(total, assigned, rate),
                        progress_log,
                    )

    elapsed = time.time() - t0
    _log(
        "BLAZE read assignment: {}/{} reads assigned in {:.1f}s "
        "({:.0f} reads/s)".format(
            assigned, total, elapsed, total / max(elapsed, 1e-6)
        ),
        progress_log,
    )
    return assigned, total


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Assign reads to BLAZE whitelist barcodes for BAM CB/UB tagging."
    )
    parser.add_argument(
        "--putative-bc",
        required=True,
        help="BLAZE putative_bc.csv from step 1",
    )
    parser.add_argument(
        "--whitelist",
        required=True,
        help="BLAZE whitelist.csv from step 2",
    )
    parser.add_argument(
        "--output",
        required=True,
        help="Output TSV: read_id, corrected CB, corrected UMI",
    )
    parser.add_argument(
        "--max-ed",
        type=int,
        default=DEFAULT_ASSIGNMENT_ED,
        help=f"Maximum subsequence edit distance (default: {DEFAULT_ASSIGNMENT_ED})",
    )
    parser.add_argument(
        "--min-q",
        type=int,
        default=0,
        help="Minimum putative barcode minQ for assignment (default: 0, match BLAZE demux)",
    )
    parser.add_argument(
        "--threads",
        type=int,
        default=max(1, (os.cpu_count() or 4) // 2),
        help="Parallel worker processes for assignment (default: half CPUs)",
    )
    parser.add_argument(
        "--chunksize",
        type=int,
        default=50000,
        help="Rows per assignment chunk (default: 50000)",
    )
    parser.add_argument(
        "--progress-log",
        default=None,
        help="Optional path to append progress lines (tail -f this while running)",
    )
    args = parser.parse_args(argv)

    if args.progress_log:
        # Truncate at start of a fresh run
        open(args.progress_log, "w", encoding="utf-8").close()

    assign_reads(
        args.putative_bc,
        args.whitelist,
        args.output,
        max_ed=args.max_ed,
        min_q=args.min_q,
        threads=args.threads,
        chunksize=args.chunksize,
        progress_log=args.progress_log,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
