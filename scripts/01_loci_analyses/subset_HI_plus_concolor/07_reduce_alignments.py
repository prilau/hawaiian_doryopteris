#!/usr/bin/env python3
"""
reduce_alignments.py

Take the per-locus FASTA alignments in an input directory and, for each one:

  1. Drop the samples D_decipiens_Hawaii_17 and D_decora_Hawaii_6 wherever they
     appear.
  2. Rename every remaining sample using metadata/combined_sample_data.csv:
     "<new_ID>_<coll_num>" with all whitespace in coll_num collapsed to single
     underscores (e.g. D_decipiens_Hawaii_12 -> D_subdecipiens_KR_Wood_10713).
  3. Keep only the wanted samples. By default that is any sample whose new_ID is
     one of the target species (D_concolor, D_decipiens, D_takeuchii, D_angelica,
     D_decora). With --keep-list FILE, keep instead exactly the renamed IDs
     listed in FILE (one per line; '/' is normalised to '_').
  4. Re-align the surviving sequences with MAFFT and write the result to the
     output directory (original files are left untouched).

Usage:
    python3 reduce_alignments.py
    python3 reduce_alignments.py --input-dir data_retrieved --output-dir data_reduced
    python3 reduce_alignments.py --output-dir data_reduced_50_with_subdecip \
                                 --keep-list metadata/samples_50_with_subdecip.txt
"""

import argparse
import csv
import re
import subprocess
import sys
import tempfile
from pathlib import Path

# --- fixed parameters from the task ------------------------------------------

DROP_SAMPLES = {"D_decipiens_Hawaii_17", "D_decora_Hawaii_6"}

KEEP_SPECIES = {
    "D_concolor",
    "D_decipiens",
    "D_takeuchii",
    "D_angelica",
    "D_decora",
}


# --- helpers ----------------------------------------------------------------

def norm(s):
    return s.replace("/", "_")


def load_rename_map(metadata_path):
    """Return {old_name: (new_name, species)} built from the metadata CSV."""
    rename = {}
    with open(metadata_path, newline="", encoding="utf-8-sig") as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            old = row["name"].strip()
            species = row["new_ID"].strip()
            coll = re.sub(r"\s+", "_", row["coll_num"].strip())
            new = f"{species}_{coll}"
            rename[old] = (new, species)
    return rename


def read_fasta(path):
    """Yield (header, sequence) pairs in file order."""
    header, chunks = None, []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(chunks)
                header, chunks = line[1:].strip(), []
            elif line:
                chunks.append(line)
    if header is not None:
        yield header, "".join(chunks)


def write_fasta(path, records):
    with open(path, "w") as fh:
        for header, seq in records:
            fh.write(f">{header}\n{seq}\n")


def run_mafft(in_path, out_path, extra_opts, threads):
    cmd = ["mafft", "--preservecase", "--thread", str(threads)]
    cmd += extra_opts
    cmd.append(str(in_path))
    with open(out_path, "w") as out:
        proc = subprocess.run(cmd, stdout=out, stderr=subprocess.PIPE, text=True)
    if proc.returncode != 0:
        sys.stderr.write(proc.stderr)
        raise RuntimeError(f"mafft failed on {in_path} (exit {proc.returncode})")


# --- main -----------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--input-dir", default="data_retrieved", type=Path)
    ap.add_argument("--output-dir", default="data_reduced", type=Path)
    ap.add_argument("--metadata", default="metadata/combined_sample_data.csv", type=Path)
    ap.add_argument("--keep-list", type=Path, default=None,
                    help="file of renamed IDs to keep (one per line); overrides "
                         "the default species filter")
    ap.add_argument("--mafft-opts", default="--auto",
                    help="options passed through to mafft (default: --auto)")
    ap.add_argument("--threads", default="-1",
                    help="mafft --thread value (default: -1 = autodetect)")
    args = ap.parse_args()

    rename = load_rename_map(args.metadata)
    keep_ids = None
    if args.keep_list:
        keep_ids = {norm(l.strip()) for l in open(args.keep_list) if l.strip()}
    extra_opts = args.mafft_opts.split()

    args.output_dir.mkdir(parents=True, exist_ok=True)
    fasta_files = sorted(args.input_dir.glob("*.fasta"))
    if not fasta_files:
        sys.exit(f"No .fasta files found in {args.input_dir}")

    n_written = n_skipped = 0
    unknown = set()

    for fp in fasta_files:
        kept = []
        seen_new = set()
        for header, seq in read_fasta(fp):
            if header in DROP_SAMPLES:                       # step 1
                continue
            if header not in rename:
                unknown.add(header)
                continue
            new_name, species = rename[header]              # step 2
            if keep_ids is not None:                        # step 3
                if norm(new_name) not in keep_ids:
                    continue
            elif species not in KEEP_SPECIES:
                continue
            if new_name in seen_new:
                sys.stderr.write(
                    f"  {fp.name}: duplicate name {new_name} after rename, "
                    f"keeping first only\n")
                continue
            seen_new.add(new_name)
            degapped = seq.replace("-", "").replace(".", "")
            kept.append((new_name, degapped))

        out_path = args.output_dir / fp.name

        if len(kept) < 2:                                   # nothing to align
            if kept:
                write_fasta(out_path, kept)
                n_written += 1
                print(f"{fp.name}: {len(kept)} seq, not realigned (n<2)")
            else:
                n_skipped += 1
                print(f"{fp.name}: 0 seqs left, skipped")
            continue

        with tempfile.NamedTemporaryFile("w", suffix=".fasta", delete=False) as tmp:
            tmp_path = Path(tmp.name)
        try:
            write_fasta(tmp_path, kept)
            run_mafft(tmp_path, out_path, extra_opts, args.threads)  # step 4
        finally:
            tmp_path.unlink(missing_ok=True)
        n_written += 1
        print(f"{fp.name}: {len(kept)} seq -> realigned")

    print(f"\nDone. {n_written} alignments written to {args.output_dir}, "
          f"{n_skipped} skipped (empty).")
    if unknown:
        sys.stderr.write("\nHeaders with no metadata match (left out):\n")
        for h in sorted(unknown):
            sys.stderr.write(f"  {h}\n")


if __name__ == "__main__":
    main()
