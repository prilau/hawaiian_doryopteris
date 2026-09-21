#!/usr/bin/env python3
"""
Subset the per-locus alignments in data_reduced/ to only the samples present in
a keep-list (the VCF sample columns), then MAFFT-realign each locus.

    python3 subset_to_vcf.py
    python3 subset_to_vcf.py --input-dir data_reduced --output-dir data_reduced_vcf \
                             --keep-list metadata/vcf_samples.txt

Sample names are compared with '/' normalised to '_' (IQ-TREE rewrites '/' in
taxon labels, so downstream names use the underscore form).
"""

import argparse
import subprocess
import sys
import tempfile
from pathlib import Path


def norm(s):
    return s.replace("/", "_")


def read_fasta(path):
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
    cmd = ["mafft", "--preservecase", "--thread", str(threads)] + extra_opts + [str(in_path)]
    with open(out_path, "w") as out:
        proc = subprocess.run(cmd, stdout=out, stderr=subprocess.PIPE, text=True)
    if proc.returncode != 0:
        sys.stderr.write(proc.stderr)
        raise RuntimeError(f"mafft failed on {in_path} (exit {proc.returncode})")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--input-dir", default="data_reduced", type=Path)
    ap.add_argument("--output-dir", default="data_reduced_vcf", type=Path)
    ap.add_argument("--keep-list", default="metadata/vcf_samples.txt", type=Path)
    ap.add_argument("--mafft-opts", default="--auto")
    ap.add_argument("--threads", default="-1")
    args = ap.parse_args()

    keep = {norm(l.strip()) for l in open(args.keep_list) if l.strip()}
    extra_opts = args.mafft_opts.split()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    fastas = sorted(args.input_dir.glob("*.fasta"))
    if not fastas:
        sys.exit(f"No .fasta files in {args.input_dir}")

    n_written = n_skipped = 0
    for fp in fastas:
        kept = [(h, seq.replace("-", "").replace(".", ""))
                for h, seq in read_fasta(fp) if norm(h) in keep]
        out_path = args.output_dir / fp.name

        if len(kept) < 2:
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
            run_mafft(tmp_path, out_path, extra_opts, args.threads)
        finally:
            tmp_path.unlink(missing_ok=True)
        n_written += 1
        print(f"{fp.name}: {len(kept)} seq -> realigned")

    print(f"\nDone. {n_written} alignments in {args.output_dir}, {n_skipped} skipped (empty).")


if __name__ == "__main__":
    main()
