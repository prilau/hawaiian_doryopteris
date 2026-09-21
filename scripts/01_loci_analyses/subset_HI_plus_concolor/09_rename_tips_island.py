#!/usr/bin/env python3
"""
Append each tip's collection locality (island) to its label in a Newick tree.

    python3 rename_tips_island.py wastral_all.tre wastral_all.island.tre

Mapping comes from metadata/combined_sample_data.csv: a tip named
"<new_ID>_<coll_num>" gets "_<locality>" appended, with every run of
non-alphanumeric characters in coll_num / locality turned into a single "_"
(matching how the names were sanitised upstream).
"""
import csv
import re
import sys

META = "metadata/combined_sample_data.csv"


def sanitize(s):
    """For the appended locality: any non-alphanumeric run -> single '_'."""
    return re.sub(r"_+", "_", re.sub(r"[^A-Za-z0-9]+", "_", s.strip())).strip("_")


def tip_name(new_id, coll_num):
    """Reproduce the tip label as it ends up in the IQ-TREE treefiles:
    reduce_alignments.py collapses whitespace to '_', then IQ-TREE rewrites '/'."""
    return new_id.strip() + "_" + re.sub(r"\s+", "_", coll_num.strip()).replace("/", "_")


def main():
    src, dst = sys.argv[1], sys.argv[2]

    loc = {}
    with open(META, newline="", encoding="utf-8-sig") as fh:
        for row in csv.DictReader(fh):
            loc[tip_name(row["new_ID"], row["coll_num"])] = sanitize(row["locality"])

    newick = open(src).read()
    tips = set(re.findall(r"[(,]([^(),:]+)(?=[,:)])", newick))

    renamed, missing = 0, []
    for tip in sorted(tips, key=len, reverse=True):
        island = loc.get(tip)
        if not island:
            missing.append(tip)
            continue
        newick = re.sub(r"(?<=[(,])" + re.escape(tip) + r"(?=[,:)])",
                        f"{tip}_{island}", newick)
        renamed += 1

    with open(dst, "w") as out:
        out.write(newick)

    print(f"{renamed} tips relabelled -> {dst}")
    if missing:
        print("no locality found for:", *missing, sep="\n  ")


if __name__ == "__main__":
    main()
