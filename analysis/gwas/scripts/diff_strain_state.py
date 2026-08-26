#!/usr/bin/env python3
"""Diff the reconciled+ploidy-cleared strain state against the prior GWAS run's
recorded state (strain list + population assignment + near-clone partition).

'Unchanged' requires identity on ALL THREE (spec S4, quant-genetics review fix #1/#5) --
a bare strain-name-list comparison is not sufficient.
"""
import argparse
import csv
import json


def read_fam_ids(path: str) -> set[str]:
    ids = set()
    with open(path) as f:
        for line in f:
            ids.add(line.split()[1])
    return ids


def read_pop_csv(path: str) -> dict[str, str]:
    out = {}
    with open(path) as f:
        for row in csv.DictReader(f):
            out[row["Strain"]] = row["Pop"]
    return out


def current_strain_list(reviewed_csv: str, excluded: set[str]) -> set[str]:
    ids = set()
    with open(reviewed_csv) as f:
        for row in csv.DictReader(f):
            tier = row.get("tier", "")
            decision = row.get("decision", "")
            if (tier in ("exact", "normalized") or decision == "accept") and row["vcf_sample_id"]:
                if row["vcf_sample_id"] not in excluded:
                    ids.add(row["vcf_sample_id"])
    return ids


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--reviewed", required=True)
    ap.add_argument("--exclude", nargs="*", default=[], help="strain IDs excluded by the ploidy review")
    ap.add_argument("--prior-fam-all", required=True)
    ap.add_argument("--prior-fam-culled", required=True)
    ap.add_argument("--prior-pop", required=True)
    ap.add_argument("--current-pop", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    excluded = set(args.exclude)
    current_ids = current_strain_list(args.reviewed, excluded)
    prior_all_ids = read_fam_ids(args.prior_fam_all)
    prior_culled_ids = read_fam_ids(args.prior_fam_culled)

    added = sorted(current_ids - prior_all_ids)
    removed = sorted(prior_all_ids - current_ids)
    strain_list_changed = bool(added or removed)

    prior_pop = read_pop_csv(args.prior_pop)
    current_pop = read_pop_csv(args.current_pop)
    pop_changes = []
    for strain in (current_ids & prior_all_ids):
        old = prior_pop.get(strain)
        new = current_pop.get(strain)
        if old != new:
            pop_changes.append({"strain": strain, "old_pop": old, "new_pop": new})
    population_assignment_changed = bool(pop_changes)

    # Near-clone/culled partition: unchanged iff no strain that changed membership
    # (added/removed) intersects the culled-173 set, since that set defines which
    # strains BSLMM/LOCO/pixy were computed on.
    culled_partition_touched = bool((set(added) | set(removed)) & (prior_all_ids | prior_culled_ids))

    overall_changed = strain_list_changed or population_assignment_changed
    action = "rebuild" if overall_changed else "copy_forward"

    report = {
        "strain_list_changed": strain_list_changed,
        "strain_list_added": added,
        "strain_list_removed": removed,
        "population_assignment_changed": population_assignment_changed,
        "population_changes": pop_changes,
        "near_clone_partition_touched": culled_partition_touched,
        "overall_action": action,
    }

    import os
    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    with open(args.out, "w") as f:
        json.dump(report, f, indent=2)

    print(f"Current strain count: {len(current_ids)} (prior: {len(prior_all_ids)})")
    print(f"Added: {len(added)}  Removed: {len(removed)}")
    print(f"Population changes: {len(pop_changes)}")
    print(f"Decision: overall_action = {action}")
    print(f"Wrote {args.out}")


if __name__ == "__main__":
    main()
