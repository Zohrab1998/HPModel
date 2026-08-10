"""
Command-line entry point for the HP model folding exercise.

Usage examples:
    python3 run.py --sequence HPPHPH --method exact
    python3 run.py --sequence HPHPPHHPHPPHPHHPPHPH --method anneal
    python3 run.py --sequence HPHPPHHPHPPHPHHPPHPH --method anneal --runs 200

--method exact   : exhaustive enumeration. Only feasible for short
                   sequences (roughly length <= 14-15 on a normal
                   machine; grows ~7x for every +2 residues).
--method anneal  : simulated annealing with multiple random restarts.
                   Reports the best energy found and the distribution
                   of per-run results, which is the evidence for
                   whether the true ground state was reached.
"""

import argparse
import time
from collections import Counter

from hp_model import energy
from enumerate_exact import find_ground_state
from anneal import multi_restart
from visualize import render


def run_exact(sequence):
    print(f"Method: exhaustive enumeration (exact)")
    t0 = time.time()
    best_e, confs = find_ground_state(sequence)
    t1 = time.time()

    print(f"Sequence: {sequence}  (length {len(sequence)})")
    print(f"Ground-state energy: {best_e}")
    print(f"Conformations achieving it (raw count, includes lattice "
          f"symmetry duplicates): {len(confs)}")
    print(f"Search time: {t1 - t0:.2f}s")
    print()
    print("(letters = residues, '-'/'|' = chain bonds, "
          "'o' = non-bonded H-H contact contributing to the energy)")
    print(render(sequence, confs[0]))


def run_anneal(sequence, n_runs):
    print(f"Method: simulated annealing, {n_runs} independent restarts")
    t0 = time.time()
    best_e, best_coords, per_run = multi_restart(sequence, n_runs=n_runs)
    t1 = time.time()

    # Independent cross-check: recompute the energy of the reported
    # best conformation directly from its coordinates, rather than
    # trusting the value the annealer tracked internally.
    verified_e = energy(sequence, best_coords)
    assert verified_e == best_e, "internal energy tracking mismatch!"

    dist = dict(sorted(Counter(per_run).items()))
    hits = dist.get(best_e, 0)

    print(f"Sequence: {sequence}  (length {len(sequence)})")
    print(f"Best energy found: {best_e}  (independently verified from coords)")
    print(f"Energy distribution across {n_runs} restarts: {dist}")
    print(f"Fraction of runs reaching the best energy: {hits}/{n_runs}")
    print(f"Search time: {t1 - t0:.1f}s")
    print()
    print("(letters = residues, '-'/'|' = chain bonds, "
          "'o' = non-bonded H-H contact contributing to the energy)")
    print(render(sequence, best_coords))


def main():
    parser = argparse.ArgumentParser(
        description="Fold an HP sequence on the 2D square lattice.")
    parser.add_argument(
        "--sequence", default="HPPHPH",
        help="HP sequence, e.g. HPPHPH (default: %(default)s)")
    parser.add_argument(
        "--method", choices=["exact", "anneal"], default="exact",
        help="'exact' = exhaustive enumeration (short sequences only), "
             "'anneal' = simulated annealing with restarts (default: %(default)s)")
    parser.add_argument(
        "--runs", type=int, default=100,
        help="number of independent annealing restarts (default: %(default)s)")
    args = parser.parse_args()

    sequence = args.sequence.upper()
    if not set(sequence) <= {"H", "P"}:
        raise ValueError("sequence must contain only 'H' and 'P' characters")

    if args.method == "exact":
        run_exact(sequence)
    else:
        run_anneal(sequence, args.runs)


if __name__ == "__main__":
    main()