# Sweep across range of tempratures measure avarage energy vs temprature

import random
import math

from hp_model import energy
from anneal import MOVES, random_initial_conformation


def run_at_temperature(sequence, temperature, equilibration_steps=5000,
                        sample_steps=5000, sample_every=20, seed=None):

    # At fixed temprature, discard the initial equilibration period and sample the energy periodically
    # And return the avarge over those samples
    if seed is not None:
        random.seed(seed)

    coords = random_initial_conformation(len(sequence))
    current_e = energy(sequence, coords)

    def metropolis_step(coords, current_e):
        move_fn = random.choice(MOVES)
        candidate = move_fn(coords)
        if candidate is None:
            return coords, current_e
        cand_e = energy(sequence, candidate)
        delta = cand_e - current_e
        if delta <= 0 or random.random() < math.exp(-delta / temperature):
            return candidate, cand_e
        return coords, current_e

    for _ in range(equilibration_steps):
        coords, current_e = metropolis_step(coords, current_e)

    samples = []
    for step in range(sample_steps):
        coords, current_e = metropolis_step(coords, current_e)
        if step % sample_every == 0:
            samples.append(current_e)

    return sum(samples) / len(samples)


def temperature_sweep(sequence, temperatures, n_repeats=5):

    # For each temperature, average <E> over several independent runs
    # smooths run-to-run noise.
    
    results = {}
    for t in temperatures:
        vals = [
            run_at_temperature(sequence, t, seed=seed)
            for seed in range(n_repeats)
        ]
        results[t] = sum(vals) / len(vals)
    return results


if __name__ == "__main__":
    seq = "HPHPPHHPHPPHPHHPPHPH"
    temperatures = [0.1, 0.3, 0.5, 0.75, 1.0, 1.5, 2.0, 3.0, 5.0]
    results = temperature_sweep(seq, temperatures, n_repeats=5)

    print(f"sequence: {seq}")
    print(f"{'T':>6} | {'<E>':>8}")
    for t in temperatures:
        print(f"{t:>6.2f} | {results[t]:>8.3f}")