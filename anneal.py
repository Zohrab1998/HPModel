# Simulated annealing search
# as exhaustive enumeration is only for short sequences
""" Added some move set :
    1. End move (NON-ergodic move)
    2. Corner move (NON-ergodic move)
    3. Pivot move (Ergodic move)
"""
import random
import math

from hp_model import DIRECTIONS, energy, is_self_avoiding


def try_end_move(coords):
    # Move one random endpoint -> returns new coords or None (no valid move exists)
    n = len(coords)
    end = random.choice([0, n - 1])
    anchor = 1 if end == 0 else n - 2
    ax, ay = coords[anchor]
    candidates = [(ax + dx, ay + dy) for dx, dy in DIRECTIONS]
    random.shuffle(candidates)
    occupied = set(coords)
    occupied.discard(coords[end])
    for cand in candidates:
        if cand not in occupied:
            new_coords = list(coords)
            new_coords[end] = cand
            return new_coords
    return None


def try_corner_move(coords):
    # Corner flip at a random interior residue.
    n = len(coords)
    if n < 3:
        return None
    i = random.randint(1, n - 2)
    x0, y0 = coords[i - 1]
    x1, y1 = coords[i]
    x2, y2 = coords[i + 1]
    if abs(x0 - x2) != 1 or abs(y0 - y2) != 1:
        return None  # i-1 and i+1 are not diagonal -> not a corner
    opt_a, opt_b = (x0, y2), (x2, y0)
    new_pos = opt_b if (x1, y1) == opt_a else opt_a
    occupied = set(coords)
    occupied.discard(coords[i])
    if new_pos in occupied:
        return None
    new_coords = list(coords)
    new_coords[i] = new_pos
    return new_coords


def try_pivot_move(coords):
    # Rotate the tail after a random residue
    n = len(coords)
    if n < 3:
        return None
    i = random.randint(1, n - 2)
    pivot_x, pivot_y = coords[i]
    angle = random.choice([90, 180, 270])

    def rotate(px, py):
        dx, dy = px - pivot_x, py - pivot_y
        if angle == 90:
            dx, dy = -dy, dx
        elif angle == 180:
            dx, dy = -dx, -dy
        else:  # 270
            dx, dy = dy, -dx
        return pivot_x + dx, pivot_y + dy

    new_coords = coords[:i + 1] + [rotate(x, y) for x, y in coords[i + 1:]]
    if not is_self_avoiding(new_coords):
        return None
    return new_coords


MOVES = [try_end_move, try_corner_move, try_pivot_move]


def random_initial_conformation(n):
    # random self-avoiding start walking by backtrace 
    while True:
        coords = [(0, 0)]
        occupied = {(0, 0)}
        stuck = False
        for _ in range(n - 1):
            x, y = coords[-1]
            options = [(x + dx, y + dy) for dx, dy in DIRECTIONS]
            options = [o for o in options if o not in occupied]
            if not options:
                stuck = True
                break
            nxt = random.choice(options)
            coords.append(nxt)
            occupied.add(nxt)
        if not stuck:
            return coords


def anneal(sequence, steps=60000, t_start=3.0, t_end=0.02,
           cooling_rate=0.998, steps_per_temp=50, seed=None):
    # Return(best_energy,best_coords, final_coords).
    rng_state = random.getstate()
    if seed is not None:
        random.seed(seed)

    coords = random_initial_conformation(len(sequence))
    current_e = energy(sequence, coords)
    best_e, best_coords = current_e, list(coords)

    temperature = t_start
    step_count = 0
    while step_count < steps and temperature > t_end:
        for _ in range(steps_per_temp):
            move_fn = random.choice(MOVES)
            candidate = move_fn(coords)
            step_count += 1
            if candidate is None:
                continue  # skip
            cand_e = energy(sequence, candidate)
            delta = cand_e - current_e
            if delta <= 0 or random.random() < math.exp(-delta / temperature):
                coords = candidate
                current_e = cand_e
                if current_e < best_e:
                    best_e = current_e
                    best_coords = list(coords)
            if step_count >= steps:
                break
        temperature *= cooling_rate

    if seed is not None:
        random.setstate(rng_state)

    return best_e, best_coords, coords


def multi_restart(sequence, n_runs=100, **anneal_kwargs):

    #  Return (overall_best_energy, overall_best_coords,list_of_per_run_best_energies). over multiple annealing trajectories
    
    overall_best_e = None
    overall_best_coords = None
    per_run_energies = []

    for seed in range(n_runs):
        best_e, best_coords, _ = anneal(sequence, seed=seed, **anneal_kwargs)
        per_run_energies.append(best_e)
        if overall_best_e is None or best_e < overall_best_e:
            overall_best_e = best_e
            overall_best_coords = best_coords

    return overall_best_e, overall_best_coords, per_run_energies


if __name__ == "__main__":
    from collections import Counter

    seq = "HPHPPHHPHPPHPHHPPHPH"
    best_e, best_coords, per_run = multi_restart(seq, n_runs=100)
    print(f"sequence: {seq}")
    print(f"best energy found across 100 restarts: {best_e}")
    print(f"energy distribution across runs: {dict(sorted(Counter(per_run).items()))}")
    print(f"best conformation: {best_coords}")