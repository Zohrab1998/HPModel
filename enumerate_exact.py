# Exhaustive enumeration of self-avoiding walks on the 2D square lattice.

from hp_model import DIRECTIONS, energy


def enumerate_walks(length):

    # Yield every self-avoiding walk of the given length, starting fixed
    # at (0, 0), as a list of (x, y) coordinate tuples.
    
    start = (0, 0)

    def backtrack(coords):
        if len(coords) == length:
            yield list(coords)
            return
        last_x, last_y = coords[-1]
        for dx, dy in DIRECTIONS:
            nxt = (last_x + dx, last_y + dy)
            if nxt in coords:
                continue
            coords.append(nxt)
            yield from backtrack(coords)
            coords.pop()  # Undo the move and try the next one

    yield from backtrack([start])


def find_ground_state(sequence):

    # Exhaustively search all self-avoiding walks for the given sequence
    # and return (min_energy, list_of_all_coords_achieving_it).
    
    best_energy = None
    best_coords = []

    for coords in enumerate_walks(len(sequence)):
        e = energy(sequence, coords)
        if best_energy is None or e < best_energy:
            best_energy = e
            best_coords = [coords]
        elif e == best_energy:
            best_coords.append(coords)

    return best_energy, best_coords


if __name__ == "__main__":
    seq = "HPPHPH"
    best_e, best_confs = find_ground_state(seq)
    print(f"sequence: {seq}")
    print(f"ground-state energy: {best_e}")
    print(f"number of conformations achieving it (raw, incl. rotations/"
          f"reflections/direction-reversal): {len(best_confs)}")
    print(f"example conformation: {best_confs[0]}")