from itertools import combinations

# The 4 moves on a 2D square lattice: right, left, up, down.
DIRECTIONS = [(1, 0), (-1, 0), (0, 1), (0, -1)]


def is_self_avoiding(coords):
    # True if no 2 are in the same site
    return len(set(coords)) == len(coords)


def is_valid_chain(coords):
    # True if the chian is right and they are in the right order
    if not is_self_avoiding(coords):
        return False
    for (x1, y1), (x2, y2) in zip(coords, coords[1:]):
        if abs(x1 - x2) + abs(y1 - y2) != 1:
            return False
    return True


def energy(sequence, coords):
    # H-H contacts | they are not chain adjacent and their positions are distance 1 (Manhattan distance)
    n = len(sequence)
    assert len(coords) == n, "sequence and coords must be the same length"
 
    contacts = 0
    for i, j in combinations(range(n), 2):
        if j == i + 1:
            continue  # chain-adjacent: bonded neighbor, not a "contact"
        if sequence[i] != 'H' or sequence[j] != 'H':
            continue
        (x1, y1), (x2, y2) = coords[i], coords[j]
        dist = abs(x1 - x2) + abs(y1 - y2)
        if dist == 1:
            # adjacency implies opposite checkerboard color, i.e. (j - i) must be odd.
            assert (j - i) % 2 == 1, (
                f"parity violation: residues {i},{j} are lattice-adjacent "
                f"but (j-i)={j-i} is even - this should be impossible for "
                f"a valid self-avoiding walk"
            )
            contacts += 1
 
    return -contacts