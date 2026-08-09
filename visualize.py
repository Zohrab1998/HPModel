# ASCII rendering of HP model
# residue is either (H or P)
# '-' Horizontal links
# '|' Vertical links
# Non-bonded H-H contacts are marked with 'o'


def render(sequence, coords):

    scaled = [(x * 2, y * 2) for x, y in coords]

    xs = [x for x, y in scaled]
    ys = [y for x, y in scaled]
    min_x, max_x = min(xs), max(xs)
    min_y, max_y = min(ys), max(ys)

    width = max_x - min_x + 1
    height = max_y - min_y + 1
    grid = [[' ' for _ in range(width)] for _ in range(height)]

    def put(x, y, ch, overwrite_space_only=False):
        row = max_y - y
        col = x - min_x
        if overwrite_space_only and grid[row][col] != ' ':
            return
        grid[row][col] = ch

    # Place residues
    for i, (x, y) in enumerate(scaled):
        put(x, y, sequence[i])

    # Place connectors between chain-adjacent residues
    for i in range(len(scaled) - 1):
        x1, y1 = scaled[i]
        x2, y2 = scaled[i + 1]
        mid_x, mid_y = (x1 + x2) // 2, (y1 + y2) // 2
        if x1 == x2:
            put(mid_x, mid_y, '|')
        else:
            put(mid_x, mid_y, '-')

    n = len(coords)
    for i in range(n):
        for j in range(i + 1, n):
            if j == i + 1:
                continue
            if sequence[i] != 'H' or sequence[j] != 'H':
                continue
            (x1, y1), (x2, y2) = coords[i], coords[j]
            if abs(x1 - x2) + abs(y1 - y2) != 1:
                continue
            mid_x = (scaled[i][0] + scaled[j][0]) // 2
            mid_y = (scaled[i][1] + scaled[j][1]) // 2
            put(mid_x, mid_y, 'o', overwrite_space_only=True)

    return '\n'.join(''.join(row) for row in grid)


if __name__ == "__main__":
    from enumerate_exact import find_ground_state

    seq = "HPPHPH"
    best_e, confs = find_ground_state(seq)
    print(f"sequence: {seq}  energy: {best_e}")
    print("(letters = residues, '-'/'|' = chain bonds, "
          "'o' = non-bonded H-H contact contributing to the energy)")
    print(render(seq, confs[0]))