# HP Model Lattice Protein Folding

Finds minimum-energy conformations of HP model sequences on the 2D
square lattice exact enumeration for short sequences, simulated
annealing for longer ones.

## Requirements

Python 3.8+, standard library only, no dependencies.

## Files

- `hp_model.py` - conformation representation, validity checks
  (self-avoidance, chain adjacency), and the energy function.
- `enumerate_exact.py` - exhaustive enumeration of self-avoiding
  walks, giving an exact ground state for short sequences.
- `anneal.py` - simulated annealing (end/corner/pivot moves,
  Metropolis acceptance, geometric cooling, multiple restarts) for
  sequences too long to enumerate.
- `temperature_sweep.py` - fixed-temperature Monte Carlo across a
  range of temperatures, used for the cooperativity question in the
  write-up.
- `visualize.py` - ASCII rendering of a conformation on the lattice.
- `run.py` - command-line entry point.
- `writeup.md` - results, design decisions, and discussion answers.

## Usage

Exact enumeration (only feasible up to roughly length 14-15):

```bash
python3 run.py --sequence HPPHPH --method exact
```

Simulated annealing with restarts (for longer sequences):

```bash
python3 run.py --sequence HPHPPHHPHPPHPHHPPHPH --method anneal --runs 100
```

`--runs` sets the number of independent restarts (default 100). More
restarts give better evidence for having found the true ground state,
at the cost of runtime - about 0.8s per restart on the 20-mer with the
default annealing settings.

Temperature sweep:

```bash
python3 temperature_sweep.py
```

## Output

Each run prints the energy found and an ASCII picture of the
conformation. Letters are residues, `-`/`|` are chain bonds, and `o`
marks a non-bonded H-H contact - the pairs that actually contribute to
the energy. Annealing runs also print the spread of energies across
restarts, which is the evidence for the ground-state claim.

## Validation

Both the energy function and the annealer are checked against exact
enumeration on the short sequence (HPPHPH, energy -2) before being
trusted on the 20-mer. Details and results are in `writeup.md`.