# Write-up: HP Model Lattice Protein Folding

## Results

**HPPHPH (6-mer), exact enumeration:** energy -2, matching the value
given in the exercise. There are 8 conformations that reach it, but
they're all the same shape once you account for rotation, reflection,
and reading the chain backwards - so the ground state is unique.

**HPHPPHHPHPPHPHHPPHPH (20-mer), simulated annealing, 100 restarts:**
best energy -9, matching the published Unger-Moult value. 14 of 100
restarts reached it; none went lower.

```
6-mer:          20-mer:
P-H-P            P-H-P
| o |            | o |
HoH-P          P-HoH
                  | o
              P-HoH-P
              | o o |
            P-HoH-HoH
            | o | | |
            P-H-P P-P
```
`o` marks a non-bonded H-H contact. Counting them gives 2 and 9,
matching the reported energies.

## Implementation notes

I represent a conformation as a list of (x, y) integer coordinates,
one per residue. Coordinates made the two operations I needed most
checking self avoidance and checking spatial adjacency, trivial set
and Manhattan-distance checks, at the cost of a slightly more fiddly
pivot move later on.

The energy function loops over all residue pairs with
`itertools.combinations` (so nothing gets double-counted), skips pairs
that are adjacent in the chain, and counts a contact when both
residues are H and one lattice step apart. I checked it by hand
against a small fold before trusting it on anything real.

For the 6-mer I used exhaustive enumeration - recursive backtracking
over every self avoiding walk. Timing it at a few lengths (284 walks
at 6, 881,500 at 14) showed the count growing roughly 7x for every two
extra residues, which puts length 20 somewhere around 10^8-10^9 walks.
That's why the 20-mer needs a different approach.

For the 20-mer I used simulated annealing with three move types: 
end moves, corner moves and pivot moves. Metropolis acceptance, and
geometric cooling. Before trusting it on the 20-mer I ran it on the
6-mer: 30 out of 30 runs matched the exact answer.

## Discussion questions

**1. Why exclude chain-adjacent pairs?** 
They're always lattice-adjacent just from being connected in the chain, 
so counting them would add a fixed amount to every conformation's energy 
no matter how it's folded. Forgetting self avoidance is the more serious
mistake: it lets impossible, overlapping conformations into the
search, which quietly produces the wrong optimum with no error or
crash to flag it.

**2. Move set and ergodicity.** 
End and corner moves on their own are a known non-ergodic move set
certain configurations get "locked" and can't be escaped with just those two.
Pivot moves on their own are ergodic (any self avoiding walk can reach any other). 
I used all three: pivots to guarantee the search can in principle reach any
conformation, end/corner moves to do the cheaper local refinement. I
left out crankshaft moves, since they're common in the literature but
add nothing once pivots are already in the mix.

**3. Confidence in the 20-mer ground state.** 
For the 6-mer I have certainty every conformation was checked.
For the 20-mer I don't have that, so the confidence comes from 
a few weaker signals stacked together: 
the method already proved itself on the 6-mer, 14 of 100
independent restarts landed on the same value, and that value matches
what's published in the literature. That's solid evidence, but it's
not a proof the way enumeration is.

**4. Unique or degenerate ground state?** 
Confirmed unique (up to symmetry) for the 6-mer through explicit
symmetry reduction. I didn't run the same check on the 20-mer
it would mean clustering the distinct shapes that reach -9 rather than 
just counting how many restarts get there. Uniqueness matters because
a real protein has to fold reliably into one functional shape,
a sequence with several equally low-energy but differently shaped ground 
states behaves more like a disordered polymer than a protein.

**5. Sublattice parity.** 
Color the lattice like a checkerboard by `(x+y) mod 2`. Every step along 
the chain flips that color, so residue i always sits on color `i mod 2`. 
Two residues can only be adjacent on the lattice if they sit on opposite colors 
meaning i and j must differ in parity. 
I added this as a runtime check inside `energy()`, so any bug that produced
a same-parity "contact" would get caught immediately. 
It could also be used to skip half the pairs in the energy loop for a rough
2x speedup, though I didn't need that at this scale.

**6. Temperature sweep.** 
I ran fixed-temperature Metropolis MC on the 20-mer (no cooling) across a 
range of temperatures and averaged the energy at each:

| T | 0.1 | 0.3 | 0.5 | 0.75 | 1.0 | 1.5 | 2.0 | 3.0 | 5.0 |
|---|---|---|---|---|---|---|---|---|---|
| avg E | -6.2 | -6.2 | -5.3 | -2.7 | -2.1 | -1.5 | -1.4 | -1.2 | -0.9 |

The energy doesn't drift down smoothly as it cools there's a sharp
drop between T=0.75 and T=0.3. That's the signature of cooperativity:
the chain collapses from extended to compact over a narrow temperature
window rather than gradually picking up a few more contacts as it
cools, similar to how real proteins fold.

**7. What the model gets right and wrong.** 
It captures the actual driving force behind real folding
burying hydrophobic residues away from solvent while a 
connected chain limits which shapes are even possible 
and it reproduces real behavior like the cooperative
transition above and the idea that some sequences fold to a unique
shape more readily than others. What it misses is almost everything
else: 20 amino acid types collapsed into 2, a 3D continuous space
flattened to a rigid 2D grid, and no hydrogen bonding, side-chain
packing, electrostatics, or real folding kinetics.

**8. What I ruled out, and what I'd try next.** 
I considered adding crankshaft moves to the search but 
skipped them pivots already cover ergodicity, 
so crankshafts would just be extra code with nothing
new to show for it. I also considered trying to verify the 20-mer
exactly with some symmetry reduction to cut down the walk count, but
the growth rate I measured makes that hopeless symmetry only removes
a constant factor, not the exponential itself. Given another day, I'd
first cluster the 20-mer's -9 conformations by symmetry to settle the
uniqueness question there too, and then use the parity fact to speed
up the energy function so I could afford more and longer annealing runs.