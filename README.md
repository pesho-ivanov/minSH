# Minimalistic A* with Seed Heuristic

minSH aligns sequences A and B optimally, i.e. computes the exact edit distance between them. minSH is inspired by [minGPT](https://github.com/karpathy/minGPT) to be small, clean, interpretable and educational re-implementation of the recent aligning approach based on A* with _seed heuristic (SH)_. minSH runs near-linearly for a limited error rate ~ O(1/k) = O(1/logn). [`astar.py`](https://github.com/pesho-ivanov/minSeedHeuristic/blob/master/astar.py) (~50 loc) implements A* with seed heuristic $h_{seed}(i,j) = \Big| \big\\{ s \in Seeds_{\geq i} \mid  s \notin B \big\\} \Big|$:

```Python
def build_seed_heuristic(A, B, k):
    """Builds the admissible seed heuristic for A and B with k-mers."""
    
    kmers = { B[i:i+k] for i in range(len(B)-k+1) }                 # O(nk): Index all kmers from B. O(n) with rolling hash
    seeds = [ A[i:i+k] for i in range(0, len(A)-k+1, k) ]           # O(n): Split A into non-overlapping seeds of length k.
    is_seed_missing = [ s not in kmers for s in seeds ]             # O(n): Is seed unseen in B, so >=1 edit must be made to align it.
    suffix_sum = np.cumsum(is_seed_missing[::-1])[::-1]             # O(n): How many of the remaining seeds have to be edited
    
    return lambda ij: suffix_sum[ceildiv(ij[0],k)]                  # O(1): How many seeds starting after the current position i have to be edited?
```

Next, we just use the seed heuristic for a starndard A* on the alignment graph `A x B`:
```Python
h_seed = build_seed_heuristic(A, B, k=log(len(A)))
astar(A, B, h_seed)
```

## Usage

Prerequisites:
```
pip install rolling
pip install numpy
pip install heapq
```

`astar.py` takes `k` and a file with two strings (`A` and `B`), and returns the exact edit distance `ed(A,B)` between them:
```
python astar.py A.fa B.fa
```

## TODO

Optimizations:
* rolling hash: for linear time precomputation
* greedy matching (aka sliding)
* pruning, using index trees

Presentation:
* visualization of the alignment (png, gif)
* interactivity for adding patches
* report stats
* benchmark

## Related work

[Detailed Publications on A* for alignment](https://pesho-ivanov.github.io/#A*%20for%20optimal%20sequence%20alignment)

[AStarix](https://github.com/eth-sri/astarix) semi-global seq-to-graph aligner:
* [Ivanov et al., (RECOMB 2020)](https://link.springer.com/chapter/10.1007/978-3-030-45257-5_7) &mdash; Introduces A* with a trie for semi-global.
* [Ivanov et al., (RECOMB 2022)](https://www.biorxiv.org/content/10.1101/2021.11.05.467453) &mdash; Introduces SH for read alignment on general graph references using trie.

[A*PA](https://github.com/RagnarGrootKoerkamp/astar-pairwise-aligner) global seq-to-seq aligner:
* [Groot Koerkamp and Ivanov (preprint 2023)](https://www.biorxiv.org/content/10.1101/2022.09.19.508631) &mdash; Applies SH to global alignment (edit distance). Generalizes SH with chaining, inexact matches, gap costs (for higher error rates). Optimizes SH with pruning (for near-linear scaling with length), and A* with diagonal transition (for faster quadratic mode).
