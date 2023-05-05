# Minimalistic A* with Seed heuristic

This is a very small optimal global aligner (computes edit distance) that runs near-linearly for a limited error rate $\sim 1/k = 1/logn$. [`astar.py`](https://github.com/pesho-ivanov/minSeedHeuristic/blob/master/astar.py) (~50 Python loc) runs A* with the admissible _seed heuristic (SH)_, which is precomputed in linear time (todo: rolling hash) and queired in constant time:

```Python
def build_seed_heuristic(A, B, k):
    """Builds the admissible seed heuristic for A and B with k-mers."""
    kmers = { B[i:i+k] for i in range(len(B)-k+1) }                 # O(nk), O(n) with rolling hash (Rabin-Karp)
    seeds = [ A[i:i+k] for i in range(0, len(A)-k+1, k) ]           # O(n)   
    is_seed_missing = [ s not in kmers for s in seeds ] + [False]*2 # O(n)
    suffix_sum = np.cumsum(is_seed_missing[::-1])[::-1]             # O(n)
    h_seed = lambda ij: suffix_sum[ceildiv(ij[0],k)]                # O(1)
    return h_seed
```

Explaination:

`astar.py` takes `k` and a file with two strings (`A` and `B`), and returns the
exact edit distance `ed(A,B)` between them in case it is `ed < |A|/k` or `Too different` otherwise. It splits `A` into seeds of length `k` and find a shortest path from `(0,0)` to `(|A|, |B|)` using A* with SH `sh(i,j) = |{ s | s.start >= i and s not in B }|`. That's it. Or, in code:

## Run

Prerequisites:
```
pip install rolling
pip install numpy
pip install heapq
```

Compute edit distance between two sequences:
```
python astar.py A.fa B.fa
```

## Related work

[Detailed Publications on A* for alignment](https://pesho-ivanov.github.io/#A*%20for%20optimal%20sequence%20alignment)

[AStarix](https://github.com/eth-sri/astarix) -- Semi-global seq-to-graph aligner:
* [Ivanov et al., (RECOMB 2020)](https://link.springer.com/chapter/10.1007/978-3-030-45257-5_7) -- Introduces A* with a trie for semi-global.
* [Ivanov et al., (RECOMB 2022)](https://www.biorxiv.org/content/10.1101/2021.11.05.467453) -- Introduces SH for read alignment on general graph references using trie.

[A*PA](https://github.com/RagnarGrootKoerkamp/astar-pairwise-aligner) -- Global seq-to-seq aligner:
* [Groot Koerkamp and Ivanov (preprint 2023)](https://www.biorxiv.org/content/10.1101/2022.09.19.508631) Applies SH to global alignment (edit distance). Generalizes SH with chaining, inexact matches, gap costs (for higher error rates). Optimizes SH with pruning (for near-linear scaling with length), and A* with diagonal transition (for faster quadratic mode).
