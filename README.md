# minSeedHeuristic

The shortest and simplest implementaion of A* with Seed heuristic (SH) for edit distance. For pedagogical purposes.

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

## Explaination

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

## Cite

[Publications on A* for optimal alignment](https://pesho-ivanov.github.io/#A*%20for%20optimal%20sequence%20alignment)

## TODO
* rolling hash for O(n) preprocessing
* stats and visualization
