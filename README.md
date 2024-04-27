# Minimalistic Seed Heuristic for A*

_minSH_ is a short working implementation that aligns sequences $A$ and $B$ end to end using a minimal number of edit operations (substitutions, insertions and deletions). As a side effect, it computes the exact edit distance $ed(A,B)$ with near-linear scaling, given limited divergence $d=ed(A,B)/max(|A|, |B|)$. [`astar.py`](https://github.com/pesho-ivanov/minSeedHeuristic/blob/master/astar.py) (~50 loc) implements A* with seed heuristic $h_{seed}(i,j) = \Big| \big\\{ s \in Seeds_{\geq i} \mid  s \notin B \big\\} \Big|$ in a short and simple way:

```python
def build_seed_heuristic(A, B, k):
    """Builds the admissible seed heuristic for A and B with k-mers."""
    
    kmers = { B[j:j+k] for j in range(len(B)-k+1) }                 # O(nk): Index all kmers from B. O(n) with rolling hash
    seeds = [ A[i:i+k] for i in range(0, len(A)-k+1, k) ]           # O(n): Split A into non-overlapping seeds of length k.
    is_seed_missing = [ s not in kmers for s in seeds ]             # O(n): Is seed unseen in B, so >=1 edit must be made to align it.
    suffix_sum = np.cumsum(is_seed_missing[::-1])[::-1]             # O(n): How many of the remaining seeds have to be edited
    
    return lambda ij, k=k: suffix_sum[ ceildiv(ij[0], k) ]          # O(1): How many seeds starting after the current position i have to be edited?
```

Next, we just use the seed heuristic for a starndard A* on the alignment graph `A x B`:

```python
h_seed = build_seed_heuristic(A, B, k=log(len(A)))
astar(A, B, h_seed)
```

# Background

## Sequence alignment as a shortest path problem

We can look at aligning sequences A and B as a process of sequentially aligning longer and longer prefixes of $A$ ($A[\dots i]$) to prefixes of $B$ ($B[\dots j]$) by matching, substituting, inserting or deleting single letters (minimizing [edit distance](https://en.wikipedia.org/wiki/Edit_distance)). Thus, finding an alignment with a minimal number of edit operations is equivalent to finding a shortest path starting from node $s=(0,0)$ and ending at node $t=(|A|, |B|)$ in a graph of prefixes (also called _edit graph_ or _alignment graph_), where each edge corresponds to one operation (with cost $0$ for a match or $1$ otherwise). This graph representation enables us to apply general shortest path algorithms.

## Dijkstra's shortest path algorithm

The simplest algorithm we can use is Dijkstra's algorithm which finds a shortest path of length $d$ by sequentially exploring nodes $u$ by increasing distance $dist(s,u)$ from the start node $s$ and until expanding the end node $t$. The problem is that the search circles around $s$ regardless of where the target $t$ is, so in our 2D lattice graph the number of explored nodes with $dist(s,u) < d$ grows quadratically in $d$. For most data (e.g. genetic) the edit distance $d$ grows proportionally to $|s|$ and $|t|$, so the whole algorithm becomes quadratic which is practically infeasible for long sequences.

## A* algorithm

The A* algirthm is a generalization of Dijkstra's algorithm that explores the nodes $u$ not just by their distance from the start $dist(s,u)$ but also adding an estimation of the remaining distance to the target $dist(s,u) + h(u,t)$. This heuristic function $h(u,t)$ allows for a potentially very direct search towards the target but it has to be designed depending on specific knowledge of the graph/task to be:
1. [admissible](https://en.wikipedia.org/wiki/Admissible_heuristic) (i.e. to never exceed the remaining distance $d(u,t)$), or otherwise the found path may not be shortest.
2. accurate in estimating $dist(s,u)$, or otherwise the search will not be directly going to $t$
3. fast to be computed for each explored node, or otherwise, the A* algorithm will be slow in practice

## Usage

Prerequisites:

```bash
pip install rolling
pip install numpy
pip install heapq
pip install fenwick
```

Run tests first:`

```bash
python test_astar.py
```

`astar.py` takes `k` and a file with two strings (`A` and `B`), and returns the exact edit distance `ed(A,B)` between them:

```bash
python astar.py data/small_A.fa data/small_B.fa
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

_minSH_ is inspired by [minGPT](https://github.com/karpathy/minGPT) to be small, clean, interpretable and educational re-implementation of the recent aligning approach based on the [A* shortest path algorithm](https://en.wikipedia.org/wiki/A*_search_algorithm).

[Detailed Publications on A* for alignment](https://pesho-ivanov.github.io/#A*%20for%20optimal%20sequence%20alignment)

[AStarix](https://github.com/eth-sri/astarix) semi-global seq-to-graph aligner:

* [Ivanov et al., (RECOMB 2020)](https://link.springer.com/chapter/10.1007/978-3-030-45257-5_7) &mdash; Introduces A* with a trie for semi-global.
* [Ivanov et al., (RECOMB 2022)](https://www.biorxiv.org/content/10.1101/2021.11.05.467453) &mdash; Introduces SH for read alignment on general graph references using trie.

[A*PA](https://github.com/RagnarGrootKoerkamp/astar-pairwise-aligner) global seq-to-seq aligner:

* [Groot Koerkamp and Ivanov (preprint 2023)](https://www.biorxiv.org/content/10.1101/2022.09.19.508631) &mdash; Applies SH to global alignment (edit distance). Generalizes SH with chaining, inexact matches, gap costs (for higher error rates). Optimizes SH with pruning (for near-linear scaling with length), and A* with diagonal transition (for faster quadratic mode).
Licensed under the Mozilla Public License, Version 2.0. In short, you are free to use and abuse, but give it back to the community.
