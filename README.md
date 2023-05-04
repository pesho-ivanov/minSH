# minSH

The simplest implementaion of A* with Seed heuristic (SH) for edit distance.
`astar.py` takes `k` and a file with two strings (`A` and `B`), and returns the
exact edit distance `ed(A,B)` between them in case it is `ed < |A|/k` or `Too different` otherwise. It splits `A` into seeds of length `k` and find a shortest path from `(0,0)` to `(|A|, |B|)` using A* with SH `sh(i,j) = |{ s | s.start >= i and s not in B }|`. That's it.
