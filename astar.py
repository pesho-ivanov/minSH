import sys, math
import numpy as np      # To compute sum[i] = num[i] + sum[i+1]
from heapq import *     # Heap for the priority queue
from utils import *

def build_seed_heuristic(A, B, k):
    """Precomputes the seed heuristic for A and B with k-mers."""
    kmers = { B[i:i+k] for i in range(len(B)-k+1) }                 # O(nk), O(n) with rolling hash (Rabin-Karp)
    seeds = [ A[i:i+k] for i in range(0, len(A)-k+1, k) ]           # O(n)   
    is_seed_missing = [ s not in kmers for s in seeds ] + [False]   # O(n)
    suffix_sum = np.cumsum(is_seed_missing[::-1])[::-1]             # O(n)
    h_seed = lambda ij: suffix_sum[ceildiv(ij[0],k)]                # O(1)
    return h_seed

def next_states(curr, target):
    """Generates up to three states that follow curr (right, down, diagonal) but not after target."""
    next_states = []
    if curr[0] < target[0]: next_states.append((curr[0] + 1, curr[1]))      # Insertion
    if curr[1] < target[1]: next_states.append((curr[0], curr[1] + 1))      # Deletion
    if curr[0] < target[0] and \
       curr[1] < target[1]: next_states.append((curr[0] + 1, curr[1] + 1))  # Match or mismatch
    return next_states

def edit_cost(curr, next, A, B):
    """Computes the cost an edit: 0 for match, 1 otherwise."""
    return not (curr[0]+1 == next[0] and curr[1]+1 == next[1] and A[curr[0]] == B[curr[1]])

def align(A, B, h):
    """Standard A* on the grid A x B using a given heuristic h."""
    start = (0, 0)              # Start state
    target = (len(A), len(B))   # Target state
    Q = []                      # Priority queue with candidate states
    heappush(Q, (0, start))     # Push start state with priority 0
    g = { start: 0 }            # Cost of getting to each state

    while Q:
        _, curr = heappop(Q)

        if curr == target:
            return g

        for next in next_states(curr, target):
            new_cost_to_next = g[curr] + edit_cost(curr, next, A, B)
            if next not in g or new_cost_to_next < g[next]:
                g[next] = new_cost_to_next
                priority = new_cost_to_next + h(next)
                heappush(Q, (priority, next))

if __name__ == "__main__":
    if len(sys.argv) == 3:
        print("Usage: python astar.py <a.fa> <b.fa>")
    else:
        A, B = map(read_fasta_file, sys.argv[1,2])
        k = math.log(len(A), 4)
        h_seed = build_seed_heuristic(A, B, k)
        g, prev = align(A, B, h_seed)
        print("Edit distance: ", len(g[(len(A), len(B))]))
