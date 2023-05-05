import sys
from heapq import *     # Heap for the priority queue
import numpy as np      # To compute sum[i] = num[i] + sum[i+1]
from utils import *

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
    if curr[0]+1 == next[0] and curr[1]+1 == next[1] and A[curr[0]] == B[curr[1]]:
        return 0    # Match
    return 1        # Mismatch or gap

def align(A, B, h):
    start = (0, 0)              # Start state
    target = (len(A), len(B))   # Target state
    
    Q = []                      # Priority queue with candidate states
    heappush(Q, (0, start))     # Push start state with priority 0
    prev = {start: None}        # Keep track of how we got to each state
    g = {start: 0}              # Cost of getting to each state

    while Q:
        _, curr = heappop(Q)

        if curr == target:
            break

        for next in next_states(curr, target):
            new_cost_to_next = g[curr] + edit_cost(curr, next, A, B)
            if next not in g or new_cost_to_next < g[next]:
                g[next] = new_cost_to_next
                priority = new_cost_to_next + h(next)
                heappush(Q, (priority, next))
                prev[next] = curr

    return g, prev

def h_dijkstra(u):
    return 0

def build_seed_heuristic(A, B, k):
    """Precomputes the seed heuristic for A and B with k-mers."""

    kmers = { B[i:i+k] for i in range(len(B)-k+1) }          # O(nk), O(n) with rolling hash (Rabin-Karp)
    print('B:', B)
    print('kmers:', kmers)
    
    seeds = [ A[i:i+k] for i in range(0, len(A)-k+1, k) ]    # O(n)
    print('A:', A)
    print('seeds:', seeds)
    
    is_seed_missing = [ s not in kmers for s in seeds ] + [False]      # O(n)
    print('is_seed_missing:', is_seed_missing)
    
    suffix_sum = np.cumsum(is_seed_missing[::-1])[::-1]      # O(n)
    print('suffix_sum:', suffix_sum)

    h_seed = lambda ij: suffix_sum[ceildiv(ij[0],k)]         # O(1)
    return h_seed

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python astar.py <a.fa> <b.fa> <k>")
        sys.exit(1)

    A, B = map(read_fasta_file, sys.argv[1,2])
    k = 14
    target = (len(A), len(B))

    h_seed = build_seed_heuristic(A, B, k)
    g, prev = align(A, B, h_seed)
    print_stats(A, B, g, prev)
