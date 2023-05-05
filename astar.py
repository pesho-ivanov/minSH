import sys, math
import numpy as np      # To compute sum[i] = num[i] + sum[i+1]
from heapq import *     # Heap for the priority queue
from utils import *

#import rolling
#roll = rolling.PolynomialHash("ACCAGAT", 3)
#print(list(roll))

def build_seed_heuristic(A, B, k):
    """Builds the admissible seed heuristic for A and B with k-mers."""
    kmers = { B[i:i+k] for i in range(len(B)-k+1) }                 # O(nk), O(n) with rolling hash (Rabin-Karp)
    seeds = [ A[i:i+k] for i in range(0, len(A)-k+1, k) ]           # O(n)   
    is_seed_missing = [ s not in kmers for s in seeds ] + [False]*2 # O(n)
    suffix_sum = np.cumsum(is_seed_missing[::-1])[::-1]             # O(n)
    h_seed = lambda ij: suffix_sum[ceildiv(ij[0],k)]                # O(1)
    return h_seed

def next_states_with_cost(u, A, B):
    """Generates three states following curr (right, down, diagonal) with cost 0
    for match, 1 otherwise."""
    return [ ((u[0] + 1, u[1]    ), 1),
             ((u[0],     u[1] + 1), 1),
             ((u[0] + 1, u[1] + 1), A[u[0]] == B[u[1]]) ]

def align(A, B, h):
    """Standard A* on the grid A x B using a given heuristic h."""
    start = (0, 0)              # Start state
    target = (len(A), len(B))   # Target state
    Q = []                      # Priority queue with candidate states
    heappush(Q, (0, start))     # Push start state with priority 0
    g = { start: 0 }            # Cost of getting to each state
    A += '!'; B += '!'          # Barrier to avoid index out of bounds

    while Q:
        _, u = heappop(Q)                                   # Pop state u with lowest priority
        if u == target: return g                            # Return cost to target
        if u[0] > target[0] or u[1] > target[1]: continue   # Skip states after target
        
        for v, edit_cost in next_states_with_cost(u, A, B): # For all edges u->v
            new_cost_to_next = g[u] + edit_cost             # Try optimal path through u->v
            if v not in g or new_cost_to_next < g[v]:       # If new path is better 
                g[v] = new_cost_to_next                     # Update cost to v
                priority = new_cost_to_next + h(v)          # Compute priority
                heappush(Q, (priority, v))                  # Push v with new priority

if __name__ == "__main__":
    if len(sys.argv) == 3:
        print("Usage: python astar.py <a.fa> <b.fa>")
    else:
        A, B = map(read_fasta_file, sys.argv[1:3])
        k = math.log(len(A), 4)
        h_seed = build_seed_heuristic(A, B, k)
        g, prev = align(A, B, h_seed)
        print("Edit distance: ", len(g[(len(A), len(B))]))
