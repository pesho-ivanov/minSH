import sys
from heapq import *
from math import sqrt

from utils import *

def euclidean_distance(a, b):
    return sqrt((a[0] - b[0])**2 + (a[1] - b[1])**2)

def next_states(A, B, curr):
    next_states = []
    if curr[0] < len(A): next_states.append((curr[0] + 1, curr[1]))
    if curr[1] < len(B): next_states.append((curr[0], curr[1] + 1))
    if curr[0] < len(A) and \
       curr[1] < len(B): next_states.append((curr[0] + 1, curr[1] + 1))
    return next_states

def edge_cost(A, B, curr, next):
    if curr[0]+1 == next[0] and curr[1]+1 == next[1] and A[curr[0]] == B[curr[1]]:
        return 0
    return 1

def align():
    Q = []                      # Priority queue with candidate states
    heappush(Q, (0, start))     # Push start state with priority 0
    prev = {start: None}        # Keep track of how we got to each state
    g = {start: 0}              # Cost of getting to each state

    while Q:
        _, curr = heappop(Q)

        if curr == target:
            break

        for next in next_states(A, B, curr):
            new_cost = g[curr] + edge_cost(A, B, curr, next)
            if next not in g or new_cost < g[next]:
                g[next] = new_cost
                priority = new_cost + euclidean_distance(target, next)
                heappush(Q, (priority, next))
                prev[next] = curr

    return prev, g

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python astar.py <a.fa> <b.fa> <k>")
        sys.exit(1)

    A, B = map(read_fasta_file, sys.argv[1,2])
    k = 14
    start = (0, 0)
    target = (len(A), len(B))

    #precompute_SH(A, B, k)
    came_from, cost_so_far = align()
    path = reconstruct_path(came_from, start, target)
    draw_exploration(start, target, came_from)

    print(f"Shortest path from {start} to {target} is:")
    print(" -> ".join(map(str, path)))
    print(f"Total cost: {cost_so_far[target]}")
