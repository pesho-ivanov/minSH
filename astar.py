import sys
import heapq
from math import sqrt

def euclidean_distance(a, b):
    return sqrt((a[0] - b[0])**2 + (a[1] - b[1])**2)

def astar(graph, start, end):
    frontier = []
    heapq.heappush(frontier, (0, start))
    came_from = {start: None}
    cost_so_far = {start: 0}

    while frontier:
        _, current = heapq.heappop(frontier)

        if current == end:
            break

        for next in graph[current]:
            new_cost = cost_so_far[current] + graph[current][next]
            if next not in cost_so_far or new_cost < cost_so_far[next]:
                cost_so_far[next] = new_cost
                priority = new_cost + euclidean_distance(end, next)
                heapq.heappush(frontier, (priority, next))
                came_from[next] = current

    return came_from, cost_so_far

def reconstruct_path(came_from, start, end):
    path = [end]
    node = end
    while node != start:
        node = came_from[node]
        path.append(node)
    path.reverse()
    return path

def main():
    if len(sys.argv) != 3:
        print("Usage: python a_star.py <start_node> <end_node>")
        sys.exit(1)

    start = tuple(map(int, sys.argv[1].split(',')))
    end = tuple(map(int, sys.argv[2].split(',')))

    # Example grid-based graph using a dictionary
    graph = {
        (0, 0): {(0, 1): 1, (1, 0): 1},
        (0, 1): {(0, 0): 1, (0, 2): 1},
        (0, 2): {(0, 1): 1, (1, 2): 1},
        (1, 0): {(1, 1): 1, (0, 0): 1},
        (1, 1): {(1, 0): 1, (1, 2): 1, (2, 1): 1},
        (1, 2): {(1, 1): 1, (0, 2): 1, (2, 2): 1},
        (2, 1): {(2, 0): 1, (1, 1): 1, (2, 2): 1},
        (2, 0): {(2, 1): 1},
        (2, 2): {(1, 2): 1, (2, 1): 1}
    }

    came_from, cost_so_far = astar(graph, start, end)
    path = reconstruct_path(came_from, start, end)

    print(f"Shortest path from {start} to {end} is:")
    print(" -> ".join(map(str, path)))
    print(f"Total cost: {cost_so_far[end]}")

if __name__ == "__main__":
    main()


