import matplotlib.pyplot as plt

def draw_exploration(start, end, came_from):
    fig, ax = plt.subplots()
    
    # Draw grid lines
    for x in range(end[1] + 1):
        ax.axvline(x, color='k', lw=1)
    for y in range(end[0] + 1):
        ax.axhline(y, color='k', lw=1)

    # Draw exploration
    for node, parent in came_from.items():
        if parent is not None:
            ax.plot([node[1] + 0.5, parent[1] + 0.5], [node[0] + 0.5, parent[0] + 0.5], color='blue', lw=1, marker='o', markersize=4)

    # Draw start and end nodes
    ax.plot(start[1] + 0.5, start[0] + 0.5, marker='o', color='green', markersize=8, label='Start')
    ax.plot(end[1] + 0.5, end[0] + 0.5, marker='o', color='red', markersize=8, label='End')

    ax.legend()
    ax.axis('equal')

    plt.savefig('a_star_exploration.png', dpi=300)
    plt.show()
