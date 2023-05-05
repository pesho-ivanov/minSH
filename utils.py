import random
import matplotlib.pyplot as plt

def read_fasta_file(file_path):
    with open(file_path, "r") as fasta_file:
        # Skip the description line
        fasta_file.readline()
        
        # Read the sequence and remove newline characters
        sequence = ""
        for line in fasta_file:
            sequence += line.strip()
    
    return sequence

def generate_random_sequence(length):
    nucleotides = ['A', 'T', 'C', 'G']
    return ''.join(random.choice(nucleotides) for _ in range(length))

def save_fasta_file(file_name, description, sequence, line_length=50):
    with open(file_name, "w") as fasta_file:
        fasta_file.write(f"> {description}\n")
        for i in range(0, len(sequence), line_length):
            fasta_file.write(sequence[i:i + line_length] + "\n")
        
def reconstruct_path(came_from, end):
    path = []
    node = end
    while node != None:
        path.append(node)
        node = came_from[node]
    path.reverse()
    return path
    
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
