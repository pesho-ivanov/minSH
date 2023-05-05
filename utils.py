import random
import matplotlib.pyplot as plt

def ceildiv(a, b):
    return -(a // -b)

def read_fasta_file(file_path):
    with open(file_path, "r") as fasta_file:
        fasta_file.readline()  # Skip the description line
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
    
def draw_exploration(target):
    start = (0, 0)
    fig, ax = plt.subplots()
    
    # Draw grid lines
    for x in range(target[1] + 1):
        ax.axvline(x, color='k', lw=1)
    for y in range(target[0] + 1):
        ax.axhline(y, color='k', lw=1)

    # Draw start and end nodes
    ax.plot(start[1] + 0.5, start[0] + 0.5, marker='o', color='green', markersize=8, label='Start')
    ax.plot(target[1] + 0.5, target[0] + 0.5, marker='o', color='red', markersize=8, label='End')

    ax.legend()
    ax.axis('equal')

    plt.savefig('a_star_exploration.png', dpi=300)
    plt.show()

def print_stats(A, B, g):
    target = (len(A), len(B))
    print(f"Aligning sequences with len(A)={len(A)} and len(B)={len(B)}:")
    #print(" -> ".join(map(str, path)))
    print(f"Edit distance: {g[target]}")
    print(f"Explored: {len(g)}")
