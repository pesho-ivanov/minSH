import random
import matplotlib.pyplot as plt

nucleotides = 'ACGT'

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
    return ''.join(random.choice(nucleotides) for _ in range(length))

def apply_errors(A, e):
    B = ''
    for i in range(len(A)):
        if random.random() > e:
            B += A[i]                               # No error
        elif random.random() < 1/3:
            B += random.choice(nucleotides)         # Substitution
        elif random.random() < 1/2:
            B += random.choice(nucleotides) + A[i]  # Insertion
        else:
            B += ''                                 # Deletion
    return B

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

def print_stats(A, B, k, g):
    target = (len(A), len(B))
    print(f"Aligning sequences with len(A)={len(A)}, k={k}:")
    print(f"Error rate: {g[target] / len(A) :.2f}")
    print(f"Explored band: {len(g) / len(A) :.2f}")
