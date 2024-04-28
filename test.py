import os
import unittest
import random
import math

import editdistance  # third-party baseline

from minsh.astar import (
    h_dijkstra,
    align,
    build_seedh,
    build_seedh_for_pruning,
    print_stats,
)
from minsh.utils import (
    read_fasta_file,
    generate_random_sequence,
    apply_errors,
    save_fasta_file,
)


class TestFastaFunctions(unittest.TestCase):
    def setUp(self):
        self.sequence_length = 100
        self.description = "Test sequence"
        self.file_name = "test_sequence.fasta"
        self.random_sequence = generate_random_sequence(self.sequence_length)
        save_fasta_file(self.file_name, self.description, self.random_sequence)

    def tearDown(self):
        if os.path.exists(self.file_name):
            os.remove(self.file_name)

    def test_read_fasta_file(self):
        read_sequence = read_fasta_file(self.file_name)
        self.assertEqual(self.random_sequence, read_sequence)

    def test_save_fasta_file(self):
        with open(self.file_name, "r") as fasta_file:
            description_line = fasta_file.readline().strip()
            self.assertEqual(f"> {self.description}", description_line)

            sequence = ""
            for line in fasta_file:
                sequence += line.strip()
            self.assertEqual(self.random_sequence, sequence)


class TestAStar(unittest.TestCase):
    def setUp(self):
        random.seed(42)
        self.A = "ACCAGTGCCATT"
        self.B = "ACTAGTGGCACT"
        self.target = (len(self.A), len(self.B))

    def test_dijkstra(self):
        g, _, _ = align(self.A, self.B, h_dijkstra)
        self.assertEqual(g[self.target], editdistance.eval(self.A, self.B))

    def test_astar_with_seed_heuristic_small(self):
        k = 3
        h_seed = build_seedh(self.A, self.B, k)
        g, _, _ = align(self.A, self.B, h_seed)
        self.assertEqual(g[self.target], editdistance.eval(self.A, self.B))

    def test_astar_with_seedh_big(self):
        n = 10000
        A = "".join(random.choices("ACGT", k=n))
        B = apply_errors(A, 0.012)

        target = (len(A), len(B))
        k = math.ceil(math.log(len(A), 4))
        h_seed = build_seedh(A, B, k)
        g, _, _ = align(A, B, h_seed)
        # print_stats(A, B, k, g)
        # self.assertEqual(g[target], editdistance.eval(A, B))

    def test_astar_with_seedh_pruning(self):
        n = 100000
        A = "".join(random.choices("ACGT", k=n))
        B = apply_errors(A, 0.115)  # ~6% edit distance

        # k = 1/error_rate ?
        target = (len(A), len(B))
        k = math.ceil(math.log(len(A), 4))
        h_seed_prune = build_seedh_for_pruning(A, B, k)
        g_prune, _, _ = align(A, B, h_seed_prune)
        print_stats(A, B, k, g_prune)


#        self.assertEqual(g_prune[target], editdistance.eval(A, B))

if __name__ == "__main__":
    unittest.main()
