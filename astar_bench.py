from enum import Enum, auto
from dataclasses import dataclass
from typing import Dict, List
from time import perf_counter
from collections import defaultdict
import math

import numpy as np
from stringzilla import File, Str

from astar import h_dijkstra, align, build_seedh, build_seedh_for_pruning


class AlgorithmType(Enum):
    WAGNER_FISCHER = auto()
    DIJKSTRA = auto()
    SEED = auto()
    SEED_PRUNING = auto()


@dataclass
class Result:
    matrix: np.ndarray
    distance: int
    comparisons: int


@dataclass
class BenchmarkResult:
    preprocessing_time: float
    run_time: float
    comparisons: int
    distance: int
    length_a: int
    length_b: int


def wagner_fisher(s1, s2) -> Result:
    # Create a matrix of size (len(s1)+1) x (len(s2)+1)
    matrix = np.zeros((len(s1) + 1, len(s2) + 1), dtype=int)

    # Initialize the first column and first row of the matrix
    for i in range(len(s1) + 1):
        matrix[i, 0] = i
    for j in range(len(s2) + 1):
        matrix[0, j] = j

    # Compute Levenshtein distance
    comparisons = 0
    for i in range(1, len(s1) + 1):
        for j in range(1, len(s2) + 1):
            substitution_cost = s1[i - 1] != s2[j - 1]
            matrix[i, j] = min(
                matrix[i - 1, j] + 1,  # Deletion
                matrix[i, j - 1] + 1,  # Insertion
                matrix[i - 1, j - 1] + substitution_cost,  # Substitution
            )
            comparisons += 1

    # Return the Levenshtein distance
    return Result(matrix, matrix[len(s1), len(s2)], comparisons)


def wrapped_dijkstra(A, B):
    return h_dijkstra


def wrapped_seed(A, B):
    k = math.ceil(math.log(len(A), 4))
    return build_seedh(A, B, k)


def wrapped_seed_prune(A, B):
    k = math.ceil(math.log(len(A), 4))
    return build_seedh_for_pruning(A, B, k)


if __name__ == "__main__":

    tokenization = "line"
    batch_size = 100

    datasets = [
        "data/enwik9.txt",
        "data/leipzig1M.txt",
        "data/xlsum.csv",
    ]

    for dataset in datasets:
        # Prepare a contain to assemble the results for the current dataset
        results_per_algo: Dict[AlgorithmType, List[BenchmarkResult]] = defaultdict(list)

        # Load the dataset and split it into whitespace or newline separated strings
        dataset_content = Str(File(dataset))
        strings = (
            dataset_content.splitlines()
            if tokenization == "line"
            else dataset_content.split()
        )

        # Random sample pairs from strings
        strings_a = strings.sample(batch_size)
        strings_b = strings.sample(batch_size)

        # Run the baseline algo, aggregating all the results for the Wagner Fisher
        for a, b in zip(strings_a, strings_b):
            start_time = perf_counter()
            result = wagner_fisher(a, b)
            end_time = perf_counter()
            results_per_algo[AlgorithmType.WAGNER_FISCHER].append(
                BenchmarkResult(
                    preprocessing_time=0,
                    run_time=end_time - start_time,
                    comparisons=result.comparisons,
                    distance=result.distance,
                    length_a=len(a),
                    length_b=len(b),
                )
            )

        for heursitic_generator, huristic_name in [
            (wrapped_dijkstra, AlgorithmType.DIJKSTRA),
            (wrapped_seed, AlgorithmType.SEED),
            (wrapped_seed_prune, AlgorithmType.SEED_PRUNING),
        ]:
            for a, b in zip(strings_a, strings_b):
                prep_time = perf_counter()
                heuristic = heursitic_generator(a, b)
                start_time = perf_counter()
                matrix, distance, comparisons = align(a, b, heuristic)
                end_time = perf_counter()
                results_per_algo[huristic_name].append(
                    BenchmarkResult(
                        preprocessing_time=0,  # TODO: Measure preprocessing time?
                        run_time=end_time - start_time,
                        comparisons=comparisons,
                        distance=distance,
                        length_a=len(a),
                        length_b=len(b),
                    )
                )

        # Print the results
        print(f"Dataset: {dataset}")
        print(results_per_algo)
