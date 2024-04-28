#!/usr/bin/env python3
''' compare_algs.py

    usage:

        ./compare_algs.py <csv_file>
'''


# ______________________________________________________________________
# Imports

import csv
import sys


# ______________________________________________________________________
# Functions

def print_row(row, field_lens):
    for i, field in enumerate(row):
        width = field_lens[i] + 2  # + 2 is just padding / buffer.
        s = field[:width].ljust(width, ' ')
        print(s, end=' ')
    print()

# ______________________________________________________________________
# Main

if __name__ == '__main__':

    if len(sys.argv) != 2:
        print(__doc__)
        sys.exit(0)

    with open(sys.argv[1]) as f:
        rows = list(csv.reader(f))

    alg_names = {row[0] for row in rows[1:]}
    num_cases = len([r for r in rows[1:] if r[0] == next(iter(alg_names))])

    field_lens = list(map(len, rows[0]))
    for i in range(num_cases):
        print()
        print_row(rows[0], field_lens)
        for j in range(len(alg_names)):
            print_row(rows[1 + i + j * num_cases], field_lens)
