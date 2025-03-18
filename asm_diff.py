#!/usr/bin/env python3
import os
import sys
import math
import json
import argparse
import subprocess
from tempfile import TemporaryDirectory


def SimpleFastaParser(handle):
    for line in handle:
        if line[0] == '>':
            title = line[1:].rstrip()
            break
    else:
        return
    lines = []
    for line in handle:
        if line[0] == '>':
            yield title, "".join(lines).replace(" ", "").replace('\r', "")
            lines = []
            title = line[1].rstrip()
            continue
        lines.append(line.rstrip())
    yield title, "".join(lines).replace(" ", "").replace('\r', "")


def assembly_accuracy(reference, query):
    with TemporaryDirectory() as tmpdir:
        prefix = os.path.join(tmpdir, 'nucmer_output')
        subprocess.run(['nucmer', '--maxmatch', '-p', prefix, reference, query])
        with open(prefix + '_filter.delta', 'w') as handle:
            child_process = subprocess.run(
                ['delta-filter', '-1', prefix + '.delta'],
                stdout=handle, text=True
            )
            child_process = subprocess.run(
                ['show-snps', '-rlTHC', prefix + '_filter.delta'],
                stdout=subprocess.PIPE, text=True
            )
            show_snps_output = child_process.stdout
            child_process = subprocess.run(
                ['show-coords', '-THrcl', prefix + '_filter.delta'],
                stdout=subprocess.PIPE, text=True
            )
            show_coords_output = child_process.stdout

    reference_bases = 0
    with open(reference) as handle:
        for seqid, seqchr in SimpleFastaParser(handle):
            reference_bases += len(seqchr)

    aligned_bases = 0
    for line in show_coords_output.splitlines():
        line = line.split()
        sstart, send = int(line[0]), int(line[1])
        aligned_bases += (send - sstart + 1)

    coverage = (aligned_bases / reference_bases) * 100

    mismatchs = 0
    indels = 0
    for line in show_snps_output.splitlines():
        line = line.split()
        ssub, qsub = line[1], line[2]
        if ssub == '.' or qsub == '.':
            indels += 1
        else:
            mismatchs += 1
    error_rate = (mismatchs + indels) / aligned_bases
    q_score = -10 * math.log10(error_rate) if error_rate > 0 else float('inf')
    return {
        'qscore': round(q_score, 2),
        'coverage': round(coverage, 2),
        'num_mismatches': mismatchs,
        'num_indels': indels
    }


def main():
    parser = argparse.ArgumentParser(
        description='Compare two assemblies.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("reference", help='Reference assembly.')
    parser.add_argument("query", help='Query assembly.')
    args = parser.parse_args()

    result = assembly_accuracy(args.reference, args.query)
    print(json.dumps(result), file=sys.stdout)


if __name__ == '__main__':
    main()