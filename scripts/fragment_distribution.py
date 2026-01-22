#!/usr/bin/env python3
import sys
import gzip
import matplotlib.pyplot as plt

def get_read_lengths(fastq_file, max_reads=None):
    # open gzip or plain text
    opener = gzip.open if fastq_file.endswith(".gz") else open
    lengths = []
    with opener(fastq_file, "rt") as f:
        for i, line in enumerate(f):
            if i % 4 == 1:  # sequence line
                lengths.append(len(line.strip()))
            if max_reads and len(lengths) >= max_reads:
                break
    return lengths

def main():
    if len(sys.argv) < 2:
        print(f"Usage: {sys.argv[0]} <fastq_file> [max_reads]")
        sys.exit(1)

    fastq_file = sys.argv[1]
    max_reads = int(sys.argv[2]) if len(sys.argv) > 2 else None

    lengths = get_read_lengths(fastq_file, max_reads=max_reads)

    if not lengths:
        print("No reads found.")
        sys.exit(1)

    average_length = sum(lengths) / len(lengths)
    print(f"Average read length: {average_length:.2f} bp (n={len(lengths)} reads)")

    plt.hist(lengths, bins=50, edgecolor="black")
    plt.xlabel("Read length (bp)")
    plt.ylabel("Count")
    plt.title("Read length distribution")
    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    main()
