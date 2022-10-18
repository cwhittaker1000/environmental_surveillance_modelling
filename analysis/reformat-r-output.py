import sys
import random
from collections import defaultdict


GENOME_LENGTH=30000
K=40
READ_LENGTH=120  # 150bp in practice is more like 120

def randomized_round(x):
    whole = int(x)
    fractional = x - whole
    if random.random() < fractional:
        return whole + 1
    return whole

#timestamp -> [(reads1, S1), (reads2, S2), ...]
data = defaultdict(list)

for line in sys.stdin:
    if not line.startswith("OUT"): continue
    line = line.strip()

    _, iteration, timestamp, reads, S = line.split()

    iteration = int(iteration)
    timestamp = int(float(timestamp))
    reads = float(reads)
    S = int(float(S))

    single_kmer_reads = randomized_round(
        (READ_LENGTH-K)/(GENOME_LENGTH-READ_LENGTH)*reads)

    data[timestamp].append([single_kmer_reads, S])

for timestamp, read_counts in sorted(data.items()):
    print("%s\t%s" % (
        timestamp,
        "\t".join("%s:%s" % (reads, S) for (reads, S) in read_counts)))
    

    
