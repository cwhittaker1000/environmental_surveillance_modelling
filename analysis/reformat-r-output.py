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

# timestamp -> [reads1, reads2, ...]
data = defaultdict(list)

for line in sys.stdin:
    if not line.startswith("OUT"): continue
    line = line.strip()

    _, iteration, timestamp, reads, _ = line.split()

    iteration = int(iteration)
    timestamp = float(timestamp)
    reads = float(reads)

    single_kmer_reads = randomized_round(
        (READ_LENGTH-K)/(GENOME_LENGTH-READ_LENGTH)*reads)

    data[timestamp].append(single_kmer_reads)

# combine all samples for one day into a single sample
collated_data = {}
for timestamp, read_counts in data.items():
    timestamp = int(timestamp)
    if timestamp not in collated_data:
        collated_data[timestamp] = read_counts
    else:
        for i, rc in enumerate(read_counts):
            collated_data[timestamp][i] += rc

for timestamp, read_counts in sorted(collated_data.items()):
    print("%s\t%s" % (
        timestamp,
        "\t".join(str(x) for x in read_counts)))

    

    
