
from Bio import SeqIO
import sys
from matplotlib import pyplot as plt
import numpy as np

input_file = sys.argv[1]
reads_file = sys.argv[2]
count_reads = int(sum(1 for line in open(reads_file))/4)
id_file = sys.argv[3]
set_names = {}

with open(id_file) as input_handle:
    for record in input_handle:
        set_names[record.split(" ")[0]] = record.split(" ", maxsplit=1)[1]

list_targets = []
list_accuracies = []
list_queries = []
count_hits = 0

with open(input_file) as input_handle:
    for record in input_handle:
        list_accuracies.append(float(record.split(" ")[9].split("\n")[0]))
        list_targets.append(record.split(" ")[5])
        list_queries.append(record.split(" ")[0])
        count_hits += 1

percentage_mapping = "%.2f" % float((count_hits/count_reads)*100)
set_queries = set(list_queries)
count_multimaps = len(list_queries) - len(set_queries)
set_targets = {}#set(list_targets)

for ele in list_targets:
    if ele not in set_targets:
        set_targets[ele] = 1

    elif ele in set_targets:
        set_targets[ele] += 1


result = "\n"
result += "number of reads to be mapped:                " + str(count_reads) + "\n \n"
result += "number of mappings:                          " + str(count_hits) + "\n \n"
result += "mapping percentage:                          " + percentage_mapping + "% \n \n"
result += "number of multimappings:                     " + str(count_multimaps) + "\n \n"
result += "average mapping nucleotide identity:         " + str(np.mean(list_accuracies)) + "\n \n"
result += "lowest mapping nucleotide identity:          " + str(min(list_accuracies)) + "\n \n"
result += "highest mapping nucleotide identity:         " + str(max(list_accuracies)) + "\n \n"
result += "number of targets:                           " + str(len(set_targets)) + "\n \n"
result += "target list with occurence: \n               " 

print(result)

print("target","       ", "#occurence", "   ", "target species")
for key in set_targets:
    print(key,"      ", set_targets[key], "        ", set_names.get(key))

