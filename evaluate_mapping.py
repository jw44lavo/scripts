
from Bio import SeqIO
import pysam
from matplotlib import pyplot as plt
import numpy as np
import sys
import argparse
import os


def get_arguments():
    parser = argparse.ArgumentParser(usage="mapping evaluation")
    parser.add_argument("-i", "--input", metavar="", default=[], type = str, help="path to input file")
    parser.add_argument("-o", "--output", metavar="", default=[], type = str, help="path to output file")
    parser.add_argument("-f", "--format", metavar="", default=[], type = str, help="input file from 'mm2' or 'mmp'")
    parser.add_argument("-c", "--count", metavar="", default=[], type = str, help="count file from 'wc -l fastq'")
    parser.add_argument("-n", "--names", metavar="", default=[], type = str, help="id file with names")
    parser.add_argument("-p", "--plot", default=False, action="store_true", help="plot average nucleotide identity")

    return parser.parse_args()


def print_content(file):
    with open(file) as f:
        for line in f:
            print(line)


def get_read_count(file):
    with open(file) as f:
        i = int(int(f.readline().split(" ")[0])/4)
    return i


def get_ids(file):
    i = {}
    with open(file) as f:
        for record in f:
            i[record.split(" ")[0]] = record.split(" ", maxsplit=1)[1].strip("\n")
    return i


def plot_accuracies(acc, count, format, path):
    ylim = count/5
    plt.hist(acc, bins = 100)
    plt.xlabel('Average Nucleotide Identity/ %', fontsize=16)
    plt.ylabel('Frequency', fontsize=16)
    plt.xlim(0,100)
    plt.ylim(0,ylim)
    path = str(path) + "plot_" + str(format) + ".png"
    plt.savefig(path)
    #plt.show()
    result = "\n \n"
    result += "for accuracy plot see:\n"
    result += path
    return result


def do_stats(heading ,hits, reads, acc, l_que, l_tar, d_ids):

    percentage_mapping = "%.2f" % float((hits/reads)*100)
    average_nt_identity = "%.2f" % float(np.mean(acc))
    lowest_nt_identity = "%.2f" % float(min(acc))
    highest_nt_identity = "%.2f" % float(max(acc))
    s_que = set(l_que)
    count_multimaps = len(l_que) - len(s_que)
    max_multimaps = 0
    d_tar = {} #dictionary for targets(keys) and number of hits onto these targets(values)
    d_que = {} #dictionary for queries(keys) and number of targets of these queries(values)

### create dictionary of targets(d_tar) from list of targets(l_tar) ###
    for ele in l_tar:
        if ele not in d_tar:
            d_tar[ele] = 1
        elif ele in d_tar:
            d_tar[ele] += 1
    
### create dictionary of queries(d_que) from list of queries(l_que) ###
    for ele in l_que:
        if ele not in d_que:
            d_que[ele] = 1
        elif ele in d_que:
            d_que[ele] += 1

### get the highest number of targets of all queries
    for ele in d_que: 
        if d_que[ele] > max_multimaps:
            max_multimaps = d_que[ele]


    result = "\n                     "+ heading + " \n \n"
    result += "number of reads to be mapped:                   " + str(reads) + "\n"
    result += "number of mappings:                             " + str(hits) + "\n"
    result += "mapping percentage:                             " + str(percentage_mapping) + "%\n"
    result += "number of reads, mapped to several targets:     " + str(count_multimaps) + "\n"
    result += "highest number of targets for one read:         " + str(max_multimaps) + "\n"
    result += "average mapping nucleotide identity:            " + str(average_nt_identity) + "\n"
    result += "lowest mapping nucleotide identity:             " + str(lowest_nt_identity) + "\n"
    result += "highest mapping nucleotide identity:            " + str(highest_nt_identity) + "\n"
    result += "number of targets:                              " + str(len(d_tar)) + " out of " + str(len(d_ids)) + "\n \n \n"
    result += "target list with occurence: \n"
    
### sort dictionary of targets by descending number of hits onto the target ###
    d_tar = {k: v for k, v in sorted(d_tar.items(), key=lambda item: item[1], reverse = True)}

### get nice print out of the target dictionary with target names ###
    result += "{} \t \t {} \t {} \n".format("target", "#occurence", "target name")
    
    for item, amount in d_tar.items():
        if len(item) == 11:
            result += "{} \t {} \t \t {} \n".format(item, amount, d_ids.get(item))
        elif len(item) == 1:
            result += "{} \t \t {} \t \t {} \n".format(item, amount, d_ids.get(item))
        elif len(item) == 12:
            result += "{} \t {} \t \t {} \n".format(item, amount, d_ids.get(item))
    
    return result


def get_data_from_minimap2(file):
    name = "MINIMAP2"
    hits = 0
    l_tar = []
    l_acc = []
    l_que = []
    
    
    with open(file) as f:
        for line in f:
            if not line.startswith("@"):
                l_que.append(line.split()[0])
                l_tar.append(line.split()[2])
                hits += 1
    
    with open(file) as f:
        samfile = pysam.AlignmentFile(f, "r")
        for read in samfile:
            tags = dict(read.get_tags())
            if "NM" in tags.keys():
                l_acc.append((1-(tags["NM"]/read.query_alignment_length))*100)
        samfile.close()

    return name, hits, l_tar, l_acc, l_que


def get_data_from_mashmap(file):
    name = "MASHMAP"
    hits = 0
    l_tar = []
    l_acc = []
    l_que = []

    with open(file) as f:
        for record in f:
            l_acc.append(float(record.split(" ")[9].split("\n")[0]))
            l_tar.append(record.split(" ")[5])
            l_que.append(record.split(" ")[0])
            hits += 1
    
    return name, hits, l_tar, l_acc, l_que


def main():
    a = get_arguments()
    
    input_file = ""
    input_format = ""
    count_file = ""
    id_file = ""
    result = ""

### check for input ###
    if not a.input:
        raise Exception("no input file given (-i)")
    else:
        if os.path.exists(a.input):
            input_file = a.input
        else:
            raise Exception("given input does not exist")

### check for file format ###
    if not a.format:
        raise Exception("no input format given (-f)")
    else:
        if a.format == "mm2":
            input_format = "mm2"
        elif a.format == "mmp":
            input_format = "mmp"
        #else:
        #    raise Exception("input format not known")

### check for count file ###
    if not a.count:
        raise Exception("no count file given (-c)")
    else:
        if os.path.exists(a.count):
            count_file = a.count
        else:
            raise Exception("given count file does not exist")

### check for id file
    if not a.names:
        raise Exception("no id file given (-n)")
    else:
        if os.path.exists(a.names):
            id_file = a.names
        else:
            raise Exception("given id file does not exist")

### get data from  id- and count- file ###
    dict_ids = get_ids(id_file)
    count_reads = get_read_count(count_file)

### get data from mm2- or mmp- output ###
    if input_format == "mm2":
        name_heading, count_hits, list_targets, list_accuracies, list_queries = get_data_from_minimap2(input_file)
    elif input_format == "mmp":
        name_heading, count_hits, list_targets, list_accuracies, list_queries = get_data_from_mashmap(input_file)

### do the statistics ###
    result = do_stats(name_heading, count_hits, count_reads, list_accuracies, list_queries, list_targets, dict_ids)

### check for plot option ###
    if a.plot is True:
        plot_accuracies(list_accuracies, count_reads, input_format, "/home/johann/OneDrive/Projects/commonplace/scripts/")

### check for output ###
    if not a.output:
        print(result)
    else:
        with open(a.output, "w+") as o:
            o.write(result)

if __name__ == "__main__":
    
    try:
        main()

    except IOError:
        raise Exception("\n \n ups, something failed tremendously... good luck finding the problem \n")

