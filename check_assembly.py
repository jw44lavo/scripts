from Bio import SeqIO
import numpy as np
import sys
import argparse
import os

def get_arguments():
    parser = argparse.ArgumentParser(usage="script to filter assembled contigs by length")
    parser.add_argument("-i", "--input",        metavar="", default=[],     type = str, help= "path to input file")
    parser.add_argument("-o", "--output",       metavar="", default=[],     type = str, help= "path to output file for contigs >= length ")
    parser.add_argument("-l", "--length",       metavar="", default=2000,   type = int, help= "set filter length; default = 2000")
    parser.add_argument("-a", "--alternative",  metavar="", default=[],     type = str, help= "path to output file for contigs < length")

    return parser.parse_args()



def main():

    a = get_arguments()

    input_file = ""
    output_file = ""
    filter_length = 0

    input_file = a.input
    output_file = a.output
    filter_length = a.length
    alternative_output_file = a.alternative

    with open(input_file) as f:
        i = 0
        k = 0
        for line in f:
            if line.startswith(">"):
                length = int(line.split(" ")[3].split("=")[1])

                if length >= filter_length:
                    i += 1 
                    if i == 1:
                        with open(output_file, "w") as o:
                            o.write(line + next(f))
                    
                    elif i > 1:
                        with open(output_file, "a") as o:
                            o.write(line + next(f))

                if length < filter_length:
                    k += 1
                    if alternative_output_file:
                        if k == 1:
                            with open(alternative_output_file, "w") as o:
                                o.write(line + next(f))
                            
                        elif k > 1:
                            with open(alternative_output_file, "a") as o:
                                o.write(line + next(f))
        if i == 0 and output_file:
            with open(output_file, "w") as o:
                o.write("")

        if k == 0 and alternative_output_file:
            with open(alternative_output_file, "w") as o:
                o.write("")



if __name__ == "__main__":
    
    try:
        main()

    except IOError:
        raise Exception("\n \n ups, something failed tremendously... good luck finding the problem \n")

