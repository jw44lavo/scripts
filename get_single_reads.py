from Bio import SeqIO
import sys 

input_file = sys.argv[1]
input_base = input_file.split("/")[-1]
input_format = input_base.split(".")[-1]
input_name = input_base.split(".")[0]


fastq_formats = ["fastq", "fq"]
fasta_formats = ["fasta", "fa"]

### decision about input file format ###
if input_format in fastq_formats:
    input_format = "fastq"

elif input_format in fasta_formats:
    input_format = "fasta"


### read input and write output
with open(input_file) as input_handle:
    for i,record in enumerate(SeqIO.parse(input_handle, input_format)):
        output_file = input_name + "_" + str(i) + "." + input_format
        output_handle = open(output_file, "w")
        SeqIO.write(record, output_handle, input_format)

output_handle.close()

