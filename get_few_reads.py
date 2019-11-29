from Bio import SeqIO
import sys 

input_file = sys.argv[1]
input_base = input_file.split("/")[-1]
input_format = input_base.split(".")[-1]
input_name = input_base.split(".")[0]

number_of_reads = int(sys.argv[2])          # number of reads to select
current_read_count = int(sys.argv[3]) -1    # counter, where to start to select

fastq_formats = ["fastq", "fq"]
fasta_formats = ["fasta", "fa"]


### decision about input file format ###
if input_format in fastq_formats:
    input_format = "fastq"

elif input_format in fasta_formats:
    input_format = "fasta"


output_file = input_name + "_" + str(number_of_reads) + "_" + str(current_read_count + 1) + "." + input_format


output_handle = open(output_file, "w")

### read input and write output
with open(input_file) as input_handle:
    for i, record in enumerate(SeqIO.parse(input_handle, input_format)):
        
        if i > current_read_count and i <= current_read_count + number_of_reads:
            SeqIO.write(record, output_handle, input_format)

        elif i > current_read_count + number_of_reads:
            break
                
output_handle.close()
