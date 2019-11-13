from Bio import SeqIO
import sys 

input_file = sys.argv[1]
output_file = sys.argv[2]
number_of_reads = int(sys.argv[3])
current_read_count = int(sys.argv[4]) -1
output_handle = open(output_file, "w")

fastq_formats = ["fastq", "fq"]
fasta_formats = ["fasta", "fa"]
file_format = ""

if input_file.split(".")[-1] in fastq_formats:
    file_format = "fastq"

elif input_file.split(".")[-1] in fasta_formats:
    file_format = "fasta"

with open(input_file) as input_handle:
    for i, record in enumerate(SeqIO.parse(input_handle, file_format)):
        
        if i > current_read_count and i <= current_read_count + number_of_reads:
            SeqIO.write(record, output_handle, file_format)

        elif i > current_read_count + number_of_reads:
            break
                
output_handle.close()

print("done")