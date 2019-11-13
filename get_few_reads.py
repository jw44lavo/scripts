from Bio import SeqIO
import sys 

input_file = sys.argv[1]
output_file = sys.argv[2]
number_of_reads = int(sys.argv[3])
current_read_count = int(sys.argv[4]) -1

output_handle = open(output_file, "w")

with open(input_file) as input_handle:
    for i, record in enumerate(SeqIO.parse(input_handle, "fastq")):
        
        if i > current_read_count and i <= current_read_count + number_of_reads:
            SeqIO.write(record, output_handle, "fastq")

        elif i > current_read_count + number_of_reads:
            break
                
output_handle.close()

print("done")
