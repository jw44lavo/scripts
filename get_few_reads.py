from Bio import SeqIO

number_of_reads = 10

input_file = "input.fastq"

output_file = "output.fastq"
output_handle = open(output_file, "w")

with open(input_file) as input_handle:
    for i, record in enumerate(SeqIO.parse(input_handle, "fastq")):
        
        if i <= number_of_reads:
            SeqIO.write(record, output_handle, "fastq")

        if i > number_of_reads:
            break
                
output_handle.close()


print("done")
