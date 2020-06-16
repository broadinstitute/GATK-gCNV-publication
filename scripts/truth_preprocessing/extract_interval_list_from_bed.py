import pandas as pd
import os

file_directory = "/Users/asmirnov/Desktop/evaluations/gCNV_common_improvements/1000G/1KGP_liftover"
bed_hg38_path = os.path.join(file_directory, "1KGP.hg38.cnv_only.bed")
hg19_interval_list_header = os.path.join(file_directory, "hg38_interval_list_header.txt")

output_interval_list_path = os.path.join(file_directory, "1KGP.hg38.cnv_only.interval_list")

with open(hg19_interval_list_header, 'r') as header_file:
    header_string = header_file.read()

with open(output_interval_list_path, 'w') as output_file:
    output_file.write(header_string + "\n")
    with open(bed_hg38_path, 'r') as input_file:
        for line in input_file:
            if line.startswith("#"):
                continue
            split_line = line.split('\t')
            chrom = split_line[0]
            start = split_line[1]
            end = split_line[2]
            record_id = split_line[3]
            output_file.write("{}\t{}\t{}\t{}\t{}\n".format(chrom, start, end, "+", record_id))
