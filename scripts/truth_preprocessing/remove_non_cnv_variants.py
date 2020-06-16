import os

file_directory = "/Users/asmirnov/Desktop/evaluations/gCNV_common_improvements/1000G/1KGP_liftover"
original_vcf_hg38_path = os.path.join(file_directory, "1KGP.hg38.vcf")
output_cnv_vcf_hg38_path = os.path.join(file_directory, "1KGP.hg38.cnv_only.vcf")

comments = []
comments_written = False
with open(original_vcf_hg38_path, 'r') as vcf_reader:
    with open(output_cnv_vcf_hg38_path, 'w') as vcf_writer:
        for line in vcf_reader:
            if line.startswith("#"):
                comments.append(line)
                continue
            elif not comments_written:
                for comment in comments:
                    vcf_writer.write(comment)
                comments_written = True
            split_line = line.split('\t')
            if "<DEL>" in split_line[4] or "<DUP>" in split_line[4]:
                vcf_writer.write(line)