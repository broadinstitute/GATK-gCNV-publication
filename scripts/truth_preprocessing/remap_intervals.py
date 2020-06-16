import os
import pandas as pd

file_directory = "/Users/asmirnov/Desktop/evaluations/gCNV_common_improvements/1000G/1KGP_liftover"
original_bed_hg38_path = os.path.join(file_directory, "1KGP.hg38.cnv_only.bed")
output_bed_hg19_path = os.path.join(file_directory, "1KGP.hg19.cnv_only.bed")
remapped_intervals_hg19_path = os.path.join(file_directory, "1KGP.hg19.cnv_only.interval_list")


original_bed_hg38_df = pd.read_csv(original_bed_hg38_path,
                                   names=["chr", "start", "end", "variant_id", "svtype", "samples"],
                                   comment="#", delimiter='\t',
                                   dtype={"chr": object, "start": int, "end": int, "variant_id": object, "svtype": object, "samples": object})
original_bed_hg38_df = original_bed_hg38_df.set_index('variant_id')

remapped_intervals_hg19_df = pd.read_csv(remapped_intervals_hg19_path,
                                         names=["chr", "start", "end", "+", "variant_id"],
                                         comment="@", delimiter='\t',
                                         dtype={"chr": object, "start": int, "end": int, "+": object, "variant_id": object})

output_bed_hg19_df = pd.DataFrame(columns=["chr", "start", "end", "name", "svtype", "samples"])


remapped_intervals_hg19_df["svtype"] = original_bed_hg38_df.loc[remapped_intervals_hg19_df["variant_id"]]["svtype"].values
remapped_intervals_hg19_df["samples"] = original_bed_hg38_df.loc[remapped_intervals_hg19_df["variant_id"]]["samples"].values
remapped_intervals_hg19_df = remapped_intervals_hg19_df.drop(columns=["+"])

remapped_intervals_hg19_df.to_csv(output_bed_hg19_path, sep='\t', index=False, header=["#chr", "start", "end", "name", "svtype", "samples"])