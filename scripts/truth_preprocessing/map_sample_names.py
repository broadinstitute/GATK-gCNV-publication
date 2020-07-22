import pandas as pd

callset_bed_file = "/Users/asmirnov/Desktop/SFARI_truth.sorted.bed"
mapping_file = "/Users/asmirnov/Desktop/SFARI_sample_id_map.tsv"
output_file = "/Users/asmirnov/Desktop/SFARI_truth.sorted.remapped_ids.bed"

bed_columns = ["#chr", "start", "end", "name", "svtype", "samples"]

callset_pd = pd.read_csv(open(callset_bed_file, "r"), sep="\t", comment="#", names=bed_columns)
mapping_table_pd = pd.read_csv(open(mapping_file, "r"), sep="\t", comment="#", names=["Person", "SSC_ID"])

def map_sample_list(samples):
    sample_list = samples.split(",")
    sample_output = []
    for ssd_id in sample_list:
        result_df = mapping_table_pd.loc[mapping_table_pd['SSC_ID'] == ssd_id]
        result_ids = result_df["Person"].tolist()
        if len(result_ids) == 0:
            sample_output.append(ssd_id)
            continue
        if len(result_ids) == 1:
            sample_output.append(result_df.iloc[0]["Person"])
        else:
            matching = [s for s in result_ids if "." in s]
            if matching:
                sample_output.append(matching[0])
            else:
                sample_output.append(ssd_id)
    return ",".join(sample_output)


callset_pd["samples"] = callset_pd["samples"].apply(map_sample_list)
callset_pd.to_csv(output_file, sep='\t', header=True, index=False)