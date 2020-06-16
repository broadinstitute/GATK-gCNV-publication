import pyranges as pr
import pandas as pd

truth_callset_bed_file = "/Users/asmirnov/Desktop/evaluations/gCNV_common_improvements/1000G/1KGP_liftover/1KGP.hg19.cnv_only.bed"
output_callset_bed_file = "/Users/asmirnov/Desktop/evaluations/gCNV_common_improvements/1000G/1KGP_liftover/1KGP.hg19.cnv_only.preprocessed.bed"

bed_columns = ["#chr", "start", "end", "name", "svtype", "samples"]
pyrange_columns = ["Chromosome", "Start", "End", "Name", "Genotype", "Samples"]


class Interval:
    """Stores a Genomic interval"""

    def __init__(self, chrom: str, start: int, end: int):
        self.chrom = chrom
        assert (start <= end)
        self.start = start
        self.end = end

    def __eq__(self, other):
        return (self.chrom == other.chrom) and (self.start == other.start) and (self.end == other.end)

    def __hash__(self):
        return hash((self.chrom, self.start, self.end))


class Event:
    """Stores an event type and call qualities for a single interval and single sample"""

    def __init__(self, interval: Interval, name: str,  genotype: str):
        self.interval = interval
        self.name = name
        self.genotype = genotype

    def __eq__(self, other):
        return self.interval == other.interval \
                 and self.name == other.name and self.genotype == other.genotype

    def __hash__(self):
        return hash((self.interval, self.name, self.genotype))


truth_callset_df = pd.read_csv(open(truth_callset_bed_file, "r"), sep="\t", comment="#", names=pyrange_columns)
truth_callset_df = truth_callset_df.astype({'Chromosome': 'str'})
truth_callset_df["Samples"] = truth_callset_df["Samples"].map(lambda s: s.split(",")).map(frozenset)
samples = list(frozenset.union(*truth_callset_df['Samples'].values))
truth_callset_pyrange = pr.PyRanges(truth_callset_df)

sample_to_sample_to_event_list_map = {}
truth_callset_length = str(len(truth_callset_pyrange))
event_index = 0
for k, df in truth_callset_pyrange:
    for index, truth_event in df.iterrows():
        event_index += 1
        print("\rAnalyzing events: " + str(event_index) + "/" + str(truth_callset_length), end="", flush=True)
        for sample in truth_event["Samples"]:
            event = {"Chromosome": truth_event["Chromosome"], "Start": truth_event["Start"], "End": truth_event["End"],
                         "Genotype": truth_event["Genotype"], "Name": truth_event["Name"]}
            sample_to_sample_to_event_list_map.setdefault(sample, []).append(event)
print("")

updated_callset_event_to_sample_map = {}
num_events_removed = 0
for index, sample in enumerate(samples):
    print("\rProcessing samples: " + str(index+1) + "/" + str(len(samples)), end="", flush=True)
    events_df = pd.DataFrame(sample_to_sample_to_event_list_map[sample])

    sample_pr = pr.PyRanges(events_df)
    sample_pr = sample_pr.count_overlaps(sample_pr, overlap_col="Count")
    num_events_removed += len(sample_pr.df[sample_pr.df["Count"] > 1])
    sample_updated_df = sample_pr.df[sample_pr.df["Count"] == 1]

    for index, row in sample_updated_df.iterrows():
        event = Event(Interval(row["Chromosome"], row["Start"], row["End"]), row["Name"], row["Genotype"])
        updated_callset_event_to_sample_map.setdefault(event, []).append(sample)

updated_callset_dict_entries = []
for event, samples in updated_callset_event_to_sample_map.items():
    samples_str = ",".join(samples)
    contig, start, end, genotype, name = event.interval.chrom, event.interval.start, event.interval.end, event.genotype, event.name
    updated_callset_dict_entries.append({"#chr": contig, "start": start, "end": end,
                                         "svtype": truth_event["Genotype"], "name": name, "samples": samples_str})

updated_callset_df = pd.DataFrame(updated_callset_dict_entries, columns=bed_columns)
updated_callset_df.to_csv(output_callset_bed_file, sep='\t', header=True, index=False)

print("\nNumber of sample level events removed: " + str(num_events_removed))
print("Done!")

