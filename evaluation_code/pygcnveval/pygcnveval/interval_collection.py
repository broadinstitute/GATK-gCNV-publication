import pandas as pd
import pyranges


class IntervalCollection:
    """Collection of genomic intervals"""
    interval_collection_columns = ["Chromosome", "Start", "End"]
    interval_collection_column_types = {"Chromosome": "category", "Start": "int32", "End": "int32"}

    def __init__(self, pyrange: pyranges.PyRanges):
        self.pyrange = pyrange

    @classmethod
    def read_interval_list(cls, interval_list_file):
        intervals_df = pd.read_csv(open(interval_list_file, 'r'), comment="@", delimiter="\t", usecols=[0, 1, 2],
                                   header=None, names=IntervalCollection.interval_collection_columns)
        intervals_df = intervals_df.astype(IntervalCollection.interval_collection_column_types)
        intervals_pyrange = pyranges.PyRanges(intervals_df)
        return cls(intervals_pyrange)
