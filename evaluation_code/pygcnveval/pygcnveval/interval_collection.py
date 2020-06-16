import pandas as pd
import pyranges


class IntervalCollection:
    """Collection of genomic intervals"""

    INTERVAL_COLLECTION_COLUMNS = ["Chromosome", "Start", "End"]
    INTERVAL_COLLECTION_COLUMN_TYPES = {"Chromosome": "category", "Start": "int32", "End": "int32"}

    def __init__(self, pyrange: pyranges.PyRanges):
        self.pyrange = pyrange

    @classmethod
    def read_interval_list(cls, interval_list_file):
        intervals_df = pd.read_csv(open(interval_list_file, 'r'), comment="@", delimiter="\t", usecols=[0, 1, 2],
                                   header=None, names=IntervalCollection.INTERVAL_COLLECTION_COLUMNS)
        intervals_df = intervals_df.astype(IntervalCollection.INTERVAL_COLLECTION_COLUMN_TYPES)
        intervals_pyrange = pyranges.PyRanges(intervals_df)
        return cls(intervals_pyrange)

    def __len__(self):
        return len(self.pyrange.df)
