import numpy as np


class Interval:
    """Stores a Genomic interval"""

    def __init__(self, chrom: str, start: int, end: int):
        self.chrom = chrom
        assert (start <= end)
        self.start = start
        self.end = end
        self.length = end - start + 1

    def intersects_with(self, interval):
        if self.chrom != interval.chrom:
            return False
        return ((self.start >= interval.start) and (self.start <= interval.end)) or \
               ((interval.start >= self.start) and (interval.start <= self.end))

    def get_overlap(self, other):
        """
        Calculate overlap (non-reciprocal) given interval with self, i.e. what fraction of self is covered by given
        interval
        :param other: interval to calculate overlap with
        :return: overlap fraction
        """
        if not self.intersects_with(other):
            return 0.
        overlap = np.minimum(self.end, other.end) - np.maximum(self.start, other.start)
        return overlap / self.length

    def get_reciprocal_overlap(self, other):
        """
        Calculate reciprocal overlap of two intervals
        :param other: interval to calculate overlap with
        :return: reciprocal overlap
        """
        if not self.intersects_with(other):
            return 0.
        return Interval.calculate_reciprocal_overlap(self.start, self.end, other.start, other.end)

    @staticmethod
    def calculate_reciprocal_overlap(start_0, end_0, start_1, end_1):
        overlap = np.minimum(end_0, end_1) - np.maximum(start_0, start_1)
        return np.minimum(overlap / (end_0 - start_0 + 1), overlap / (end_1 - start_1 + 1))

    def __eq__(self, other):
        return (self.chrom == other.chrom) and (self.start == other.start) and (self.end == other.end)

    def __str__(self):
        return str(self.chrom) + "_" + str(self.start) + "_" + str(self.end)

    def __hash__(self):
        return hash((self.chrom, self.start, self.end))
