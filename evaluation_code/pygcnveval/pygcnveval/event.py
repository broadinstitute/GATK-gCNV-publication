from enum import Enum
from interval import Interval
from typing import List, Set


class EventType(Enum):
    """Enumeration of possible alleles"""

    REF = 0
    DEL = 1
    DUP = 2
    NO_CALL = 3

    @classmethod
    def gcnv_call_to_event_type(cls, gcnv_call: str):
        if gcnv_call == '.':
            return EventType.NO_CALL
        return cls(int(gcnv_call))

    @staticmethod
    def get_event_type_from_svtype(sv_type: str):
        """
        This method will return None if Structural Variation event type is not a Copy Number Variant
        """
        cnv_type_events = {"DUP", "DEL"}
        if sv_type not in cnv_type_events:
            return None
        else:
            return EventType[sv_type]


class Event:
    """Stores an event type and call qualities for a single interval and single sample"""

    def __init__(self, interval: Interval, sample: str, event_type: EventType, call_attributes: map,
                 overlapping_target_set: set):
        self.interval = interval
        self.sample = sample
        self.event_type = event_type
        self.call_attributes = call_attributes
        self.overlapping_target_set = overlapping_target_set

    def __eq__(self, other):
        return self.interval == other.interval and self.sample == other.sample \
                 and self.event_type == other.event_type and self.call_attributes == other.call_attributes

    def compare_to(self, other, minimum_target_overlap: float) -> bool:
        """

        :param other: event to validate against
        :param minimum_target_overlap: minimum fraction of overlapping targets required to validate event
        :return: whether self validates against provided event given conditions
        """
        assert self.sample == other.sample
        if self.event_type != other.event_type:
            return False
        number_overlapping_targets = len(self.overlapping_target_set.intersection(other.overlapping_target_set))
        if number_overlapping_targets / len(self.overlapping_target_set) < minimum_target_overlap:
            return False
        return True

    def find_event_with_largest_overlap(self, event_list: List):
        """
        From list of events select one with largest overlap (not reciprocal) with interval in self
        :param event_list: given interval list
        :return: event with overlap
        """
        if not event_list:
            return None
        events_overlaps = [self.interval.get_overlap(event.interval) for event in event_list]
        return event_list[events_overlaps.index(max(events_overlaps))]

    def find_event_with_largest_reciprocal_overlap(self, event_list: List):
        """
        From list of events select one with largest reciprocal overlap with interval in self
        :param event_list: given interval list
        :return: event with largest reciprocal overlap
        """
        if not event_list:
            return None
        events_reciprocal_overlaps = [self.interval.get_reciprocal_overlap(event.interval) for event in event_list]
        return event_list[events_reciprocal_overlaps.index(max(events_reciprocal_overlaps))]


class Allele:
    """Stores an event type and call qualities for a single interval and single sample"""

    def __init__(self, interval: Interval, event_type: EventType, overlapping_target_set: Set[Interval],
                 sample_to_call_attributes_map: dict):
        self.interval = interval
        self.event_type = event_type
        self.overlapping_target_set = overlapping_target_set
        self.sample_to_call_attributes_map = sample_to_call_attributes_map

