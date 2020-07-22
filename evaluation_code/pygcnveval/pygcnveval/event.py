from enum import Enum
from interval import Interval
from typing import List


class EventType(Enum):
    """Enumeration of possible alleles"""

    REF = 0
    DEL = 1
    DUP = 2
    NO_CALL = 3

    @classmethod
    def gcnv_call_to_event_type(cls, gcnv_call: int):
        return cls(gcnv_call)

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

    def __init__(self, interval: Interval, sample: str, event_type: EventType, call_attributes: map):
        self.interval = interval
        self.sample = sample
        self.event_type = event_type
        self.call_attributes = call_attributes

    def __eq__(self, other):
        return self.interval == other.interval and self.sample == other.sample \
                 and self.event_type == other.event_type and self.call_attributes == other.call_attributes

    def compare_to(self, other, minimum_reciprocal_overlap: float) -> bool:
        """

        :param other: event to validate against
        :param minimum_reciprocal_overlap: minimum ro required to validate event
        :return: whether two events are equivalent given conditions
        """
        assert self.sample == other.sample
        if self.event_type != other.event_type:
            return False
        if self.interval.get_reciprocal_overlap(other.interval) < minimum_reciprocal_overlap:
            return False
        return True

    def find_event_with_largest_reciprocal_overlap(self, event_list: List):
        """
        From list of events select one with largest reciprocal overlap with interval in self
        :param event_list: given interval list
        :return: event with largest reciprocal overlap
        """
        if not event_list:
            return None
        events_reciprocal_ovelaps = [self.interval.get_reciprocal_overlap(event.interval) for event in event_list]
        return event_list[events_reciprocal_ovelaps.index(max(events_reciprocal_ovelaps))]
