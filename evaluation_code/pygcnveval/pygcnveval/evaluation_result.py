
class EvaluationResult:
    DEFAULT_EVENT_SIZE_BINS = list(range(25))

    def __init__(self, event_size_bins=DEFAULT_EVENT_SIZE_BINS):
        self.event_size_bins = event_size_bins
        self.precision_size_to_tp = {i: 0 for i in event_size_bins}
        self.precision_size_to_fp = {i: 0 for i in event_size_bins}

        self.recall_size_to_tp = {i: 0 for i in event_size_bins}
        self.recall_size_to_fn = {i: 0 for i in event_size_bins}

    def update_precision(self, num_bins: int, event_valid: bool):
        if event_valid:
            self.precision_size_to_tp[num_bins - 1] += 1
        else:
            self.precision_size_to_fp[num_bins - 1] += 1

    def update_recall(self, num_bins: int, event_valid: bool):
        if event_valid:
            self.recall_size_to_tp[num_bins - 1] += 1
        else:
            self.recall_size_to_fn[num_bins - 1] += 1

    def __eq__(self, other):
        if isinstance(other, EvaluationResult):
            return self.precision_size_to_tp == other.precision_size_to_tp\
                   and self.precision_size_to_fp == other.precision_size_to_fp\
                   and self.recall_size_to_tp == other.recall_size_to_tp\
                   and self.recall_size_to_fn == other.recall_size_to_fn
        return False
