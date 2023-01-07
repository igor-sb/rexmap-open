#include "MergedAlignment.h"
#include "QualityMergingMap.h"

MergedAlignment merge_paired_sequence_and_quality(
    PairedString &paired_sequence,
    PairedString &paired_quality,
    std::vector<int> &path,
    QualityMergingMap quality_merging_map
);