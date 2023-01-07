#include "SeqAndQual.h"

SeqAndQualChar::SeqAndQualChar(char sequence, char quality) :
  sequence(sequence), quality(quality) {}

SeqAndQualChar SeqAndQualChar::set(char sequence, char quality) {
  this->sequence = sequence;
  this->quality = quality;
  return *this;
};

PairedSeqAndQualChar::PairedSeqAndQualChar(
    char sequence_forward_char,
    char sequence_reverse_char,
    char quality_forward_char,
    char quality_reverse_char
  ) : forward(sequence_forward_char, quality_forward_char),
  reverse(sequence_reverse_char, quality_reverse_char) {
  is_equal_seq = sequence_forward_char == sequence_reverse_char;
}

SeqAndQualChar PairedSeqAndQualChar::merge(
    QualityMergingMap merge_map
) {
  SeqAndQualChar merged_char;
  merged_char.sequence =
    forward.quality >= reverse.quality ? forward.sequence : reverse.sequence;
  merged_char.quality = merge_map.merge_quality_chars(
    forward.quality,
    reverse.quality,
    is_equal_seq
  );
  return merged_char;
};