#include "QualityMergingMap.h"

class SeqAndQualChar {
public:
  char sequence;
  char quality;
  SeqAndQualChar() {}
  SeqAndQualChar(char sequence, char quality);
  SeqAndQualChar set(char sequence, char quality);
};

class PairedSeqAndQualChar {
public:
  SeqAndQualChar forward;
  SeqAndQualChar reverse;
  bool is_equal_seq;
  PairedSeqAndQualChar() {}
  PairedSeqAndQualChar(
    char sequence_forward_char,
    char sequence_reverse_char,
    char quality_forward_char,
    char quality_reverse_char
  );
  SeqAndQualChar merge(QualityMergingMap merge_map);
};