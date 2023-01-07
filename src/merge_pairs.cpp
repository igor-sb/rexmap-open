#include <Rcpp.h>
#include "PairedString.h"

class SeqAndQualChar {
public:
    char sequence;
    char quality;
    SeqAndQualChar() {}
    SeqAndQualChar(char sequence, char quality) :
      sequence(sequence), quality(quality) {}
    SeqAndQualChar set(char sequence, char quality) {
      this->sequence = sequence;
      this->quality = quality;
      return *this;
    }
};

class PairedSeqAndQualChar {
public:
  SeqAndQualChar forward;
  SeqAndQualChar reverse;
  PairedSeqAndQualChar(
    char sequence_forward_char,
    char sequence_reverse_char,
    char quality_forward_char,
    char quality_reverse_char
  ) : forward(sequence_forward_char, quality_forward_char),
      reverse(sequence_reverse_char, quality_reverse_char) {}
};

PairedSeqAndQualChar locate_sequence_and_quality_char(
    unsigned int row,
    unsigned int column,
    PairedString &sequences,
    PairedString &qualities
) {
  return PairedSeqAndQualChar(
    sequences.forward[column],
    sequences.reverse[row],
    qualities.forward[column],
    qualities.reverse[row]
  );
}

class QualityMergingMap {
private:
  std::vector<std::vector<unsigned int>> &match;
  std::vector<std::vector<unsigned int>> &mismatch;
public:
  QualityMergingMap(
    std::vector<std::vector<unsigned int>> &match,
    std::vector<std::vector<unsigned int>> &mismatch
  ) : match(match), mismatch(mismatch) {}
  
  std::vector<std::vector<unsigned int>> &get_map(bool is_match) {
    return is_match ? match : mismatch;
  }
  
  int get_merged_quality
  
  char merge_quality_chars(
    char &quality_forward_char,
    char &quality_reverse_char,
    bool is_match,
    int PHRED_OFFSET = 33
  ) {
    unsigned int quality_forward_int = int(quality_forward_char) - PHRED_OFFSET;
    unsigned int quality_reverse_int = int(quality_reverse_char) - PHRED_OFFSET;
    unsigned int quality_merged_int;
    if (quality_reverse_int > quality_forward_int) {
      std::swap(quality_forward_int, quality_reverse_int);
    }
    quality_merged_int = quality_merging_map
      [quality_forward_int - 1][quality_reverse_int - 1] + PHRED_OFFSET;
    return char(quality_merged_int);
    
  }
};

// [[Rcpp::export]]
char merge_aligned_quality_char(
    char &quality_forward_char,
    char &quality_reverse_char,
    std::vector<std::vector<unsigned int>> &quality_merging_map,
    int PHRED_OFFSET = 33
) {
  unsigned int quality_forward_int = int(quality_forward_char) - PHRED_OFFSET;
  unsigned int quality_reverse_int = int(quality_reverse_char) - PHRED_OFFSET;
  unsigned int quality_merged_int;
  if (quality_reverse_int > quality_forward_int) {
    std::swap(quality_forward_int, quality_reverse_int);
  }
  quality_merged_int = quality_merging_map
    [quality_forward_int - 1][quality_reverse_int - 1] + PHRED_OFFSET;
  return char(quality_merged_int);
}

SeqAndQualChar merge_forward_and_reverse_chars(
    PairedSeqAndQualChar aligned_char,
    std::vector<std::vector<unsigned int>> &matching_quality_merging_map,
    std::vector<std::vector<unsigned int>> &mismatching_quality_merging_map
) {
  char merged_quality_char;
  bool is_match;
  SeqAndQualChar merged_char;
  
  if (aligned_char.forward.quality >= aligned_char.reverse.quality) {
    merged_char.sequence = aligned_char.forward.sequence;
  } else {
    merged_char.sequence = aligned_char.reverse.sequence;
  }
  
  is_seq_match = aligned_char.forward.sequence == aligned_char.reverse.sequence;
  merged_char.quality = merge_aligned_quality_char(
    aligned_char.forward.quality,
    aligned_char.reverse.quality,
    matching_quality_merging_map
  );
  return merged_char;
}

bool is_same_sequence_char(
    std::unordered_map<std::string, char> &unmerged_chars
) {
  return unmerged_chars["sequence_forward"] == 
    unmerged_chars["sequence_reverse"];
}

bool is_path_start(unsigned int &row, unsigned int &col) {
  return col == 0 && row == 0;
}

MergedAlignment merge_by_path_backtrack(
    std::vector<int> &path,
    PairedString &sequences,
    PairedString &qualities,
    Vector2d<unsigned int> &merged_qualities_match,
    Vector2d<unsigned int> &merged_qualities_mismatch
) {
  unsigned int ncol = sequences.forward.size() + 1;
  unsigned int nrow = sequences.reverse.size() + 1;
  unsigned int flat_id, col = ncol - 1, row = nrow - 1;
  MergedAlignment aligned;
  PairedSeqAndQualChar unmerged_chars;
  SeqAndQualChar merged_char;

  while (!is_path_start(row, col)) {
    flat_id = flatten_index(col, row, ncol);
    switch(path[flat_id]) {
      case int('d'):
        col--; row--;
        unmerged_chars = find_char_at_row_column(
          row, col, sequences, qualities
        );
        merged_char = merge_forward_and_reverse_chars(
          unmerged_chars, merged_qualities_match, merged_qualities_mismatch
        );
        aligned.sequence.push_back(merged_chars["sequence"]);
        aligned.quality.push_back(merged_chars["quality"]);
        aligned.overlap_length++;
        if (is_same_sequence_char(unmerged_chars)) aligned.overlap_matches++;
        break;
      case int('l'):
        col--;
        aligned.sequence.push_back((sequences.forward)[col]);
        aligned.quality.push_back((qualities.forward)[col]);
        if (row > 0) aligned.overlap_length++;
        break;
      case int('u'):
        row--;
        aligned.sequence.push_back((sequences.reverse)[row]);
        aligned.quality.push_back((qualities.reverse)[row]);
        if (col < ncol - 1) aligned.overlap_length++;
        break;
      default:
        Rcpp::stop("Invalid backtracking value.");
    }
  }
  std::reverse(aligned.sequence.begin(), aligned.sequence.end());
  std::reverse(aligned.quality.begin(), aligned.quality.end());
  return aligned;
}

// [[Rcpp::export]]
MergedAlignment align_seqs_and_quals(
    Rcpp::List &sequences,
    Rcpp::List &qualities,
    Rcpp::IntegerVector &alignment_scores,
    Vector2d<unsigned int> &merged_qualities_match,
    Vector2d<unsigned int> &merged_qualities_mismatch
) {
  MergedAlignment alignment;
  Vector2d<int> scoring_matrix;
  Hashmap<std::vector<int>> score_and_path;

  PairedString sequences_pair(sequences["forward"], sequences["reverse"]);
  PairedString qualities_pair(qualities["forward"], qualities["reverse"]);
  scoring_matrix = create_scoring_matrix(
    alignment_scores["match"], alignment_scores["mismatch"]
  );
  score_and_path = find_best_scoring_path(
    sequences_pair, scoring_matrix, alignment_scores["gap_penalty"]
  );
  alignment = merge_by_path_backtrack(
    score_and_path["path"],
    sequences_pair,
    qualities_pair,
    merged_qualities_match,
    merged_qualities_mismatch
  );
  return(alignment);
}
