#include "mergepairs.h"
#include <unordered_map>
#include <string>
// #include <cstdio>

std::unordered_map<std::string, char> get_current_chars(
    unsigned int row,
    unsigned int col,
    std::unordered_map<std::string, std::string> &sequences,
    std::unordered_map<std::string, std::string> &qualities
) {
  return {
    {"sequence_forward", sequences["forward"][col]},
    {"quality_forward", qualities["forward"][col]},
    {"sequence_reverse", sequences["reverse"][row]},
    {"quality_reverse", qualities["reverse"][row]}
  };
}

// [[Rcpp::export]]
char get_merged_qualities(
    char &q1char,
    char &q2char,
    std::vector<std::vector<unsigned int>> &merged_qualities
) {
  unsigned int q1 = int(q1char) - PHRED_OFFSET;
  unsigned int q2 = int(q2char) - PHRED_OFFSET;
  unsigned int q_merged;
  if (q2 > q1) std::swap(q1, q2);
  q_merged = merged_qualities[q1 - 1][q2 - 1] + PHRED_OFFSET;
  return char(q_merged);
}

std::unordered_map<std::string, char> merge_forward_and_reverse_chars(
    std::unordered_map<std::string, char> chars,
    std::vector<std::vector<unsigned int>> merged_qualities_match,
    std::vector<std::vector<unsigned int>> merged_qualities_mismatch
) {
  std::unordered_map<std::string, char> merged;
  if (chars["sequence_forward"] == chars["sequence_reverse"]) {
    merged["sequence"] = chars["sequence_forward"];
    merged["quality"] = get_merged_qualities(
      chars["quality_forward"],
      chars["quality_reverse"],
      merged_qualities_match
    );
  } else {
    if (chars["quality_forward"] < chars["quality_reverse"]) {
      merged["sequence"] = chars["sequence_reverse"];
    } else {
      merged["sequence"] = chars["sequence_forward"];
    }
    merged["quality"] = get_merged_qualities(
      chars["quality_forward"],
      chars["quality_reverse"],
      merged_qualities_mismatch
    );
  }
  return merged;
}

std::unordered_map<std::string, std::string> merge_by_path_backtrack(
    std::vector<int> &path,
    std::unordered_map<std::string, std::string> &sequences,
    std::unordered_map<std::string, std::string> &qualities,
    std::vector<std::vector<unsigned int>> &merged_qualities_match,
    std::vector<std::vector<unsigned int>> &merged_qualities_mismatch
) {
  unsigned int ncol = sequences["forward"].size() + 1;
  unsigned int nrow = sequences["reverse"].size() + 1;
  unsigned int col = ncol - 1;
  unsigned int row = nrow - 1;
  unsigned int flat_id;
  std::unordered_map<std::string, std::string> aligned;
  std::unordered_map<std::string, char> unmerged_chars, merged_chars;
  std::string merged_sequence, merged_quality;
  
  while (col > 0 || row > 0) {
    flat_id = flat_index(col, row, ncol);
    switch(path[flat_id]) {
      case int('d'):
        col--; row--;
        unmerged_chars = get_current_chars(row, col, sequences, qualities);
        merged_chars = merge_forward_and_reverse_chars(
          unmerged_chars, merged_qualities_match, merged_qualities_mismatch
        );
        aligned["sequence"].push_back(merged_chars["sequence"]);
        aligned["quality"].push_back(merged_chars["quality"]);
        break;
      case int('l'):
        col--;
        aligned["sequence"].push_back(sequences["forward"][col]);
        aligned["quality"].push_back(qualities["forward"][col]);
        break;
      case int('u'):
        row--;
        aligned["sequence"].push_back(sequences["reverse"][row]);
        aligned["quality"].push_back(qualities["reverse"][row]);
        break;
      default:
        Rcpp::stop("Invalid backtracking value.");
    }
  }
  std::reverse(aligned["sequence"].begin(), aligned["sequence"].end());
  std::reverse(aligned["quality"].begin(), aligned["quality"].end());
  return aligned;
}

// [[Rcpp::export]]
std::unordered_map<std::string, std::string> align_seqs_and_quals(
    std::string &sequence_forward,
    std::string &sequence_reverse,
    std::string &quality_forward,
    std::string &quality_reverse,
    std::vector<int> &alignment_scores,
    std::vector<std::vector<unsigned int>> &merged_qualities_match,
    std::vector<std::vector<unsigned int>> &merged_qualities_mismatch
) {
  std::unordered_map<std::string, std::string> alignment;
  std::vector<std::vector<int>> scoring_matrix;
  std::unordered_map<std::string, std::vector<int>> score_and_path;
  std::unordered_map<std::string, std::string> sequences_umap = {
    {"forward", sequence_forward}, {"reverse", sequence_reverse}
  };
  std::unordered_map<std::string, std::string> qualities_umap = {
    {"forward", quality_forward}, {"reverse", quality_reverse}
  };
  std::unordered_map<std::string, int> alignment_scores_umap = {
    {"match", alignment_scores[0]},
    {"mismatch", alignment_scores[1]},
    {"gap_penalty", alignment_scores[2]}
  };;
  
  scoring_matrix = create_scoring_matrix(
    alignment_scores_umap["match"], alignment_scores_umap["mismatch"]
  );
  score_and_path = find_best_scoring_path(
    sequences_umap, scoring_matrix, alignment_scores_umap["gap_penalty"]
  );
  alignment = merge_by_path_backtrack(
    score_and_path["path"],
    sequences_umap,
    qualities_umap,
    merged_qualities_match,
    merged_qualities_mismatch
  );
  return(alignment);
}
