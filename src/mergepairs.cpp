// #include <Rcpp.h>
#include "rexmap.h"
using namespace Rcpp;

// Precomputer posterior q scores for matching and mismatching bases
std::vector< std::vector<int> > post_q_match;
std::vector< std::vector<int> > post_q_mismatch;
std::string posterior_match_file;
std::string posterior_mismatch_file;

// Scoring matrix
int c_score[5][5];

// Load pre-computed posterior quality scores for matching and mismatching
// bases.
// [[Rcpp::export]]
std::vector< std::vector<int> > load_posterior (std::string filename) {
  // Load pre-computer posterior Q scores into a vector of vector of ints.
  std::vector< std::vector<int> > out;
  std::vector<int> row;
  std::ifstream fin(filename.c_str());
  int prev_q1 = 0, q1, q2, q;
  while (fin >> q1 >> q2 >> q) {
    if (q1 == prev_q1) {
      // same first index
      row.push_back(q);
    } else {
      // new first index
      out.push_back(row);
      row.clear();
      row.push_back(q);
    }
    prev_q1 = q1;
  }
  // Add last row
  out.push_back(row);
  return(out);
}

// Modified dada2 nwalign_endsfree to allow Ns and
void nt2int(char *oseq, const char *iseq) {
  long unsigned int i, len = strlen(iseq);
  for (i = 0; i < len; i++, iseq++, oseq++) {
    switch (*iseq) {
    case 'A':
      *oseq = 1;
      break;
    case 'C':
      *oseq = 2;
      break;
    case 'G':
      *oseq = 3;
      break;
    case 'T':
      *oseq = 4;
      break;
    case 'N':
      *oseq = 5;
      break;
    case '-':
      *oseq = '-';
      break;
    default:
      //Rprintf("invalid character in input:%c.\n",*iseq);
      *oseq = '\0';
    }
  }
  *oseq = '\0';
  return;
}

void qs2int(char *oseq, const char *iseq) {
  long unsigned int i, len = strlen(iseq);
  for (i = 0; i < len; i++, iseq++, oseq++) {
    *oseq = int(*iseq);
  }
  *oseq = '\0';
  return;
}

void int2qs(char *oseq, const char *iseq) {
  long unsigned int i, len = strlen(iseq);
  for (i = 0; i < len; i++, iseq++, oseq++) {
    *oseq = char(*iseq);
  }
  *oseq = '\0';
  return;
}

void int2nt(char *oseq, const char *iseq) {
  long unsigned int i, len = strlen(iseq);
  for (i = 0; i < len; i++, iseq++, oseq++) {
    switch (*iseq) {
    case 1:
      *oseq = 'A';
      break;
    case 2:
      *oseq = 'C';
      break;
    case 3:
      *oseq = 'G';
      break;
    case 4:
      *oseq = 'T';
      break;
    case 5:
      *oseq = 'N';
      break;
    case '-':
      *oseq = '-';
      break;
    default:
      break;
    }
  }
  *oseq = '\0';
  return;
}

char **himap_nwalign_endsfree(const char *s1, const char *s2,
                              const char *q1, const char *q2,
                              int score[5][5], int gap_p) {
  static size_t nnw = 0;
  long unsigned int i, j;
  long unsigned int l, r;
  // int band = -1;
  long unsigned int len1 = strlen(s1);
  long unsigned int len2 = strlen(s2);
  int diag, left, up;

  long unsigned int nrow = len1+1;
  long unsigned int ncol = len2+1;
  int *d = (int *) malloc(nrow * ncol * sizeof(int)); //E
  int *p = (int *) malloc(nrow * ncol * sizeof(int)); //E
  if(d == NULL || p == NULL) Rcpp::stop("Memory allocation failed.");

  // Fill out left columns of d, p.
  for (i = 0; i <= len1; i++) {
    d[i*ncol] = 0; // ends-free gap
    p[i*ncol] = 3;
  }

  // Fill out top rows of d, p.
  for (j = 0; j <= len2; j++) {
    d[j] = 0; // ends-free gap
    p[j] = 2;
  }

  // Calculate left/right-bands in case of different lengths
  //int lband, rband;
  // if(len2 > len1) {
  //   lband = band;
  //   rband = band+len2-len1;
  // } else if(len1 > len2) {
  //   lband = band+len1-len2;
  //   rband = band;
  // } else {
  //   lband = band;
  //   rband = band;
  // }

  // Fill out band boundaries of d.
  // if(band>=0 && (band<len1 || band<len2)) {
  //   for(i=0;i<=len1;i++) {
  //     if(i-lband-1 >= 0) { d[i*ncol + i-lband-1] = -9999; }
  //     if(i+rband+1 <= len2) { d[i*ncol + i+rband+1] = -9999; }
  //   }
  // }

  // Fill out the body of the DP matrix.
  for (i = 1; i <= len1; i++) {
    // if(band>=0) {
    //   l = i-lband; if(l < 1) { l = 1; }
    //   r = i+rband; if(r>len2) { r = len2; }
    //} else {
    l=1; r=len2;
    //  }

    for (j = l; j <= r; j++) {
      // Score for the left move.
      if (i == len1) {
        left = d[i*ncol + j-1]; // Ends-free gap.
      } else {
        left = d[i*ncol + j-1] + gap_p;
      }

      // Score for the up move.
      if (j == len2) {
        up = d[(i-1)*ncol + j]; // Ends-free gap.
      } else {
        up = d[(i-1)*ncol + j] + gap_p;
      }

      // Score for the diagonal move.
      diag = d[(i-1)*ncol + j-1] + score[s1[i-1]-1][s2[j-1]-1];

      // Break ties and fill in d,p.
      if (up >= diag && up >= left) {
        d[i*ncol + j] = up;
        p[i*ncol + j] = 3;
      } else if (left >= diag) {
        d[i*ncol + j] = left;
        p[i*ncol + j] = 2;
      } else {
        d[i*ncol + j] = diag;
        p[i*ncol + j] = 1;
      }
    }
  }

  char *al0 = (char *) malloc((len1+len2+1) * sizeof(char));
  char *al1 = (char *) malloc((len1+len2+1) * sizeof(char));
  if(al0 == NULL || al1 == NULL) Rcpp::stop("Memory allocation failed.");

  char *qs0 = (char *) malloc((len1+len2+1) * sizeof(char));
  char *qs1 = (char *) malloc((len1+len2+1) * sizeof(char));
  if(qs0 == NULL || qs1 == NULL) Rcpp::stop("Memory allocation failed.");

  // Trace back over p to form the alignment.
  size_t len_al = 0;
  i = len1;
  j = len2;

  while ( i > 0 || j > 0 ) {
    switch ( p[i*ncol + j] ) {
    case 1:
      i = i - 1;
      j = j - 1;
      al0[len_al] = s1[i];
      al1[len_al] = s2[j];
      qs0[len_al] = q1[i];
      qs1[len_al] = q2[j];
      break;
    case 2:
      j = j - 1;
      al0[len_al] = '-';
      al1[len_al] = s2[j];
      qs0[len_al] = ' ';
      qs1[len_al] = q2[j];
      break;
    case 3:
      i = i - 1;
      al0[len_al] = s1[i];
      al1[len_al] = '-';
      qs0[len_al] = q1[i];
      qs1[len_al] = ' ';
      break;
    default:
      Rcpp::stop("N-W Align out of range.");
    }
    len_al++;
  }
  al0[len_al] = '\0';
  al1[len_al] = '\0';
  qs0[len_al] = '\0';
  qs1[len_al] = '\0';

  // Allocate memory to alignment strings.
  char **al = (char **) malloc( 4 * sizeof(char *) ); //E
  if (al == NULL)  Rcpp::stop("Memory allocation failed.");
  al[0] = (char *) malloc(len_al+1); //E
  al[1] = (char *) malloc(len_al+1); //E
  al[2] = (char *) malloc(len_al+1); //E
  al[3] = (char *) malloc(len_al+1); //E
  // char *merged_seq = (char *) malloc(len_al+1);
  // char *merged_qual = (char *) malloc(len_al+1);
  if (al[0] == NULL || al[1] == NULL)  Rcpp::stop("Memory allocation failed.");

  // Reverse the alignment strings (since traced backwards).
  for (i=0;i<len_al;i++) {
    al[0][i] = al0[len_al-i-1];
    al[1][i] = al1[len_al-i-1];
    al[2][i] = qs0[len_al-i-1];
    al[3][i] = qs1[len_al-i-1];
  }
  al[0][len_al] = '\0';
  al[1][len_al] = '\0';
  al[2][len_al] = '\0';
  al[3][len_al] = '\0';

  // Free allocated memory
  free(d);
  free(p);
  free(al0);
  free(al1);
  free(qs0);
  free(qs1);

  nnw++;
  return al;
}

char **himap_merge_alignment(char** al) {
  // al[0] seq1, al[1] seq2, al[2] qual1, al[3] qual 2
  long unsigned int i, len = strlen(al[0]);
  int q1, q2, q;
  char **out = (char **) malloc(2*sizeof(char*));
  // char tmp;
  out[0] = (char *) malloc(len+1); // merged sequence
  out[1] = (char *) malloc(len+1); // merged qual scores

  for (i=0; i<len; i++) {
    q1 = int(al[2][i]) - PHRED_OFFSET;
    q2 = int(al[3][i]) - PHRED_OFFSET;
    // std::cout << "q1: " << q1 << ", q2: " << q2 << "\n";
    //q = post_q_match[q1][q2] + PHRED_OFFSET;
    // std::cout << "q: " << q << "\n";
    if (al[0][i] == al[1][i]) { // matching letters, so neither is -
      if (al[0][i] == 'N') {
        out[0][i] = 'N';
        out[1][i] = '\"';
      }
      else {
        out[0][i] = al[0][i];
        q = post_q_match[q1][q2] + PHRED_OFFSET;
        // std::cout << "q: " << q << "\n";
        out[1][i] = (char)q;
      }
    }
    else { // mismatch
      if (al[0][i] == '-') { // fwd is missing, use just reverse
        out[0][i] = al[1][i];
        out[1][i] = al[3][i];
      } else if (al[1][i] == '-') { // rev is missing, use forward
        out[0][i] = al[0][i];
        out[1][i] = al[2][i];
      } else { // mismatch and not indel
        q = post_q_mismatch[q1][q2] + PHRED_OFFSET;
        out[1][i] = (char)q;
        if (q1 >= q2) {
          out[0][i] = al[0][i];
        } else {
          out[0][i] = al[1][i];
        }
      }
    }
  }
  out[0][len] = '\0';
  out[1][len] = '\0';
  return(out);
}

// [[Rcpp::export]]
Rcpp::CharacterVector nwalign_endsfree_test (std::string s1, std::string s2, std::string q1, std::string q2,
                                        int match=5, int mismatch=-2) {
  // Make  c-style 2d score array
  int i, j;
  int c_score[5][5];
  for(i=0;i<5;i++) {
    for(j=0;j<5;j++) {
      if(i==j) {
        if (i==4) {
          // both i and j are N, declare it mismatch
          c_score[i][j] = mismatch;
        } else {
          c_score[i][j] = match;
        }
      }
      else {
        c_score[i][j] = mismatch;
      }
    }
  }
  char *seq1 = (char *) malloc(s1.size()+1);
  char *seq2 = (char *) malloc(s2.size()+1);
  if (seq1 == NULL || seq2 == NULL)  Rcpp::stop("Memory allocation failed.");
  nt2int(seq1, s1.c_str());
  nt2int(seq2, s2.c_str());
  const char *qs1 = q1.c_str();
  const char *qs2 = q2.c_str();
  int gap_p = -7;
  char** al = himap_nwalign_endsfree(seq1, seq2, qs1, qs2, c_score, gap_p);
  int2nt(al[0], al[0]);
  int2nt(al[1], al[1]);
  Rcpp::CharacterVector rval;
  rval.push_back(std::string(al[0]));
  rval.push_back(std::string(al[1]));
  return(rval);
}


// [[Rcpp::export]]
Rcpp::CharacterVector C_mergepairs(std::string s1, std::string s2,
                                std::string q1, std::string q2,
                                std::string posterior_match_file,
                                std::string posterior_mismatch_file,
                                int match, int mismatch, int gap_p,
                                double min_pct_sim, int min_aln_len) {

  // Alignment filter parameters
  int aln_matches = 0; // number of matching letters in the alignment
  int left_end = 0, right_end = 0; // lengths of ---- ends
  bool still_left = TRUE;
  int len_al, len_overlap;
  double pct_sim;

  // Load pre-computed posterior quality scores
  post_q_match = load_posterior(posterior_match_file);
  post_q_mismatch = load_posterior(posterior_mismatch_file);

  // Make  c-style 2d score array
  int i, j;
  int c_score[5][5];
  for(i=0;i<5;i++) {
    for(j=0;j<5;j++) {
      if(i==j) {
        if (i==4) {
          // both i and j are N, declare it mismatch
          c_score[i][j] = mismatch;
        } else {
          c_score[i][j] = match;
        }
      }
      else {
        c_score[i][j] = mismatch;
      }
    }
  }

  // Make integer-ized c-style sequence strings
  char **al;
  char **merged;
  char *seq1 = (char *) malloc(s1.size()+1);
  char *seq2 = (char *) malloc(s2.size()+1);
  if (seq1 == NULL || seq2 == NULL)  Rcpp::stop("Memory allocation failed.");
  nt2int(seq1, s1.c_str());
  nt2int(seq2, s2.c_str());
  const char *qs1 = q1.c_str();
  const char *qs2 = q2.c_str();

  // Perform alignment and convert back to ACGT
  al = himap_nwalign_endsfree(seq1, seq2, qs1, qs2, c_score, gap_p);
  int2nt(al[0], al[0]);
  int2nt(al[1], al[1]);

  len_al = static_cast<int>(strlen(al[0]));

  // Rcpp::stop("Checkpoint1.");

  // Calculate the number of matching letters in the alignment
  // Since it is ends-free - vs anything will be a mismatch automatically.
  for (i=0; i<len_al; i++) {
    // "letter" (int) of the first alignment string
    if (al[0][len_al-i-1] == al[1][len_al-i-1]) {
      aln_matches++;
    }
    if (still_left) {
      if (al[0][len_al-i-1] == '-') {
        left_end++;
      } else {
        still_left = FALSE;
      }
    }
    if (al[1][len_al-i-1] == '-') {
      right_end++;
    } else {
      right_end = 0;
    }
  }

  // std::cout << "Here.\n";
  // Now find the length of the aligned region
  len_overlap = len_al - left_end - right_end;

  // Calculate percentage similarity
  pct_sim = aln_matches/(double)len_overlap;

  // Clean up
  free(seq1);
  free(seq2);
  // Now merge the aligned sequences
  merged = himap_merge_alignment(al);

  free(al[0]);
  free(al[1]);
  free(al[2]);
  free(al[3]);
  free(al);
  // Generate R-style return vector
  Rcpp::CharacterVector rval;
  rval.push_back(std::string(merged[0]));
  rval.push_back(std::string(merged[1]));
  rval.push_back(std::string(pct_sim >= min_pct_sim ? "TRUE" : "FALSE"));
  rval.push_back(std::string(len_overlap >= min_aln_len ? "TRUE" : "FALSE"));
  // rval.push_back(std::string(al[0]));
  // rval.push_back(std::string(al[1]));
  free(merged[0]);
  free(merged[1]);
  return(rval);
}



