#include <Rcpp.h>
// using namespace std;

// [[Rcpp::plugins(cpp11)]]

// Nucleotide-to-integer conversion map
static std::map<char, int> nt2int = {{'A',0}, {'C', 1}, {'G', 2}, {'T', 3}, {'N', 4}, {'-', 5}, {' ', 5}};
static std::map<int, char> int2nt = {{0, 'A'}, {1, 'C'}, {2, 'G'}, {3, 'T'}, {4, 'N'}, {5, '-'}};

// [[Rcpp::export]]
std::map<int, std::vector<int> > value_map (std::vector<int> x) {
  // input [1, 2, 3, 4, 1, 2, 5, 1, 6, 4, 7, 8, 8]
  //        0  1  2  3  4  5  6  7  8  9 10 11 12
  // output: 1-1=0: [0, 4, 7],
  //         2-1=1, [1, 5],
  //         3-1=2, [2],
  //         4-1=3, [3, 9],
  //         5-1=4, [6],
  //         6-1=7, [8],
  //         7-1=6, [10]
  // actually subtract 1 in the input
  std::map<int, std::vector<int> > xmap;
  int xmap_size = static_cast<int>(xmap.size());
  int len = x.size();
  for (int i=0; i<len; ++i) {
    if (x[i] == xmap_size) { // add new element
      xmap[x[i]] = {i};
    } else { // this element already exists. append.
      xmap[x[i]].push_back(i);
    }
  }
  return xmap;
}

// [[Rcpp::export]]
std::map<int, std::vector<int> > partid_to_fastqid (std::vector<int> dada_map, std::vector<int> derep_map) {
  // Generate map partitionID -> fastqID
  // inputs are maps from paritionID -> uniqueIDs and uniqueID -> fastqIDs
  std::map<int, std::vector<int> > pid_to_uid = value_map(dada_map);
  std::map<int, std::vector<int> > uid_to_fid = value_map(derep_map);
  // Now combine these into paritionID -> fastqIDs
  std::map<int, std::vector<int> > pid_to_fid;
  for (std::map<int, std::vector<int> >::const_iterator ptu_it = pid_to_uid.begin(); ptu_it != pid_to_uid.end(); ++ptu_it) {
    pid_to_fid[ptu_it->first] = {};
    for (std::vector<int>::const_iterator u_it = ptu_it->second.begin(); u_it != ptu_it->second.end(); ++u_it) {
      pid_to_fid[ptu_it->first].insert(end(pid_to_fid[ptu_it->first]),
                                       begin(uid_to_fid[*u_it]),
                                       end(uid_to_fid[*u_it]));
    }
  }
  return pid_to_fid;
}

// [[Rcpp::export]]
std::string consensus_sequence (std::vector<std::string> seqs) {
  std::string cons;
  // Number of sequences in a vector.
  int n = seqs.size();
  // Find the longest sequence and use that to preallocate vector.
  int len = 0;
  std::vector<int> lens(n);
  int ni = 0;
  for (auto s : seqs) {
    lens[ni] = s.length();
    if (lens[ni] > len) len = lens[ni];
    ni++;
  }
  // Preallocate output vector as a flattened array, "row" = ACGTN- (6), "col"=position
  // transpose for a conventional form. this is faster for finding max.
  int* cons_mat = new int[6*len];
  // To do later: Check for successful memory allocation from the heap.
  int i, j, k, s;
  for (k=0; k<6*len; ++k) {
    cons_mat[k] = 0;
  }
  // i = seqs index, j = index over elements in each seq
  for (i=0; i<n; ++i) {
    for (j=0; j<lens[i]; ++j) { // iterate seqs[i] character [j]
      s = nt2int[seqs[i][j]]; // integer value of a nt
      k = j*6 + s; // k = cons_mat flattened index
      cons_mat[k]++;
    }
    // Fill up remaining elements with 5 ('-' or ' ' symbols)
    for (j=lens[i]; j<len; ++j) cons_mat[j*6 + 5]++;
  }
  // Debug: print consensus matrix
  // cout << "Consensus matrix" << "\n";
  // cout << ".: A C G T N -\n";
  // for (j=0; j<len; ++j) {
  //   cout << j << ": ";
  //   for (i=0; i<6; ++i) {
  //     cout << cons_mat[j*6+i] << " ";
  //   }
  //   cout << "\n";
  // }

  // Find maximum values for each position
  int max_id;
  for (j=0; j<len; ++j) {
    max_id = std::distance(cons_mat+j*6, std::max_element(cons_mat+j*6, cons_mat+j*6+6));
    // cout << "max value for pos " << j << " = " << max_id << "\n";
    cons += int2nt[max_id]; // append character to string
  }
  // Release memory for dynamically allocated array
  free(cons_mat);

  // Return consensus string
  return cons;
}

// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically
// run after the compilation.
//

/*** R
*/
