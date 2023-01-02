#include <RcppCommon.h>

#include "rexmap.h"

namespace Rcpp {
  template <>
  SEXP wrap(const MergedAlignment& x);
}

#include <Rcpp.h>

namespace Rcpp {
  template <>
  SEXP wrap(const MergedAlignment& x) {
    return Rcpp::wrap(
      Rcpp::List::create(
        Rcpp::Named("sequence") = x.sequence,
        Rcpp::Named("quality") = x.quality,
        Rcpp::Named("overlap_length") = x.overlap_length,
        Rcpp::Named("overlap_matches") = x.overlap_matches
      )
    );
  };
}
using namespace Rcpp;
