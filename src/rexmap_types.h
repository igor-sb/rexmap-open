#include <RcppCommon.h>

#include "rexmap.h"

namespace Rcpp {
  template <> SEXP wrap(const MergedAlignment& x);
}

#include <Rcpp.h>

namespace Rcpp {
  template <>
  SEXP wrap(const MergedAlignment& x) {
    return wrap(
      List::create(
        _["sequence"] = x.sequence,
        _["quality"] = x.quality,
        _["overlap_length"] = x.overlap_length,
        _["overlap_matches"] = x.overlap_matches
      )
    );
  };
}


using namespace Rcpp;
