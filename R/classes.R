#' The object class returned by \code{\link{blast}}
#'
#' A multi-item list with the following named values:
#' \itemize{
#'  \item{$alignments: }{Best scored alignments with the reference database.}
#'  \item{$cp: }{Copy number table for 16S genes in the reference database.}
#'  \item{$parameters: }{Parameters used for the BLAST alignment.}
#' }
#'
#' @seealso \code{\link{blast}}
#'
#' @name blast-class
#' @rdname blast-class
setClass('blast', contains = 'list')
setMethod('show', 'blast', function(object){
  cat('blast-class: object with BLAST results', fill=T)
  cat('This class is used as an input to rexmap::abundance() function.', fill=T)
  cat('Contains data.tables:', fill=T)
  cat('* BLAST reference database alignments: $alignments', fill=T)
  cat('* Copy number table: $cp', fill=T)
  cat('* Sequence table: $sequences', fill=T)
  cat('* BLAST alignment parameters: $parameters', fill=T)
})
