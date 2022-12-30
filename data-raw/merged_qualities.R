get_merged_qualities_vector <- function(type) {
  q_dt <- read.delim(system.file(
    "merge_tables",
    paste0("himap_mergepairs_", type, "_qs.txt"),
    package = "rexmap"
  )) |> as.data.table()
  colnames(q_dt) <- c("q1", "q2", "qm")
  q_filt_dt <- q_dt[q1 > 0 & q2 > 0 & q1 >= q2][order(q1, q2)]
  q_filt_dt[, qm]
}

get_merged_qualities_2d <- function(type) {
  merged_quality_vector <- get_merged_qualities_vector(type)
  merged_matching_qualities <- vector("list", 41)
  k <- 1
  for (i in 1:41) {
    merged_matching_qualities[[i]] <- vector("integer", i)
    for (j in 1:i) {
      merged_matching_qualities[[i]][[j]] <- merged_quality_vector[k]
      k <- k + 1
    }
  }
  merged_matching_qualities
}

merged_qualities <- list(
  "match" = get_merged_qualities_2d("match"),
  "mismatch" = get_merged_qualities_2d("mismatch")
)

# For matrix display use:
#  xmat <- matrix(x.dt[q1 > 0 & q2 > 0, qm], ncol = 41)

usethis::use_data(merged_qualities, overwrite = TRUE, internal = TRUE)
