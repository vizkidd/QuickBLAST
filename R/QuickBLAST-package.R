#' @keywords internal
"_PACKAGE"

utils::globalVariables(c("seq_info", "seqids", "seqs", "hsps", "."))

## usethis namespace: start
#' @importFrom magrittr %>%
#' @importFrom dplyr %>%
#' @import RcppProgress
#' @import dplyr
#' @import BH
#' @import compiler
#' @import fs
#' @import utils
#' @import future
#' @import future.apply
#' @import future.callr
## usethis namespace: end