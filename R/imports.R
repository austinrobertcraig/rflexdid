#' @importFrom stats as.formula coef lm.fit lm.wfit model.matrix model.response
#' @importFrom stats pf pt qnorm qt terms vcov weighted.mean printCoefmat
#' @importFrom Matrix sparseMatrix
#' @importFrom sandwich vcovCL vcovHC
#' @importFrom methods as
#' @keywords internal
NULL

# Quiet R CMD check on `.data` (the ggplot2 / rlang pronoun used in plot()).
utils::globalVariables(c(".data"))
