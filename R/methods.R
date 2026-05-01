#' @export
print.flexdid <- function(x, digits = 4, ...) {
  cat("Flexible difference-in-differences regression\n")
  cat("  Specification:", x$specification, "\n")
  cat("  Outcome:      ", x$yvar, "\n")
  cat("  Treatment:    ", x$tx_var,
      " | Group:", x$group_var,
      " | Time:", x$time_var, "\n")
  cat("  Observations: ", x$nobs, "\n")
  if (x$vcov_type == "cluster") {
    cat("  VCE:           cluster (", x$cluster_var, "; ",
        x$n_clusters, " clusters)\n", sep = "")
  } else {
    cat("  VCE:           robust\n")
  }
  cat(sprintf("  R-squared:     %.4f   (adj. %.4f)\n", x$r2, x$r2_adj))
  if (!is.na(x$f)) {
    cat(sprintf("  F(%d, %d) = %.3f\n", x$f_df1, x$f_df2, x$f))
  }
  invisible(x)
}

#' @export
coef.flexdid <- function(object, ...) {
  object$coefficients
}

#' @export
vcov.flexdid <- function(object, ...) {
  object$vcov
}

#' @export
summary.flexdid <- function(object, digits = 4, ...) {
  keep <- object$pivot_keep
  beta <- object$coefficients[keep]
  V <- object$vcov[keep, keep, drop = FALSE]
  se <- sqrt(diag(V))
  tval <- beta / se
  df <- object$df_residual
  pval <- 2 * stats::pt(-abs(tval), df = df)
  ci_lo <- beta - stats::qt(0.975, df) * se
  ci_hi <- beta + stats::qt(0.975, df) * se
  table <- cbind(
    Estimate = beta,
    `Std. Error` = se,
    `t value` = tval,
    `Pr(>|t|)` = pval,
    `[95% CI lo]` = ci_lo,
    `[95% CI hi]` = ci_hi
  )
  out <- list(
    coefficients = table,
    call         = object$call,
    specification = object$specification,
    vcov_type    = object$vcov_type,
    n_clusters   = object$n_clusters,
    cluster_var  = object$cluster_var,
    nobs         = object$nobs,
    df_residual  = object$df_residual,
    r2           = object$r2,
    r2_adj       = object$r2_adj,
    f            = object$f,
    f_df1        = object$f_df1,
    f_df2        = object$f_df2
  )
  class(out) <- "summary.flexdid"
  out
}

#' @export
print.summary.flexdid <- function(x, digits = 4, ...) {
  cat("flexdid (", x$specification, ")\n", sep = "")
  cat("Observations:", x$nobs, " | Residual df:", x$df_residual, "\n")
  if (x$vcov_type == "cluster") {
    cat("VCE: cluster (", x$cluster_var, "; ", x$n_clusters,
        " clusters)\n", sep = "")
  } else {
    cat("VCE: robust\n")
  }
  cat(sprintf("R-squared: %.4f   adj.: %.4f\n", x$r2, x$r2_adj))
  if (!is.na(x$f)) {
    cat(sprintf("F(%d, %d) = %.3f\n", x$f_df1, x$f_df2, x$f))
  }
  cat("\nCoefficients (first 20 shown; use coef() / vcov() for full output):\n")
  rows <- min(20L, nrow(x$coefficients))
  printCoefmat(x$coefficients[seq_len(rows), , drop = FALSE],
               digits = digits, has.Pvalue = TRUE, signif.legend = FALSE)
  invisible(x)
}

#' @export
print.flexdid_atet <- function(x, digits = 4, ...) {
  type_label <- switch(x$type,
    overall    = "Overall ATET",
    byget      = "ATET by group X exposure time",
    byexposure = "ATET by exposure time",
    bycalendar = "ATET by calendar time",
    bycohort   = "ATET by treated cohort",
    bygroup    = "ATET by treated group"
  )
  if (!is.null(x$for_expr)) {
    type_label <- paste0(type_label, " for `", deparse(x$for_expr), "`")
  }
  cat(type_label, "\n", sep = "")
  cat("Observations: ", x$n_full,
      " | ATET sample: ", x$n_sub, "\n", sep = "")
  cat("Aggregation weight: ", x$aggregation_weight, "\n", sep = "")
  cat("\n")
  format_atet_table(x$tidy_table, digits = digits)
  if (!is.null(x$test_result)) {
    tr <- x$test_result
    cat("\n", tr$title, "\n", sep = "")
    cat("  H0: ", tr$h0, "\n", sep = "")
    cat(sprintf("  F(%d, %d) = %.3f   Prob > F = %.4f\n",
                tr$df1, tr$df2, tr$F, tr$p))
  }
  invisible(x)
}

# Minimal table formatter that doesn't depend on printCoefmat's column-order
# conventions. Each column is formatted independently with the chosen digits.
#' @keywords internal
#' @noRd
format_atet_table <- function(mat, digits = 4) {
  fmt <- function(z) {
    out <- formatC(z, digits = digits, format = "f", flag = " ")
    out[is.na(z)] <- "      NA"
    out
  }
  cols <- vapply(seq_len(ncol(mat)), function(j) fmt(mat[, j]),
                 FUN.VALUE = character(nrow(mat)))
  if (!is.matrix(cols)) cols <- matrix(cols, ncol = ncol(mat))
  colnames(cols) <- colnames(mat)
  rownames(cols) <- rownames(mat)
  print(noquote(cols))
}

#' @export
as.data.frame.flexdid_atet <- function(x, row.names = NULL, optional = FALSE, ...) {
  d <- as.data.frame(x$tidy_table)
  d$label <- rownames(d)
  rownames(d) <- NULL
  d <- d[, c("label", setdiff(colnames(d), "label"))]
  d
}

#' @export
plot.flexdid_atet <- function(x, ...) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for plot(); install it or use as.data.frame() instead.",
         call. = FALSE)
  }
  d <- as.data.frame(x)
  d$lo <- d[, "Estimate"] - 1.96 * d[, "Std. Error"]
  d$hi <- d[, "Estimate"] + 1.96 * d[, "Std. Error"]
  # Numeric x for line/scatter; falls back gracefully for non-numeric labels.
  d$.label_num <- suppressWarnings(as.numeric(
    gsub("[^0-9.\\-]", "", d$label)
  ))

  if (x$type == "overall") {
    message("plot() is not produced for type = 'overall' (single ATET).")
    return(invisible(NULL))
  }

  xlab <- switch(x$type,
    byexposure = "Exposure time",
    bycalendar = "Time period",
    bycohort   = "Treated cohort",
    bygroup    = "Treated group",
    byget      = "Group X exposure time"
  )
  ylab <- "Average treatment effect"

  p <- ggplot2::ggplot(d, ggplot2::aes(x = .data[[".label_num"]],
                                        y = .data[["Estimate"]]))
  if (x$type == "byexposure") {
    p <- p +
      ggplot2::geom_ribbon(ggplot2::aes(ymin = .data[["lo"]],
                                        ymax = .data[["hi"]]),
                           alpha = 0.2) +
      ggplot2::geom_line() +
      ggplot2::geom_point()
  } else {
    p <- p +
      ggplot2::geom_errorbar(ggplot2::aes(ymin = .data[["lo"]],
                                          ymax = .data[["hi"]]),
                             width = 0.2) +
      ggplot2::geom_point()
  }
  p +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed") +
    ggplot2::labs(x = xlab, y = ylab) +
    ggplot2::scale_x_continuous(
      breaks = unique(d$.label_num[!is.na(d$.label_num)])
    ) +
    ggplot2::theme_minimal()
}

`%||%` <- function(a, b) if (is.null(a)) b else a
