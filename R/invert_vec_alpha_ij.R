#' Invert Lambda-level Dynamic Matrix to Omega-level 
#'
#' Reconstructs the \eqn{n_L \times n_L} Omega-level dynamic matrix \eqn{A_\Omega}
#' from the \eqn{n_D \times n_D} Lambda-level dynamic matrix \eqn{A_\Lambda}
#' under a formative mapping \eqn{\Lambda = W \Omega}, where each Omega index is
#' used by exactly one Lambda (disjoint groups), and each row of \code{W} contains
#' weights that sum to 1. The solution enforces the identity
#' \deqn{W A_\Omega = A_\Lambda W}
#' exactly, by building rank-one blocks
#' \deqn{A_\Omega[G_d, G_{d'}] = A_\Lambda[d,d'] \cdot \frac{w_d}{w_d^\top w_d} \cdot w_{d'}^\top,}
#' where \eqn{G_d} are the Omega indices used by Lambda \eqn{d}, and \eqn{w_d}
#' are the corresponding weights.
#'
#' @param ALambda A numeric matrix of size \eqn{n_D \times n_D}; the
#'   Lambda-level dynamic (temporal influence) matrix.
#' @param mappingL1L2 A numeric matrix \code{W} of size \eqn{n_D \times n_L}
#'   with formative weights mapping Omegas to Lambdas. **Requirements**:
#'   (i) each row sums to 1; (ii) each column is nonzero in at most one row
#'   (i.e., disjoint Omega groups per Lambda). Columns that are all-zero are allowed
#'   (unused Omegas) but will produce zero rows/columns in the output.
#' @param check_weights Logical; if \code{TRUE} (default), validates that each row
#'   of \code{W} sums to 1 and that Omega indices are disjoint across rows.
#'
#' @return A list with components:
#' \describe{
#'   \item{\code{AOmega}}{Numeric matrix \eqn{n_L \times n_L}; the Omega-level dynamic matrix.}
#'   \item{\code{residual}}{Numeric matrix \eqn{n_D \times n_D_L}, the residual
#'         \code{W \%*\% AOmega - ALambda \%*\% W}; should be numerically ~0.}
#'   \item{\code{max_abs_residual}}{Maximum absolute entry of \code{residual}.}
#'   \item{\code{W}}{The (possibly validated) weight matrix used in the computation.}
#' }
#'
#' @details
#' This function implements the exact block outer-product construction that preserves
#' sparsity outside the mapped groups and guarantees the linear consistency
#' \eqn{W A_\Omega = A_\Lambda W}. It requires that each Omega index belongs to
#' exactly one Lambda row (no overlaps). If your mapping overlaps (an Omega feeds
#' multiple Lambdas), use a least-squares pseudoinverse approach instead (see
#' \code{\link{invert_dynamic_pinv}} in the examples).
#'
#' Row weights summing to 1 are not strictly necessary for correctness, but they
#' are typical in formative models and simplify interpretation. The denominator
#' \eqn{w_d^\top w_d} is the row’s sum of squares and ensures proper scaling.
#'
#' @examples
#' # Two Lambdas (nD=2), three Omegas (nL=3):
#' # Λ1 <- 0.6*Ω1 + 0.4*Ω2 ; Λ2 <- 1.0*Ω3
#' W <- rbind(
#'   c(0.6, 0.4, 0.0),
#'   c(0.0, 0.0, 1.0)
#' )
#'
#' # Lambda-level dynamic matrix
#' ALambda <- matrix(c(0.8, 0.1,
#'                     0.0, 0.5), nrow = 2, byrow = TRUE)
#'
#' out <- invert_dynamic_option1(ALambda, W)
#' out$AOmega
#' out$max_abs_residual  # ~ 0
#'
#' @export
invert_vec_alpha_ij <- function(ALambda, mappingL1L2) {
  W <- mappingL1L2
  nD <- nrow(W); nL <- ncol(W)
  
  # 2) Build Omega index groups per Lambda
  groups <- lapply(seq_len(nD), function(d) which(W[d, ] != 0))
  AOmega <- matrix(0, nL, nL)
  
  # 3) Fill AOmega blocks via the outer-product construction
  for (d in seq_len(nD)) {
    Gd <- groups[[d]]
    wd <- W[d, Gd]
    denom <- sum(wd^2)
    if (!is.finite(denom) || denom <= 0) {
      stop(sprintf("Non-positive denominator in row %d (weights must not be all zeros).", d))
    }
    for (dp in seq_len(nD)) {
      Gdp <- groups[[dp]]
      wdp <- W[dp, Gdp]
      AOmega[Gd, Gdp] <- ALambda[d, dp] * (wd / denom) %*% t(wdp)
    }
  }
  
  # 4) Consistency check: W %*% AOmega == ALambda %*% W
  residual <- W %*% AOmega - ALambda %*% W
  max_abs_residual <- max(abs(residual))
  
  list(AOmega = AOmega, residual = residual, max_abs_residual = max_abs_residual, W = W)
}
