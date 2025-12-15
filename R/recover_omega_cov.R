#' Recover Omega-level covariance from Lambda-level Cholesky
#'
#' @param Llambda_int nD x nD lower-triangular Cholesky for intercept block (Lambda)
#' @param Llambda_slope nD x nD lower-triangular Cholesky for slope block (Lambda). Can be NULL if no slopes.
#' @param mappingL1L2 nD x nL matrix; entry (d, l) is the weight ω_{d,l} of Omega l in Lambda d.
#                   If you only have a binary mapping, put 1's and set weights_normalize=TRUE.
#' @param method  "structured" (default) or "pseudoinverse"
#' @param rho_int scalar or length nD vector of correlations within Omega block of each Lambda (intercepts)
#' @param rho_slope scalar or length nD vector (slopes). Ignored if Llambda_slope is NULL.
#' @param intercept_var_fixed if TRUE, force Var(Λ_d intercept)=1 even if Llambda_int implies something else.
#' @param weights_normalize if TRUE, normalize each row of mapping to sum to 1 (formative "weights sum to 1").
#'
#' @returns  A list with BOmega_int, BOmega_slope (if slopes), W, BLambda_int, BLambda_slope, diagnostics.
#' @export
#'
#' @examples  Example: nD = 2 Lambdas, nL = 3 Omegas
#' Λ1 uses Ω1, Ω2; Λ2 uses Ω3
#' mappingL1L2 <- rbind(
#'   c(1, 1, 0),   #' Λ1 <- Ω1, Ω2
#'   c(0, 0, 1)    #' Λ2 <- Ω3
#' )
#' 
#' #' Cholesky (lower-triangular) for Λ intercepts and slopes
#' Llambda_int <- matrix(c(
#'   1.0, 0.0,
#'   0.0, 0.7
#' ), nrow = 2, byrow = TRUE)   #' diag(1, sqrt(0.7)) -> Var(Λ1_int)=1 (standardized), Var(Λ2_int)=0.49
#' 
#' Llambda_slope <- matrix(c(
#'   0.5, 0.0,
#'   0.1, 0.6
#' ), nrow = 2, byrow = TRUE)
#' #' BLambda_slope = L L^T => will have off-diagonals since l21 != 0
#' 
#' #' Structured method: rho within Λ blocks (intercepts: 0; slopes: 0.2)
#' res_struct <- recover_omega_cov(
#'   Llambda_int   = Llambda_int,
#'   Llambda_slope = Llambda_slope,
#'   mappingL1L2   = mappingL1L2,
#'   method        = "structured",
#'   rho_int       = 0,
#'   rho_slope     = 0.2,
#'   intercept_var_fixed = TRUE,
#'   weights_normalize   = TRUE
#' )
#' 
#' res_struct$BOmega_int   #' 3x3 Ω-intercept covariance
#' res_struct$BOmega_slope #' 3x3 Ω-slope covariance (if provided)
#' 
#' #' Pseudoinverse fallback (matches Λ off-diagonals as much as possible)
#' res_pinv <- recover_omega_cov(
#'   Llambda_int   = Llambda_int,
#'   Llambda_slope = Llambda_slope,
#'   mappingL1L2   = mappingL1L2,
#'   method        = "pseudoinverse",
#'   intercept_var_fixed = TRUE,
#'   weights_normalize   = TRUE
#' )
#' 
#' res_pinv$BOmega_int
#' res_pinv$BOmega_slope

recover_omega_cov <- function(Llambda_int,
                              Llambda_slope = NULL,
                              mappingL1L2,
                              method = c("structured", "pseudoinverse"),
                              rho_int = 0,
                              rho_slope = 0,
                              intercept_var_fixed = TRUE,
                              weights_normalize = TRUE) {
  
  method <- match.arg(method)
  
  nD <- nrow(mappingL1L2)
  nL <- ncol(mappingL1L2)
  
  # -- helper: build BLambda from Cholesky
  chol_to_cov <- function(L) {
    if (is.null(L)) return(NULL)
    S <- L %*% t(L)
    # defensive symmetrization
    S <- (S + t(S)) / 2
    S
  }
  
  BLambda_int   <- chol_to_cov(Llambda_int)
  BLambda_slope <- chol_to_cov(Llambda_slope)
  
  # -- weights (W) for the chosen block (same mapping for both blocks)
  W <- mappingL1L2
  if (weights_normalize) {
    row_sums <- rowSums(abs(W))
    # avoid division by zero
    if (any(row_sums == 0)) stop("A Lambda row has no mapped Omegas (zero weights).")
    W <- W / row_sums
  }
  
  # checks: overlapping Omegas across Lambdas (for structured method we assume disjoint groups)
  omega_usage <- colSums(abs(W) > 0)
  overlapped  <- which(omega_usage > 1)
  
  # -- helper: pinv via SVD (no external packages)
  pinv <- function(A, tol = .Machine$double.eps) {
    sv <- svd(A)
    d <- sv$d
    Dplus <- diag(ifelse(d > tol, 1/d, 0), nrow = length(d))
    sv$v %*% Dplus %*% t(sv$u)
  }
  
  # -- helper: nearest PSD
  nearest_psd <- function(S, eps = 1e-10) {
    S <- (S + t(S)) / 2
    ed <- eigen(S, symmetric = TRUE)
    vals <- pmax(ed$values, eps)
    Spsd <- ed$vectors %*% diag(vals) %*% t(ed$vectors)
    (Spsd + t(Spsd)) / 2
  }
  
  # -- compute BOmega by chosen method
  if (method == "structured") {
    
    if (length(overlapped) > 0) {
      stop("Structured method assumes disjoint Omega groups per Lambda. Found overlapping Omega indices: ",
           paste(overlapped, collapse = ", "),
           ". Use method='pseudoinverse' or redesign mapping.")
    }
    
    # rho can be scalar or per-Lambda vector
    if (length(rho_int) == 1)  rho_int  <- rep(rho_int,  nD)
    if (!is.null(BLambda_slope) && length(rho_slope) == 1) rho_slope <- rep(rho_slope, nD)
    
    # initialize Omega covariance matrices
    BOmega_int   <- matrix(0, nL, nL)
    BOmega_slope <- if (!is.null(BLambda_slope)) matrix(0, nL, nL) else NULL
    
    # indices of Omega per Lambda d (disjoint assumption)
    groups <- lapply(seq_len(nD), function(d) which(abs(W[d, ]) > 0))
    
    # ---- Intercepts block
    for (d in seq_len(nD)) {
      G  <- groups[[d]]
      wd <- W[d, G]
      # Lambda variance for intercepts with optional constraint
      varLambda_d <- if (intercept_var_fixed) 1 else BLambda_int[d, d]
      
      # sums needed: sum of squares, and sum over pairs
      sum_w2   <- sum(wd^2)
      sum_pair <- ((sum(wd))^2 - sum_w2) / 2
      
      denom <- sum_w2 + 2 * rho_int[d] * sum_pair
      if (denom <= 0) stop(sprintf("Non-positive denominator for Lambda %d (intercepts). Check rho_int/weights.", d))
      
      # common variance within Omega block d
      sigma2_d <- varLambda_d / denom
      
      # build block: equal variance & equal correlation
      block <- matrix(rho_int[d] * sigma2_d, nrow = length(G), ncol = length(G))
      diag(block) <- sigma2_d
      
      BOmega_int[G, G] <- block
    }
    
    # ---- Slopes block (if provided)
    if (!is.null(BLambda_slope)) {
      for (d in seq_len(nD)) {
        G  <- groups[[d]]
        wd <- W[d, G]
        
        varLambda_d <- BLambda_slope[d, d]  # slopes typically not constrained
        sum_w2   <- sum(wd^2)
        sum_pair <- ((sum(wd))^2 - sum_w2) / 2
        
        denom <- sum_w2 + 2 * rho_slope[d] * sum_pair
        if (denom <= 0) stop(sprintf("Non-positive denominator for Lambda %d (slopes). Check rho_slope/weights.", d))
        
        sigma2_d <- varLambda_d / denom
        
        block <- matrix(rho_slope[d] * sigma2_d, nrow = length(G), ncol = length(G))
        diag(block) <- sigma2_d
        
        BOmega_slope[G, G] <- block
      }
    }
    
    # Diagnostics: off-diagonals in BLambda must be ~0 under independence across Omega-blocks
    offdiag_int <- BLambda_int - diag(diag(BLambda_int))
    diag(offdiag_int) <- 0
    warn_int <- any(abs(offdiag_int) > 1e-8)
    
    warn_slope <- FALSE
    if (!is.null(BLambda_slope)) {
      offdiag_slope <- BLambda_slope - diag(diag(BLambda_slope))
      diag(offdiag_slope) <- 0
      warn_slope <- any(abs(offdiag_slope) > 1e-8)
    }
    
    out <- list(
      BOmega_int   = BOmega_int,
      BOmega_slope = BOmega_slope,
      W            = W,
      BLambda_int  = BLambda_int,
      BLambda_slope= BLambda_slope,
      method       = "structured",
      diagnostics  = list(
        overlapping_omegas = overlapped,
        BLambda_int_offdiag_nonzero   = warn_int,
        BLambda_slope_offdiag_nonzero = warn_slope
      )
    )
    return(out)
    
  } else {  # method == "pseudoinverse"
    
    # Build Omega covariance via pseudoinverse mapping, then PSD projection
    
    # Optional: enforce intercept variance constraint by replacing BLambda_int diag with 1
    if (intercept_var_fixed) {
      d <- diag(BLambda_int)
      BLambda_int <- BLambda_int
      diag(BLambda_int) <- 1
      # note: off-diagonals are kept; if they are nonzero, they imply cross-Omega correlations in least-squares solution
    }
    
    Wplus <- pinv(W)
    
    BOmega_int   <- nearest_psd(Wplus %*% BLambda_int   %*% t(Wplus))
    BOmega_slope <- if (!is.null(BLambda_slope)) nearest_psd(Wplus %*% BLambda_slope %*% t(Wplus)) else NULL
    
    out <- list(
      BOmega_int   = BOmega_int,
      BOmega_slope = BOmega_slope,
      W            = W,
      BLambda_int  = BLambda_int,
      BLambda_slope= BLambda_slope,
      method       = "pseudoinverse",
      diagnostics  = list(
        rank_W = qr(W)$rank,
        note   = "Solution may not be unique; nearest PSD projection applied."
      )
    )
    return(out)
  }
}
