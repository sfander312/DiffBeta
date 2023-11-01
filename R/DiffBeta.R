# Delta Method CIs for Difference Between Standardized Regression Coefficients

#' Computes delta method-based confidence intervals for the difference
#' between pairs of standardized regression coefficients
#'
#' @description \code{Diff.Beta} returns the lower limit and upper limit of the
#'   delta method based confidence interval for the difference between
#'   standardized regression coefficients. The estimate of the difference itself
#'   and a p-value indicating whether the difference is statistically significant
#'   are also returned. In the case of multiple differences of interest, such as
#'   when there are more than 2 predictor variables, multiple confidence intervals
#'   will be returned. In this case, users can also request intervals corrected for
#'   multiple comparisons using the Bonferroni approach, or p-values corrected
#'   using the Holm or false discovery rate methods. This approach shows better
#'   confidence interval coverage and Type I error rates when compared to an approach
#'   that relies on the standard "textbook" formula to compute the covariance matrix.
#'
#' @details Researchers often compare standardized regression coefficients informally
#'   when comparing predictors in a regression model. However, the informal approach
#'   can be misleading or ineffective.
#'
#'   The approach implemented in \code{Diff.Beta} uses the delta method to
#'   appropriately calculate the standard error for the difference between two
#'   standardized regression coefficients. The function assumes that predictor
#'   variables are random as opposed to fixed, which is typical of most predictor
#'   variables outside of experimental conditions.
#'
#'   The calculation of a confidence interval for the difference in standardized
#'   regression coefficients, \eqn{\beta_{j}-\beta_{k}} requires three elements
#'   from the delta method covariance matrix: the variance of \eqn{\beta_{j}}, the
#'   variance of \eqn{\beta_{k}}, and their covariance. These elements could
#'   also be taken from the standard covariance matrix that duplicates the
#'   formula for unstandardized regression coefficients using z-scores. However,
#'   the standard method has been shown to be biased, producing poorly
#'   performing confidence intervals for single standardized regression
#'   coefficients (see Yuan & Chan, 2011; Jones & Waller, 2013). Intervals
#'   for the difference in standardized coefficients formed using the standard
#'   method can also behave poorly.
#'
#'   In addition to the estimates of \eqn{\beta_{j}-\beta_{k}}, the limits of
#'   the confidence interval, and the p-value, standard output also prints
#'   out the full delta method covariance matrix and the standard error for
#'   each difference.
#'
#'   Regarding function inputs, researchers can enter one of two sets of
#'   information. First, users can enter raw data directly. Predictors
#'   \code{X} are entered separately from the outcome variable \code{y}, and
#'   sample size is derived from the raw data. The relevant covariance
#'   matrices are then computed internally. Alternatively, users can enter
#'   covariance input and sample size. Specifically, the covariance matrix among
#'   predictors (\code{cov.x}), the covariances among the predictors and outcome
#'   (\code{cov.xy}), and the variance of the outcome (\code{var.y}) are entered,
#'   along with the sample size (\code{Nobs}).
#'
#'   When calculating confidence intervals for multiple pairs of \eqn{\beta}
#'   coefficients, users have several options for familywise error adjustment.
#'   First, users can request unadjusted intervals using \code{MCP="no"}. Second
#'   users can request a Bonferroni adjustment using \code{MCP="bonf"}. When using
#'   Bonferroni, users can also set \code{no.diffs} equal to the number of
#'   differences of interest when fewer than all possible differences are desired.
#'   Bonferroni produces "simultaneous" intervals, in which \eqn{1-\alpha}% will
#'   contain the population values of all differences of interest. Appropriate
#'   interval formation is not possible with the Holm or FDR methods, but users
#'   can request these methods to calculate p-values that correct for multiple
#'   tests, by entering \code{MCP="holm"} and \code{MCP="fdr"}, respectively.
#'   When only one pair of coefficients is to be compared, users should specify
#'   \code{MCP="no"}.
#'
#'   When researchers have flexibility in sample size planning, sequential accuracy
#'   in parameter estimation (AIPE) can be conducted to ensure the interval width
#'   is as narrow as desired. Users can set \code{aipe="yes"} for the function to
#'   only print the width of the interval. Otherwise, set \code{aipe="no"}.
#'
#' @param X Raw data of predictor variables
#'   size for a planned study
#' @param y Raw data of outcome variable
#' @param cov.x Covariance matrix among predictor variables
#' @param cov.xy Covariance vector between predictor variables and outcomes
#' @param var.y Variance of the outcome variable
#' @param Nobs Total sample size
#' @param alpha Alpha-level (\eqn{\alpha}) assumed for the planned study
#' @param MCP Type of multiple comparison adjustment. Enter "no" for no adjustment
#' or if only one pair of coefficients of interest. Enter "bonf" for Bonferroni,
#' "holm" for Holm, or "fdr" for the false discovery rate adjustment. The latter
#' two options only adjust p-value, not confidence interval.
#' @param no.diffs Number of pairs of coefficients to be tested. Only used with
#' Bonferroni method. If this value is missing, it is assumed all pairs are tested.
#' @param aipe Whether sequential accuracy in parameter estimation is desired for
#' sample size planning. Enter "yes" if aipe is desired and only interval width
#' will be printed. Otherwise, enter "no".
#' @param digits Number of digits to round confidence limits and p-values to.
#'
#' @return Delta method confidence intervals and p-values for \eqn{\beta_{i}-\beta_{j}}
#' Delta method covariance matrix
#' Delta method standard error
#' Confidence interval width (if \code{aipe="yes"})
#'
#' @export
#' @import stats
#'
#' @examples
#' library(alr4)
#' d <- Rateprof
#' Diff.Beta(X = d[,9:11], y =d$quality, MCP = "bonf", aipe = "no", digits = 3)
#'
#' @author Samantha F. Anderson \email{samantha.f.anderson@asu.edu}
#'
#' @references
#' Anderson, S. F. (2023). A Confidence Interval for the Difference Between
#' Standardized Regression Coefficients.
#'
#' Jones, J. A., & Waller, N. G. (2013). Computing Confidence Intervals for
#' Standardized Regression Coefficients. Psychological Methods, 18, 453.
#'
#' Yuan, K. H., & Chan, W. (2011). Biases and Standard Errors of Standardized
#' Regression Coefficients. Psychometrika 76, 670–690.
#'

Diff.Beta <- function(X = NULL, y = NULL,
                      cov.x = NULL, cov.xy = NULL,
                      var.y = NULL, Nobs = NULL,
                      alpha = .05, MCP = c("no", "bonf", "holm", "fdr"),
                      no.diffs = NULL, aipe = c("no", "yes"),
                      digits = 3) {

  # rounding function
  round_df <- function(df, digits) {
    nums <- vapply(df, is.numeric, FUN.VALUE = logical(1))
    df[,nums] <- round(df[,nums], digits = digits)
    (df)
  }

  # vech function
  vech <- function(x) t(x[!upper.tri(x)])

  # transition/duplicator matrix function
  Dn <- function(x) {
    mat <- diag(x)
    index <- seq(x * (x+1) / 2)
    mat[lower.tri(mat, TRUE)] <- index
    mat[upper.tri(mat)] <- t(mat)[upper.tri(mat)]
    outer(c(mat), index, function(x, y) ifelse(x ==y, 1, 0))
  }

  # Error checking
  if(alpha > 1 | alpha <= 0)
    stop("\n Please specify alpha between 0 and 1 (the default is .05) \n")
  if(missing(aipe))
    stop("\n You must specify aipe = 'yes' or aipe = 'no' \n")
  if(aipe!="no" & aipe!="yes")
    stop("\n You must specify aipe = 'yes' or aipe = 'no' \n")
  if(missing(MCP))
    stop("\n You must specify mcp = 'no', 'bonf', 'holm', or 'fdr' \n")
  if(MCP!="no" & MCP!="bonf" & MCP!="holm" & MCP!="fdr")
    stop("\n You must specify mcp = 'no', 'bonf', 'holm', or 'fdr' \n")

  if(is.null(X) & !is.null(y))
    stop("\n y is not defined\n Need to specify both X and y \n")
  if(!is.null(X) & is.null(y))
    stop("\n X is not defined\n Need to specify both X and y \n")
  if(is.null(X) & is.null(y)) {
    if(is.null(cov.x) | is.null(cov.xy) | is.null(var.y) | is.null(Nobs))
      stop("\n You need to specify all covariances and sample size \n")

    scov <- rbind(cbind(cov.x, cov.xy), c(cov.xy, var.y))
    N <- Nobs
    p <- nrow(cov.x)
  }
  else {
    scov <- cov(cbind(X, y))
    N <- length(y)
    p <- ncol(X)
  }

  # Create covariance matrix of covariances assuming normality
  Kp.lft <- solve(t(Dn(p + 1)) %*% Dn(p + 1)) %*% t(Dn(p + 1)) # D_p+ eqn: p. 675 Yuan & Chan 2011
  cov.cov <- 2 * Kp.lft %*% (scov %x% scov) %*% t(Kp.lft)  # Q: Eqn 17 Yuan & Chan 2011 (p. 675)
  param <- c(vech(scov))
  ncovs <- length(param)

  # Find the vector element numbers for the variances of predictors
  v.x.pl <- c(1, rep(0, p - 1))
  for(i in 2:p) v.x.pl[i] <- v.x.pl[i - 1] + p - (i - 2)

  # Store covariances and variances to use in calculating derivatives
  cx <- scov[1:p, 1:p] # Sigma2_xx submatrix of Sigma_v
  cxy <- scov[1:p, p+1] # sigma2_xy vector of Sigma_v
  vy <- scov[p+1, p+1] # sigma2_y element of Sigma_v
  sx <- sqrt(diag(cx)) # sigma_xx submatrix
  sy <- sqrt(vy) # sigma_y element
  bu <- solve(cx) %*% cxy # unstandardized b coefficient vector
  ncx <- length(vech(cx))

  # Derivatives of standardized regression coefficients w.r.t the covariances
  # Yuan & Chan's (2011) Equation 13 (Jacobian)
  db <- matrix(0, p, ncovs)
  V <- matrix(0, p, ncx)
  V[as.matrix(cbind(1:p, v.x.pl))] <- 1

  # h1(sigma)
  db[, 1:ncx] <- ( diag(c(solve(diag(2 * sx * sy)) %*% bu)) %*% V -
                     diag(sx / sy) %*% (t(bu) %x% solve(cx)) %*% Dn(p) ) # error in JW code? (changed to kronecker prod)
  # h2(sigma)
  db[, (ncx+1):(ncx+p)] <- diag(sx / sy) %*% solve(cx)
  # h3(sigma)
  db[, ncovs] <- -diag(sx / (2 * sy^3)) %*% bu

  # Re-order the derivatives
  cx.nms <- matrix(0, p, p)
  cxy.nms <- c(rep(0, p), "Var_y")
  for(i in 1:p) for(j in 1:p) cx.nms[i,j] <- paste("cov_x", i, "x", j, sep='')
  for(i in 1:p) cxy.nms[i] <- paste("cov_x", i, "y", sep='')

  old.ord <- c(vech(cx.nms), cxy.nms)
  new.ord <- vech(rbind(cbind(cx.nms, cxy.nms[1:p]), c(cxy.nms)))
  db <- db[, match(new.ord, old.ord)]

  # Compute covariance matrix of standardized
  # regression coefficients using the Delta Method
  # includes N-3 denominator (Yuan & Chan , p. 676)
  DEL.cmat <- db %*% cov.cov %*% t(db) / (N - 3)
  b.nms <- NULL

  for(i in 1:p) b.nms[i] <- paste("beta_", i, sep='')
  rownames(DEL.cmat) <- colnames(DEL.cmat) <- b.nms

  # compute differences in coefficients, SEs, and CIs
  beta <- diag(sx) %*% bu * sy^-1
  tc <- qt(alpha / 2, N - p - 1, lower = FALSE) # critical value at selected 1-alpha CI percentage

  # MCP adjustment
  if(MCP=="bonf") {
    if(!is.null(no.diffs)) {
      if(no.diffs > p*(p-1)/2)
        stop("\n no.diffs cannot be more than p(p-1)/2 \n")
      if(no.diffs <= p*(p-1)/2) {no.diffs <- no.diffs}
    }
    if(is.null(no.diffs)) {no.diffs <- (p*(p-1))/2}
    tc <- qt(((alpha/no.diffs) / 2), N - p - 1, lower = FALSE) }

  # create storage matrices for differences
  diff <- matrix(NA,p,p)
  DELvar <- matrix(NA,p,p)
  DELse <- matrix(NA,p,p)

  CI.L <- matrix(NA, p, p)
  CI.U <- matrix(NA, p, p)
  width <- matrix(NA, p, p)
  p.val <- matrix(NA, p, p)
  mat.index <- matrix(NA, p, p)

  # repeat for possible pairs of coefficients
  for(i in 1:p) {
    for(j in 1:p) {

      if(i != j) {

        diff[i,j] <- beta[i]-beta[j] # center of CI
        DELvar[i,j] <- DEL.cmat[i,i] + DEL.cmat[j,j] - (2*DEL.cmat[i,j]) # Delta method variance beta_i-beta_j
        DELse[i,j] <- sqrt(DELvar[i,j]) # Delta method SE beta_i-beta_j

        DELvar[upper.tri(DELvar, diag=TRUE)] <- NA
        DELse[upper.tri(DELse, diag=TRUE)] <- NA # remove duplicates (opposed in sign)
      }}}

  for(i in 1:p) {
    for(j in 1:p) {
      if(i != j) {

        if(!is.na(DELvar[i,j])) {

          CI.L[i,j] <- diff[i,j] - tc * DELse[i,j] # lower limit CI
          CI.U[i,j] <- diff[i,j] + tc * DELse[i,j] # uppder limit CI
          width[i,j] <- CI.U[i,j]-CI.L[i,j] # full CI width
          p.val[i,j] <- 2*pt(q = abs(diff[i,j]/DELse[i,j]), df = N - p - 1, lower.tail=FALSE)
          mat.index[i,j] <- paste(i, ",", j, sep='')
        }
      }
    }}

  # Create output data frames
  CIs <- as.data.frame(cbind(c(mat.index), c(CI.L), c(diff), c(CI.U), c(p.val)))
  DELse <- as.data.frame(cbind(c(mat.index), c(DELse)))
  Width <- as.data.frame(cbind(c(mat.index), c(width)))

  colnames(CIs) <- c("index", "lbound", "estimate", "ubound", "pval")
  colnames(Width) <- c("index", "width")
  colnames(DELse) <- c("index", "se")
  CIs <- na.omit(CIs)
  Width <- na.omit(Width)
  DELse <- na.omit(DELse)

  # corrected p-values for MCP procedures
  if(MCP=="bonf"){
    CIs$pval <- p.adjust(CIs$pval, method = "bonferroni", n = length(CIs$pval))
  }

  if(MCP=="holm"){
    CIs$pval <- p.adjust(CIs$pval, method = "holm", n = length(CIs$pval))
  }

  if(MCP=="fdr"){
    CIs$pval <- p.adjust(CIs$pval, method = "fdr", n = length(CIs$pval))
  }

  # Simplify output data frames
  CIs[, 2:5] <- sapply(CIs[, 2:5], as.numeric)
  CIs <- round_df(CIs, digits=digits)
  Width[, 2] <- sapply(Width[, 2], as.numeric)
  Width <- round_df(Width, digits=digits)

  # Change output depending on whether conducting AIPE sample size planning or not
  if(aipe == "no") { return(list(CIs = CIs, DEL.cmat = DEL.cmat, DELse = DELse)) }
  if(aipe == "yes") { return(list(Width = Width)) }

}
