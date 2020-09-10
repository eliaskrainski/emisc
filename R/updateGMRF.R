#' Computes the posterior Gaussian distribution
#' given the (mixed) design and precision matrices.
#' @param y the outcome
#' @param Qe the error precision matrix
#' @param A the design matrix
#' @param Q the prior precision matrix for the
#' latent random (Markov) random field
#' @export
#' @examples
#' ## Consider the Orange dataset
#' ## intercept plus linear trend over age with
#' ## additional trend for each tree j
#' ##  y_ij = a0 + b0*age_ij + b_j age_ij + e_ij
#' y <- Orange$circumference
#' z0 <- model.matrix(~Tree-1, Orange)
#' A <- cbind(1, Orange$age, z0*Orange$age)
#' Qx <- Diagonal(7, c(0,0.001, rep(2843.56, 5)))
#' Qe <- Diagonal(nrow(Orange), 0.010493)
#' up <- updateGMRF(y, Qe, A, Qx)
#' cbind(x=up$mu, sd=sqrt(diag(solve(up$Q))))
#'
#' \dontrun{
#' ## compare with INLA
#'   library(INLA)
#'   ff <- circumference ~ age +
#'      f(Tree, age, model='iid')
#'   res <- inla(ff, data=Orange,
#'      control.inla=list(int.strategy='eb'))
#'   res$summary.fixed[, 1:2]
#'   res$summary.random$Tree[, 1:3]
#'   round(exp(res$mode$theta), 6)
#' }
updateGMRF <- function(y, Qe, A, Qx) {
    aqe <- Matrix::crossprod(A, Qe)
    Qx.new <- Qx + aqe%*%A
    L <- chol(Qx.new)
    x <- Matrix::solve(Qx.new, aqe%*%y)
    return(list(mu=x, Q=Qx.new, L=L,
                sldL=sum(log(diag(L)))))
}
