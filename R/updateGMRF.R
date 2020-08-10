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
#' ## one tred for each tree
#' ##   y_{ij} = a + b_j age_i + error_{ij}
#' y <- Orange$circumference
#' z0 <- model.matrix(~Tree-1, Orange)
#' A <- cbind(1, Orange$age, z0*Orange$age)
#' Qx <- Diagonal(2 + 5, c(0, 0, rep(3000, 5)))
#' Qe <- Diagonal(nrow(Orange), 0.01)
#' result <- updateGMRF(y, Qe, A, Qx)
#' print(cbind(x=result$mu,
#'    sd=sqrt(diag(chol2inv(result$L)))))
#'
#' ## compare with INLA result
#' \dontrun{
#'   ff <- circumference ~ f(Tree, age, model='iid')
#'   res <- inla(ff, data=Orange)
#'   res$summary.fixed[, 1:2]
#'   res$summary.random$Tree[, 1:3]
#' }
updateGMRF <- function(y, Qe, A, Qx) {
    aqe <- Matrix::crossprod(A, Qe)
    a2qe <- aqe%*%A
    Qx.new <- Qx + a2qe
    L <- chol(Qx.new)
    W <- aqe%*%y
    x <- Matrix::solve(
        L, Matrix::solve(t(L), W))
    return(list(mu=x, Q=Qx.new, L=L,
                sldL=sum(log(diag(L)))))
}
