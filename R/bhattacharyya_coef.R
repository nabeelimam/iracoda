#' Bhattacharyya coefficient for two multinomial distributions
#'
#' Calculate the Bhattacharyya coefficient (BC) for two multinomial distributions
#' @param p Vector of probabilities for multinomial distribution 1, should sum to 1.
#' @param q Vector of probabilities for multinomial distribution 2, should sum to 1.
#' @return The value of BC
#' @examples
#' high_agreement <- bhattacharyyaC_MN(c(0.5, 0.4, 0.1), c(0.4, 0.4, 0.2));
#' low_agreement <- bhattacharyyaC_MN(c(0.5, 0.4, 0.1), c(0.1, 0.5, 0.4));
#' @export
bhattacharyyaC_MN <- function (p, q) {

  if (sum(p) > 1 || sum(q) > 1)
    stop("One or both probability vectors do not add up to 1")

  sum(sqrt(p * q))
}

#' Bhattacharyya coefficient for two Dirichlet distributions
#'
#' Calculate the Bhattacharyya coefficient (BC) for two Dirichlet distributions
#' @param a1 Vector of alpha parameters for Dirichlet distribution 1. All alpha must be > 0.
#' @param a2 Vector of alpha parameters for Dirichlet distribution 2. All alpha must be > 0.
#' @return The value of BC
#' @examples
#' high_agreement <- bhattacharyyaC_D(c(5, 4, 1), c(4, 4, 2));
#' low_agreement <- bhattacharyyaC_D(c(5, 4, 1), c(1, 5, 4));
#' @export
bhattacharyyaC_D <- function(a1, a2) {

  if (any(a1 <= 0) || any(a2 <= 0))
    stop("All Dirichlet paramters should be > 0")

  # formula from Rauber et al. 2008 (https://ieeexplore.ieee.org/document/4604388)

  exp(
    sum(lgamma(0.5 * (a1 + a2))) +
      (0.5 * (lgamma(sum(a1)) + lgamma(sum(a2)))) -
      lgamma(0.5 * sum(a1 + a2)) -
         (0.5 * (sum(lgamma(a1)) + sum(lgamma(a2))))
  )
}

