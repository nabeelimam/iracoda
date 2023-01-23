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
bhattacharyyaC_MN <- function (p, q) sum(sqrt(p * q))

#' Bhattacharyya coefficient for two Dirichlet distributions
#'
#' Calculate the Bhattacharyya coefficient (BC) for two Dirichlet distributions
#' @param p Vector of mean probabilities for Dirichlet distribution 1, should sum to 1.
#' @param q Vector of mean probabilities for Dirichlet distribution 2, should sum to 1.
#' @return The value of BC
#' @examples
#' high_agreement <- bhattacharyyaC_D(c(0.5, 0.4, 0.1), c(0.4, 0.4, 0.2));
#' low_agreement <- bhattacharyyaC_D(c(0.5, 0.4, 0.1), c(0.1, 0.5, 0.4));
#' @export
bhattacharyyaC_D <- function(p, q) {
  # formula from Rauber et al
  tl <- prod(gamma((p/2) + (q/2)))
  bl <- gamma(0.5 * sum(p + q))
  br <- sqrt(prod(gamma(p))) * sqrt(prod(gamma(q)))
  tr <- sqrt(gamma(sum(p)) * gamma(sum(q)))
  return((tl * tr)/(bl * br))

}
