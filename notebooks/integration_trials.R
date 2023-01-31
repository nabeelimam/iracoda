require(calculus)
require(MGLM)

set.seed(321)

D <- 4
n <- 50
size <- 1000

a1 <- rgamma(D, 4)
b1 <- rgamma(D, 14)
Y1 <- rgdirmn(n, size, a1, b1)

a2 <- rgamma(D, 5)
b2 <- rgamma(D, 7)
Y2 <- rgdirmn(n, size, a2, b2)

dist1 <- MGLMfit(Y1, dist = "GDM")
dist2 <- MGLMfit(Y2, dist = "GDM")

sample1 <- Y1[1,]
sample2 <- Y2[1,]

dist_wrapper <- function (dist_pf, Y, ...) {

  params <- unlist(list(...))

  alpha <- params[grep("alpha", names(params))]
  beta <- params[grep("beta", names(params))]

  as.call(list(match.fun(dist_pf), Y, alpha, beta))
}

expr <- dist_wrapper("dgdirmn", Y1[1,],
                     alpha1 = 7.42, alpha2 = 2.29, alpha3 = 2.78, alpha4 = 3.27,
                     beta1 = 12.2, beta2 = 14.4, beta3 = 14.8, beta4 = 11.5)

dens_prod <- function (s1, a1, b1, s2, a2, b2) {
  sqrt(exp(dgdirmn(s1, a1, b1) + dgdirmn(s2, a2, b2)))
}

dens_prod(sample1, a1, b1, sample2, a2, b2)
integral(dens_prod, list(a1 = c(0, 10), a2 = c(0, 10), b1 = c(0, 10), b2 = c(0, 10)
  # a1 = list(c(0, 10), c(0, 10), c(0, 10), c(0, 10)),
  #                         b1 = list(c(0, 10), c(0, 10), c(0, 10), c(0, 10)),
  #                         a2 = list(c(0, 10), c(0, 10), c(0, 10), c(0, 10)),
  #                         b2 = list(c(0, 10), c(0, 10), c(0, 10), c(0, 10))
                         ),
         params = list(s1 = sample1, s2 = sample2)
         )

to_params <- function(l) {
  z <- as.list(l)
  setNames(lapply(names(z), function(x) bquote(args[[.(x)]])), names(z))
}

to_params(a1)
