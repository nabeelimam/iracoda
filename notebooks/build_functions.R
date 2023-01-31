require(calculus)
require(MGLM)

set.seed(321)

D <- 4
n <- 5
size <- 1000

a1 <- rgamma(D, 4)
b1 <- rgamma(D, 14)
Y1 <- rgdirmn(n, size, a1, b1)

a2 <- rgamma(D, 5)
b2 <- rgamma(D, 7)
Y2 <- rgdirmn(n, size, a2, b2)

# dist1 <- MGLMfit(Y1, dist = "GDM")
# dist2 <- MGLMfit(Y2, dist = "GDM")

sample1 <- Y1[1,]
sample2 <- Y2[1,]

to_params <- function(l) {
  z <- as.list(l)
  setNames(lapply(names(z), function(x) bquote(args[[.(x)]])), names(z))
}

add_exprs <- function(...) {
  x <- list(...)
  Reduce(function(a,b) bquote(.(a) + .(b)), x)
}

get_densities <- function(f) {
  lapply(paste0("d", f), as.name)
}

weight_expr <- function(w, e) {
  bquote(.(w) * .(e))
}
add_params <- function(x, p) {
  as.call(c(as.list(x), p))
}
call_with_x <- function(fn) {
  as.call(list(fn, quote(x)))
}

fitmix <- function(data, dist, params, weights) {
  fb <- Reduce( add_exprs, Map(function(d, p, w) {
    weight_expr(w, add_params(call_with_x(d), to_params(p)))
  }, get_densities(dist), params, weights))
  f <- function(x, args) {}
  body(f) <- fb
  f
}

add_exprs(1, 2)

Map(function(d, p, w) {
  weight_expr(w, add_params(call_with_x(d), to_params(p)))
}, get_densities(c("norm","chisq")), list(c(mean=0,sd=3),c(df=2)), c(0.5,0.5))

call_with_x("dnorm")

to_params(c(mean=0,sd=3))

add_params(call_with_x("dnorm"), to_params(c(mean=0,sd=3)))

  #### Tests ####
ff <- fitmix(data, dist=c("norm","chisq"), params=list(c(mean=0,sd=3),c(df=2)),
             weights=c(0.5,0.5))

ff(0, list(mean=3, sd=2, df=2))

0.5 * dnorm(0, mean = 3, sd = 2) + 0.5 * dchisq(0, df = 2)

#### Build ####

multroot <- function(...) {
  x <- list(...)
  n <- length(x)
  logsum <- Reduce(function(a,b) bquote(.(a) + .(b)), x)
  return(bquote(exp(.(logsum)/.(n))))
}

# eval(multroot(1, 2, 3, 4))

flatten_params <- function (param_list) {

  params_flat <- list()

  for (p in 1:length(params)) {

    alpha_names <- setNames(params[[p]][[1]],
                            paste0("alpha", p, "_", 1:length(params[[p]][[1]])))

    beta_names <- setNames(params[[p]][[2]],
                           paste0("beta", p, "_", 1:length(params[[p]][[2]])))
    params_flat <- append(params_flat, alpha_names)
    params_flat <- append(params_flat, beta_names)

  }

  return(params_flat)
}

params <- list(list(a1, b1), list(a2, b2))
flatten_params(params)

make_args <- function(x) {
  setNames(lapply(names(x), function(y) bquote(args[[.(y)]])), names(x))
}

all_args <- make_args(flatten_params(params))

dist_wrapper <- function (dist_pf, Y, dindex, ...) {

  params <- unlist(list(...))

  alpha <- params[grep(paste0("alpha", dindex), names(params))]
  beta <- params[grep(paste0("beta", dindex), names(params))]

  as.call(list(match.fun(dist_pf), Y, alpha, beta))
}

dgd <- function (Y, a, b) {

  m <- sum(Y)
  z <- rev(cumsum(rev(Y)))
  logC <- lgamma(m + 1) - sum(lgamma(Y + 1))

  l <- logC

  for (i in 1:length(a)) {
    l <- l + (
      lgamma(a[[i]] + Y[i]) + lgamma(b[[i]] + z[i+1]) + lgamma(a[[i]] + b[[i]]) -
        (lgamma(a[[i]]) + lgamma(b[[i]]) + lgamma(a[[i]] + b[[i]] + z[[i]]))
    )
  }

  return(l)
}

microbenchmark(
  dgd(Y1[1,], list(1, 2, 3, 4), list(4, 3, 2, 1)),
  dgdirmn(Y1[1,], c(1, 2, 3, 4), c(4, 3, 2, 1)),
  times = 100,
  check = "equal"
)


dist_wrapper("dgd", Y1[1,], dindex = 1, all_args)

integral_fun <- function(data, dist_fun, dindices, params) {

  fb <- Reduce(multroot, Map(function(d, s, i, p) {
    dist_wrapper(d, s, i, make_args(flatten_params(p)))
  }, dist_fun, data, dindices, params))
  f <- function(Y, args) {}
  body(f) <- fb
  f
}

integral_fun(
  data = list(Y1[1,], Y2[1,]),
  dist_fun = rep("dgd", 2),
  dindices = 1:2,
  params = list(c(a1, b1), c(a2, b2)))
