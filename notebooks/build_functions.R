require(calculus)
require(MGLM)
require(microbenchmark)

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

#### Example ####

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

ff

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

flatten_params <- function (param_list, dindex) {

  params_flat <- list()
  param_names <- c("Y", "alpha", "beta")

  for (p in 1:length(param_list)) {

    ls <- setNames(param_list[[p]],
                   paste0(param_names[p], dindex, "_", 1:length(param_list[[p]])))

    params_flat <- append(params_flat, ls)

  }

  return(params_flat)
}

params <- list(Y1[1,], a1, b1)
flatten_params(params, 1)

make_args <- function(x) {
  lapply(names(x), function(y) bquote(args[[.(y)]]))
}

make_args(flatten_params(params, 1))

dist_wrapper <- function (dist_pf, ...) {

  params <- list(...)

  # alpha <- params[grep(paste0("alpha", dindex), names(params))]
  # beta <- params[grep(paste0("beta", dindex), names(params))]

  as.call(list(match.fun(dist_pf), ...))
}

dgd <- function (args) {

  Y <- args[grep("Y", names(args))]
  a <- args[grep("alpha", names(args))]
  b <- args[grep("beta", names(args))]

  m <- Reduce(sum, Y)
  z <- rev(Y)

  for (j in 2:length(z)) {
    z[[j]] <- z[[j]] + z[[j - 1]]
  }
  z <- rev(z)
  logC <- lgamma(m + 1) -
    Reduce(sum, Map(function (x) lgamma(x + 1), Y))

  l <- logC

  for (i in 1:length(a)) {

    l <- l + (
      lgamma(a[[i]] + Y[[i]]) + lgamma(b[[i]] + z[[i+1]]) + lgamma(a[[i]] + b[[i]]) -
        (lgamma(a[[i]]) + lgamma(b[[i]]) + lgamma(a[[i]] + b[[i]] + z[[i]]))
    )
  }

  return(l)
}

microbenchmark(
  dgd(flatten_params(
    list(Y1[1,], 1:4, 4:1), 1)),
  dgdirmn(Y1[1,], 1:4, 4:1),
  times = 100,
  check = "equal"
)


dist_wrapper("dgd", flatten_params(params, 1))

generate_fun <- function(dist_fun, params, dindices) {

  dist_expr <- Map(function(d, p, i) {
    dist_wrapper(d, flatten_params(p, i))
  }, dist_fun, params, dindices)

  fb <- Reduce(multroot, dist_expr)
  f <- function(args, ...) {}
  body(f) <- fb
  f
}

integral_fun <- generate_fun(
  dist_fun = "dgd",
  params = list(list(Y1[1,],  a1, b1), list(Y2[1,], a2, b2)),
  dindices = 1:2)

eval(integral_fun(NULL))

# bounds <- append(
#   rep(list(c(0, Inf)), length(a1) * 4),
#   rep(list(1), length(Y1[1,]) * 2)
# )

bounds <- rep(list(c(0, 10)), length(a1) * 4)
bound_names <- rbind.data.frame(
  expand.grid(
    c("alpha", "beta"),
    1:2,
    1:length(a1))
  # ,
  # expand.grid(
  #   "Y",
  #   1:2,
  #   1:length(Y1[1,])
  # )
)
bound_names <- apply(bound_names, 1,
                     function (x) paste0(x[1], x[2], "_", x[3]))
names(bounds) <- bound_names

integral(integral_fun, bounds, verbose = T)

