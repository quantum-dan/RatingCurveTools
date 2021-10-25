# Create rating curves for the restoration scenarios with improved
# setup compared to previous rating curve infrastructure.  Assumes that
# there is only one model geometry (not HF/LF), but its set of functions is capable
# of easily bridging the LFC jump.
#
# Written beginning September 20, 2021 by Daniel Philippus for the Los Angeles
# River Environmental Flows Project at Colorado School of Mines.

library(tidyverse)
library(hydroGOF)
library(readxl)

fit.curves <- function(data,
                       plot = F,
                       save.plot = NA,
                       na.zero = T) {
  # Fit polynomial (+ a few other functions) rating curves to the data.
  # data: ID, x, <ys...> (ID identifies the curve)
  # plot: whether to plot curves
  # na.zero: replace NAs in y with 0
  #
  # returns: tibble(ID,y,nse,rmse,X,X2...X10,LnX,SqrtX,Rt3X) where x<n> means the coefficient
  
  out <- tibble()
  
  ynames <- names(data)[-(1:2)]
  ids <- unique(data$ID)
  
  if (!is.na(save.plot))
    pdf(save.plot)
  
  for (id in ids) {
    data.set <- filter(data, ID == id)
    xs <- data.set$x
    
    for (ix in 1:length(ynames)) {
      yn <- ynames[ix]
      ys <- data.set[[2 + ix]]
      if (na.zero)
        ys <- map_dbl(ys, function(x) if (is.na(x)) 0 else x)
      
      if (sum(!is.na(ys)) > 0) {
        fit <- fit.single(xs, ys)
        print(paste(id, yn, "NSE:", fit$nse, "RMSE:", fit$rmse))
        
        out.row <- tibble(
          ID = id,
          Function = yn,
          nse = fit$nse,
          rmse = fit$rmse,
          X1 = fit$coeffs[[1]],
          X2 = fit$coeffs[[2]],
          X3 = fit$coeffs[[3]],
          X4 = fit$coeffs[[4]],
          X5 = fit$coeffs[[5]],
          X6 = fit$coeffs[[6]],
          X7 = fit$coeffs[[7]],
          X8 = fit$coeffs[[8]],
          X9 = fit$coeffs[[9]],
          X10 = fit$coeffs[[10]],
          LnX = fit$coeffs[[11]],
          SqrtX = fit$coeffs[[12]],
          Rt3X = fit$coeffs[[13]],
          MinQ = min(xs[ys > 0]),
          MaxQ = max(xs)
        )
        
        out <- rbind(out, out.row)
        
        if (plot) {
          sims <- apply.fit(xs, fit$coeffs)
          plt <- tibble(x = xs,
                        Observed = ys,
                        Simulated = sims) %>%
            ggplot(aes(x = xs)) +
            geom_point(aes(y = Observed, color = "Observed")) +
            geom_line(aes(y = Simulated, color = "Simulated")) +
            scale_x_log10() +
            labs(
              x = "X",
              y = yn,
              title = paste(yn, "Observed vs Simulated for", id),
              color = "Series"
            )
          print(plt)
          if (is.na(save.plot) && readline("Enter to continue, Q to quit: ") == "Q")
            stop("User initiated quit")
        }
      } else print(paste("Skipping", id, yn, "due to all-NA data"))
    }
  }
  if (!is.na(save.plot))
    dev.off()
  out
}

apply.fit <- function(x, coeffs) {
  # coeffs: x, x^2, ..., x^10, log(x), sqrt(x), cube rt(x)
  narm <- function(v) if (is.na(v)) 0 else v
  
  coeffs[[1]] * x +
    coeffs[[2]] * x ^ 2 +
    coeffs[[3]] * x ^ 3 +
    coeffs[[4]] * x ^ 4 +
    coeffs[[5]] * x ^ 5 +
    coeffs[[6]] * x ^ 6 +
    coeffs[[7]] * x ^ 7 +
    coeffs[[8]] * x ^ 8 +
    coeffs[[9]] * x ^ 9 +
    coeffs[[10]] * x ^ 10 +
    coeffs[[11]] * log(x) +
    coeffs[[12]] * sqrt(x) +
    coeffs[[13]] * x^(1/3)
}

fit.single <- function(x, y) {
  # Fit y to x for a 10th-order polynomial plus log, sqrt, cube rt
  # returns: list(coeffs, nse, rmse)
  fit <- lm(y ~ 0 +
              x +
              I(x ^ 2) +
              I(x ^ 3) +
              I(x ^ 4) +
              I(x ^ 5) +
              I(x ^ 6) +
              I(x ^ 7) +
              I(x ^ 8) +
              I(x ^ 9) +
              I(x ^ 10) +
              I(log(x)) +
              I(x^(1/2)) +
              I(x^(1/3)))
  coeff <- map_dbl(fit$coefficients, function(x) if (is.na(x)) 0 else x)
  sims <- apply.fit(x, coeff)
  list(
    coeffs = coeff,
    nse = NSE(sims, y),
    rmse = rmse(sims, y)
  )
}