## Self-contained comparison of the exact svpChange2 AR1 statistic and the
## unified_focus ARP/FOCuS statistic. No figures or output files are created.

if (!requireNamespace("svpChange2", quietly = TRUE))
  stop("Install svpChange2 first")
if (!requireNamespace("focus", quietly = TRUE))
  stop("Install with remotes::install_github('gtromano/unified_focus', subdir='focus')")

simulate_ar1 <- function(n = 1000, tau = NULL, mu1 = 0, mu2 = 0,
                         rho = 0.7, sigma = 1) {
  mu <- if (is.null(tau)) rep(mu1, n) else c(rep(mu1, tau), rep(mu2, n-tau))
  y <- numeric(n)
  y[1] <- mu[1] + stats::rnorm(1, sd = sigma / sqrt(1-rho^2))
  for (i in 2:n)
    y[i] <- mu[i] + rho * (y[i-1] - mu[i-1]) + stats::rnorm(1, sd=sigma)
  y
}

exact_ar1_curve <- function(y, rho = 0.7, sigma2 = 1) {
  n <- length(y); z <- y[-1] - rho * y[-n]
  rss0 <- sum((z - mean(z))^2)
  taus <- 2:(n - 2)
  stat <- vapply(taus, function(tau) {
    design <- matrix(0, nrow=n-1, ncol=2)
    design[seq_len(tau-1),1] <- 1-rho
    design[tau,] <- c(-rho,1)
    if (tau+1 <= n-1) design[(tau+1):(n-1),2] <- 1-rho
    rss1 <- sum(lm.fit(design,z)$residuals^2)
    max(0, (rss0-rss1)/(2*sigma2))
  }, numeric(1))
  list(tau=taus, statistic=stat)
}

## Online exact counterpart of focus_offline(): for each prefix 1:t, return
## the maximum exact AR1 statistic over all candidate changes in that prefix.
exact_ar1_online_curve <- function(y, rho = 0.7, sigma2 = 1) {
  n <- length(y); statistic <- rep(NA_real_, n); changepoint <- rep(NA_integer_, n)
  for (t in seq_len(n)) {
    if (t >= 4L) {
      fit <- svpChange2::AR1_single_change(y[seq_len(t)], gamma=Inf,
                                            rho=rho, sigma2=sigma2)
      ## Put the exact statistic on the 2*log-likelihood-ratio scale.
      statistic[t] <- 2 * fit$statistic
      changepoint[t] <- fit$changepoint
    }
  }
  list(statistic=statistic, changepoint=changepoint)
}

compare_statistics <- function(y, gamma = 8, rho = 0.7, sigma2 = 1,
                               mu0_arp = NULL) {
  exact <- svpChange2::AR1_single_change(
    y, gamma = gamma, rho = rho, sigma2 = sigma2,
    profile_sigma = FALSE
  )
  focus_fit <- if (is.null(mu0_arp)) {
    focus::focus_offline(y, threshold=gamma, type="arp", family="arp", rho=rho)
  } else {
    focus::focus_offline(y, threshold=gamma, type="arp", family="arp",
                         rho=rho, mu0_arp=mu0_arp)
  }
  focus_stat <- as.numeric(focus_fit$stat)
  keep <- is.finite(focus_stat)
  focus_stat <- focus_stat[keep]
  focus_time <- which(keep)[which.max(focus_stat)]
  focus_cp <- as.integer(focus_fit$changepoint[focus_time])
  data.frame(
    exact_statistic = exact$statistic,
    exact_changepoint = exact$changepoint,
    focus_statistic = max(focus_stat, 0),
    focus_changepoint = focus_cp,
    statistic_ratio = max(focus_stat, 0) / exact$statistic,
    same_changepoint = identical(as.integer(exact$changepoint), focus_cp)
  )
}

plot_statistic_panels <- function(seed = 10, n = 1000, rho = 0.7,
                                  sigma = 1, gamma = 1e6) {
  set.seed(seed)
  scenarios <- list(
    no_change = simulate_ar1(n=n, rho=rho, sigma=sigma),
    weak_change = simulate_ar1(n=n, tau=n/2, mu2=1, rho=rho, sigma=sigma),
    strong_change = simulate_ar1(n=n, tau=n/2, mu2=3, rho=rho, sigma=sigma)
  )
  oldpar <- par(mfrow=c(length(scenarios),2), mar=c(3,4,2,1)); on.exit(par(oldpar))
  for (scenario in names(scenarios)) {
    y <- scenarios[[scenario]]; ex <- exact_ar1_online_curve(y,rho,sigma^2)
    fu <- focus::focus_offline(y, threshold=gamma, type="arp",
                               family="arp", rho=rho)
    fs <- as.numeric(fu$stat); fs[!is.finite(fs)] <- NA_real_
    plot(seq_along(ex$statistic), ex$statistic, type="l", xlim=c(1,n),
         main=paste(scenario,"- exact AR1"), xlab="time",
         ylab="statistic"); abline(h=gamma,lty=2)
    abline(v=which.max(ex$statistic),col=2,lty=2)
    plot(seq_along(fs), fs, type="l", xlim=c(1,n),
         main=paste(scenario,"- unified FOCuS"), xlab="time",
         ylab="statistic"); abline(h=gamma,lty=2)
    abline(v=which.max(fs),col=2,lty=2)
  }
  invisible(NULL)
}

compare_curve_differences <- function(seed = 10, n = 1000, rho = 0.7,
                                      sigma = 1, gamma = 1e6) {
  set.seed(seed)
  scenarios <- list(
    no_change = simulate_ar1(n=n, rho=rho, sigma=sigma),
    weak_change = simulate_ar1(n=n, tau=n/2, mu2=1, rho=rho, sigma=sigma),
    strong_change = simulate_ar1(n=n, tau=n/2, mu2=3, rho=rho, sigma=sigma)
  )
  do.call(rbind, lapply(names(scenarios), function(scenario) {
    y <- scenarios[[scenario]]
    exact <- exact_ar1_online_curve(y, rho, sigma^2)$statistic
    unified <- as.numeric(focus::focus_offline(
      y, threshold=gamma, type="arp", family="arp", rho=rho
    )$stat)
    unified <- c(unified, rep(NA_real_, max(0, length(exact)-length(unified))))
    m <- min(length(exact), length(unified))
    ok <- is.finite(exact[seq_len(m)]) & is.finite(unified[seq_len(m)])
    differences <- abs(exact[seq_len(m)][ok] - unified[seq_len(m)][ok])
    data.frame(scenario=scenario,
               max_absolute_difference=if (length(differences)) max(differences) else NA_real_,
               rmse=if (length(differences)) sqrt(mean(differences^2)) else NA_real_,
               stringsAsFactors=FALSE)
  }))
}

## To compute differences for all three experiments:
## curve_differences <- compare_curve_differences()
## print(curve_differences)

run_comparison <- function(seed = 10, n = 1000, rho = 0.7, sigma = 1,
                           gamma = 8) {
  set.seed(seed)
  scenarios <- list(
    no_change = simulate_ar1(n=n, rho=rho, sigma=sigma),
    weak_change = simulate_ar1(n=n, tau=n/2, mu2=1, rho=rho, sigma=sigma),
    strong_change = simulate_ar1(n=n, tau=n/2, mu2=3, rho=rho, sigma=sigma)
  )
  result <- do.call(rbind, lapply(names(scenarios), function(s) {
    cbind(scenario=s, compare_statistics(scenarios[[s]], gamma, rho, sigma^2))
  }))
  rownames(result) <- NULL
  result
}

## Run explicitly:
print(run_comparison())
plot_statistic_panels()
print(compare_curve_differences())
