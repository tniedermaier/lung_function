###########################################################################
# COPD Markov progression parameters from longitudinal FEV1 data
# Pure base R: no tidyverse required.
#
# Goal:
#   Estimate annual transition probabilities for a 1-year Markov cycle:
#     prog1 = GOLD 1 -> GOLD 2
#     prog2 = GOLD 2 -> GOLD 3
#     prog3 = GOLD 3 -> GOLD 4
#
# Main idea:
#   Use the observed pairs of consecutive visits. For every interval, determine
#   the GOLD stage at the beginning and whether the patient crossed the next
#   GOLD threshold by the next visit.
#
#   Because visit intervals are not exactly 1.5 years, estimate a constant
#   transition intensity lambda per year:
#       P(transition over dt years) = 1 - exp(-lambda * dt)
#   Then convert to the annual Markov-cycle probability:
#       prog = 1 - exp(-lambda)
#
# Uncertainty:
#   Non-parametric bootstrap by patient, not by row/interval, to preserve the
#   within-patient dependence between repeated intervals.
###########################################################################

rm(list = ls(all.names = TRUE))
set.seed(12345)

## ---- User settings ------------------------------------------------------
file_name <- "dat_long.csv"
B <- 10000                    # bootstrap replicates; reduce to e.g. 1000 for a quick test

id_col    <- "id"
time_col  <- "time_years"
fev1_col  <- "fev1_percent_predicted"
visit_col <- "visit"

## ---- Read data ----------------------------------------------------------
dat <- read.csv2(file_name, stringsAsFactors = FALSE, check.names = FALSE)

# If the CSV contains an exported row-number column with an empty or unnamed
# column name, remove it. This is cosmetic, but avoids confusion.
if (names(dat)[1] == "" || grepl("^Unnamed", names(dat)[1])) {
  dat <- dat[, -1, drop = FALSE]
}

# Ensure numeric columns are numeric. read.csv2 should already handle decimal
# commas, but this makes the script more robust if columns were imported as text.
dat[[time_col]] <- as.numeric(dat[[time_col]])
dat[[fev1_col]] <- as.numeric(dat[[fev1_col]])
dat[[visit_col]] <- as.numeric(dat[[visit_col]])

# Keep only rows with the minimum required information.
ok <- !is.na(dat[[id_col]]) & is.finite(dat[[time_col]]) & is.finite(dat[[fev1_col]])
dat <- dat[ok, , drop = FALSE]

## ---- GOLD classification ------------------------------------------------
make_gold <- function(fev1) {
  gold <- rep(NA_integer_, length(fev1))
  gold[fev1 > 80]              <- 1L
  gold[fev1 > 50 & fev1 <= 80] <- 2L
  gold[fev1 > 30 & fev1 <= 50] <- 3L
  gold[fev1 <= 30]             <- 4L
  gold
}

## ---- Build consecutive visit intervals --------------------------------
make_intervals <- function(d) {
  # Sort within patient by time, then visit number.
  ord <- order(d[[id_col]], d[[time_col]], d[[visit_col]])
  d <- d[ord, , drop = FALSE]
  row.names(d) <- NULL

  d$gold_start <- make_gold(d[[fev1_col]])

  # Next-observation information within each patient.
  d$next_id   <- c(d[[id_col]][-1], NA)
  d$next_time <- c(d[[time_col]][-1], NA)
  d$next_fev1 <- c(d[[fev1_col]][-1], NA)

  same_patient <- !is.na(d$next_id) & d[[id_col]] == d$next_id
  d$dt <- d$next_time - d[[time_col]]

  intervals <- d[same_patient & is.finite(d$dt) & d$dt > 0 &
                   is.finite(d$next_fev1) & !is.na(d$gold_start),
                 c(id_col, time_col, visit_col, fev1_col,
                   "gold_start", "next_time", "next_fev1", "dt"),
                 drop = FALSE]

  row.names(intervals) <- NULL
  intervals
}

intervals <- make_intervals(dat)

## ---- Estimate one transition ------------------------------------------
# Estimate lambda and the corresponding 1-year probability.
# Weights are patient-bootstrap weights mapped to interval rows.
estimate_one_transition <- function(intervals, start_gold, threshold, weights = NULL) {
  ii <- which(intervals$gold_start == start_gold)

  if (length(ii) == 0L) {
    return(c(prob_1y = NA_real_, lambda = NA_real_,
             n_intervals = 0, n_events = 0))
  }

  dt <- intervals$dt[ii]
  y  <- as.integer(intervals$next_fev1[ii] <= threshold)

  if (is.null(weights)) {
    w <- rep(1, length(ii))
  } else {
    w <- weights[ii]
  }

  keep <- is.finite(dt) & dt > 0 & is.finite(y) & is.finite(w) & w > 0
  dt <- dt[keep]
  y  <- y[keep]
  w  <- w[keep]

  n_intervals <- sum(w)
  n_events <- sum(w * y)

  if (n_intervals == 0) {
    return(c(prob_1y = NA_real_, lambda = NA_real_,
             n_intervals = 0, n_events = 0))
  }

  if (n_events == 0) {
    return(c(prob_1y = 0, lambda = 0,
             n_intervals = n_intervals, n_events = n_events))
  }

  # Negative weighted log likelihood:
  # y_i ~ Bernoulli(1 - exp(-lambda * dt_i))
  neg_loglik <- function(lambda) {
    p <- 1 - exp(-lambda * dt)
    p <- pmin(pmax(p, 1e-12), 1 - 1e-12)
    -sum(w * (y * log(p) + (1 - y) * log1p(-p)))
  }

  # A very high upper bound is unnecessary here. Lambda=10/year means an annual
  # probability of 0.99995. If your data are extreme, increase upper_bound.
  opt <- optimize(neg_loglik, interval = c(1e-12, 10))
  lambda_hat <- opt$minimum
  prob_1y <- 1 - exp(-lambda_hat)

  c(prob_1y = prob_1y, lambda = lambda_hat,
    n_intervals = n_intervals, n_events = n_events)
}

## ---- Estimate all three progression probabilities ----------------------
estimate_all <- function(intervals, interval_weights = NULL) {
  r1 <- estimate_one_transition(intervals, start_gold = 1L, threshold = 80,
                                weights = interval_weights)
  r2 <- estimate_one_transition(intervals, start_gold = 2L, threshold = 50,
                                weights = interval_weights)
  r3 <- estimate_one_transition(intervals, start_gold = 3L, threshold = 30,
                                weights = interval_weights)

  data.frame(
    parameter = c("prog1", "prog2", "prog3"),
    transition = c("GOLD 1 -> GOLD 2", "GOLD 2 -> GOLD 3", "GOLD 3 -> GOLD 4"),
    estimate = c(r1["prob_1y"], r2["prob_1y"], r3["prob_1y"]),
    lambda_per_year = c(r1["lambda"], r2["lambda"], r3["lambda"]),
    n_intervals = c(r1["n_intervals"], r2["n_intervals"], r3["n_intervals"]),
    n_events = c(r1["n_events"], r2["n_events"], r3["n_events"]),
    row.names = NULL,
    check.names = FALSE
  )
}

point <- estimate_all(intervals)

## ---- Patient-level bootstrap ------------------------------------------
bootstrap_by_patient <- function(intervals, B = 10000) {
  ids <- unique(intervals[[id_col]])
  id_index <- match(intervals[[id_col]], ids)
  n_ids <- length(ids)

  boot <- matrix(NA_real_, nrow = B, ncol = 3)
  colnames(boot) <- c("prog1", "prog2", "prog3")

  for (b in seq_len(B)) {
    sampled <- sample.int(n_ids, size = n_ids, replace = TRUE)

    # Patient weights: how often each patient was sampled.
    patient_weights <- tabulate(sampled, nbins = n_ids)

    # Map patient weights to interval rows.
    interval_weights <- patient_weights[id_index]

    est_b <- estimate_all(intervals, interval_weights = interval_weights)
    boot[b, ] <- est_b$estimate
  }

  boot
}

boot <- bootstrap_by_patient(intervals, B = B)

ci <- t(apply(boot, 2, function(z) {
  quantile(z, probs = c(0.025, 0.50, 0.975), na.rm = TRUE, names = FALSE)
}))
colnames(ci) <- c("ci_low", "boot_median", "ci_high")

out <- cbind(point, ci, boot_sd = apply(boot, 2, sd, na.rm = TRUE))

## ---- Final parameters ---------------------------------------------------
prog1 <- out$estimate[out$parameter == "prog1"]
prog2 <- out$estimate[out$parameter == "prog2"]
prog3 <- out$estimate[out$parameter == "prog3"]

print(out, digits = 4)
cat("\nMarkov parameters for 1-year cycle:\n")
print(c(prog1 = prog1, prog2 = prog2, prog3 = prog3), digits = 6)

write.csv(out, "copd_progression_from_long_results.csv", row.names = FALSE)
saveRDS(list(point = point, summary = out, bootstrap = boot),
        "copd_progression_from_long_bootstrap.rds")

###########################################################################
# Optional diagnostic: simple observed interval proportions before annualising
###########################################################################
cat("\nObserved crude interval proportions, before annualisation:\n")
for (k in 1:3) {
  thr <- c("1" = 80, "2" = 50, "3" = 30)[as.character(k)]
  tmp <- intervals[intervals$gold_start == k, , drop = FALSE]
  ev <- tmp$next_fev1 <= as.numeric(thr)
  cat("\nGOLD", k, ": intervals =", nrow(tmp),
      ", events =", sum(ev),
      ", crude interval proportion =", round(mean(ev), 4),
      ", mean dt =", round(mean(tmp$dt), 3), "years")
}
cat("\n")
