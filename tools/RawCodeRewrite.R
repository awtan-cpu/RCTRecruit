# 1. Prepare data  #############################################################
grips <- readr::read_csv(file = "tools/data/GRIPS_log_by_day_with_gaps.csv")
grips$Screening_Date  <- as.Date(grips$Screening_Date, format = "%m/%d/%Y")
grips$week            <- lubridate::isoweek(grips$Screening_Date)
grips$year            <- lubridate::isoyear(grips$Screening_Date)
grips$isholiday       <- tis::isHoliday(grips$Screening_Date)
grips$notactiveenrolm <- (grips$`N Rows` == 1)
gripsweek <- aggregate(
  cbind(
    Sum_Enrolled,
    Sum_Meets_Crit,
    Sum_Prescreened,
    notactiveenrolm,
    isholiday) ~ week + year, data = grips, FUN = sum
)
grips1sty <- gripsweek[2:53, ]
grips2ndy <- gripsweek[54:105, ]
gripsy2   <- gripsweek[2:105, ]
gripsr    <- gripsweek[106:143, ]

# 2. Prepare sampling weight distributions #####################################
indxWt <- \(t) ((0L:51L + 26L - t) %% 52L) + 1L
binomWt <- stats::dbinom(0L:51L, 51L, 0.5) |> (\(x) x / sum(x))()
cauchyWt <- seq(-3L, 3L, 6 / 51) |> stats::dcauchy() |> (\(x) x / sum(x))()
probwts <- list()
probwts[["binomial"]] <- lapply(1L:52L, \(t) binomWt[indxWt(t)])
probwts[["cauchy"]]   <- lapply(1L:52L, \(t) cauchyWt[indxWt(t)])
probwts[["uniform"]]  <- lapply(1L:52L, \(t) rep(1/52, 52L))

# 3. Define helper functions to handle gap weeks and holiday weeks         #####
fillEmptyWk <- \(train, probs, EmptyWkID) {
  gapIdx <- which(EmptyWkID)
  out <- train
  for (i in gapIdx) {
    wt <- probs[[i]]
    wt[gapIdx] <- 0
    wt <- wt / sum(wt)
    out[i] <- sum(train * wt)
  }
  out
}

ignoreHolidays <- \(probs, holidayID) {
  includedWeeks <- 1 * !holidayID
  lapply(probs, \(x) (\(y) y / sum(y))(x * includedWeeks))
}

# 4. The three main functions                                              #####
wk2NSub <- \(train, probs, EmptyWkID = NULL, holidayID = NULL, efCoef = 1, 
             nSim = 1e4L, nSub = 50) {
  if (!is.null(holidayID)) probs <- ignoreHolidays(probs, holidayID)
  if (!is.null(EmptyWkID)) train <- fillEmptyWk(train, probs, EmptyWkID)
  train <- train * efCoef
  y  <- vector("numeric", nSim)
  qt <- c(.025, .5, .975)
  for (i in seq_len(nSim)) {
    n <- 0L
    k <- 0L
    while (n < nSub) {
      p <- probs[[k %% 52L + 1L]]
      n <- n + train[[.Internal(sample(52L, 1L, FALSE, p))]]
      k <- k + 1L
    }
    y[[i]] <-  k;
  }
  list(weeks = y, CI = quantile(y, qt))
}

cumPred <- \(train, probs, EmptyWkID = NULL, holidayID = NULL, efCoef = 1, 
             nSim = 1e4L) {
  if (!is.null(holidayID)) probs <- ignoreHolidays(probs, holidayID)
  if (!is.null(EmptyWkID)) train <- fillEmptyWk(train, probs, EmptyWkID)
  train <- train * efCoef
  y  <- vector("numeric", nSim)
  qt <- c(.025, .5, .975)
  outPutNames <- c("Week", "LCI", "MED", "HCI")
  out <- list()
  for (i in (seq_len(52L) - 1L)) {
    p  <- probs[[i %% 52L + 1L]]
    y  <- y + train[.Internal(sample(52L, nSim, TRUE, p))]
    wk <- i + 1L
    out[[wk]] <- c(wk, quantile(y, qt)) |> setNames(outPutNames)
  }
  do.call(rbind, out) |> as.data.frame() |> round()
}

getDist <- \(train, target, probs, EmptyWkID = NULL, holidayID = NULL, 
             efCoef = 1, nSim = 1e4L) {
  if (!is.null(holidayID)) probs <- ignoreHolidays(probs, holidayID)
  if (!is.null(EmptyWkID)) train <- fillEmptyWk(train, probs, EmptyWkID)
  train <- train * efCoef
  cumTarget <- cumsum(target)
  y  <- vector("numeric", nSim)
  qt <- c(.025, .5, .975)
  for (k in seq_len(nSim)) {
    pred <- vector("numeric", 53L)
    for (j in 2L:53L) {
      i <- j - 1L;
      p <- probs[[i %% 52L + 1L]]
      pred[[j]] <-  pred[[i]] + train[[.Internal(sample(52L, 1L, FALSE, p))]]
    }
    y[[k]] <-  sqrt(sum((pred[-1L] - cumTarget)^2));
  }
  list(EDist = y, CI = quantile(y, qt))
}

# 5. Scenarios simulation                                                  #####
scL <- lapply(1:6, list)
x <- grips1sty$Sum_Enrolled
y <- grips2ndy$Sum_Enrolled

# Scenario 1 - Bootstrap (uniform weights)                                 #####
i <- 1L
wt <- probwts$uniform
scL[[i]][["desc"]]   <- "Bootstrapping"
# Number of weeks to recruit 50 subjects
scL[[i]][["wkTo50"]] <- wk2NSub(x,    wt)
# Predicted median cumulative recruitment with CI over the next 52 weeks
scL[[i]][["pred"]] <- cumPred(x,    wt)
# Euclidean Distance
scL[[i]][["EDist"]]  <- getDist(x, y, wt)

# Scenario 2 - Binomial                                                    #####
i <- 2L
wt <- probwts$binomial
scL[[i]] <- list(
  desc   = "Binomial",
  wkTo50 = wk2NSub(x,    wt),
  pred   = cumPred(x,    wt),
  EDist  = getDist(x, y, wt)
)

# Scenario 3 - Binomial; gaps filled                                       #####
i <- 3L
gapID <- grips1sty$Sum_Prescreened == 0
scL[[i]] <- list(
  desc   = "Binomial; gaps filled",
  wkTo50 = wk2NSub(x,    wt, gapID),
  pred   = cumPred(x,    wt, gapID),
  EDist  = getDist(x, y, wt, gapID)
)

# Scenario 4 - Binomial; holidays accounted for                            #####
i <- 4L
hdays <- grips1sty$isholiday == 1
scL[[i]] <- list(
  desc   = "Binomial; holidays accounted for",
  wkTo50 = wk2NSub(x,    wt, holidayID = hdays),
  pred   = cumPred(x,    wt, holidayID = hdays),
  EDist  = getDist(x, y, wt, holidayID = hdays)
)

# Scenario 5 - Binomial; gaps and holidays                                 #####
i <- 5L
scL[[i]] <- list(
  desc   = "Binomial; gaps and holidays",
  wkTo50 = wk2NSub(x,    wt, gapID, hdays),
  pred   = cumPred(x,    wt, gapID, hdays),
  EDist  = getDist(x, y, wt, gapID, hdays)
)

# Scenario 6 - Cauchy; gaps and holidays                                   #####
i <- 6L
wt <- probwts$cauchy
scL[[i]] <- list(
  desc   = "Cauchy; gaps and holidays",
  wkTo50 = wk2NSub(x,    wt, gapID, hdays),
  pred   = cumPred(x,    wt, gapID, hdays),
  EDist  = getDist(x, y, wt, gapID, hdays)
)

# Rough plots                                                              #####
oldPar <- par(no.readonly = TRUE)

ptCI <- \(x, low, high, l = length(x), col = "gray90", border = NA, ...) {
  Larg <- list(..., col = col, border = border)
  Larg[["x"]] <- c(x, x[l], rev(x[-l]))
  Larg[["y"]] <- c(high, low[l], rev(low[-l]))
  do.call(polygon, Larg)
}

Y <- cumsum(y)
X <- cumsum(x)
l <- length(X)
wk <- seq_len(l) - 1L
yMax <- sapply(scL, \(x) x$pred$HCI[l]) |> c(X[l], Y[l]) |> max()
yl <- c(0, yMax + 1)

par(mfrow = c(2, 3), mar = c(3, 2, 3, 1), oma = c(2, 2, 0, 0), lwd = 2)
for (j in 1L:6L) {
  pred <- scL[[j]]$pred
  labTxt <- paste(j, scL[[j]]$desc, sep = ": ")
  plot(wk, type = "n", xlab = "", ylab = "", las = 1, ylim = yl)
  ptCI(wk, pred$LCI, pred$HCI)
  lines(wk, pred$MED, col = "black")
  lines(wk, X, col = "blue")
  lines(wk, Y, col = "red")
  mtext(labTxt, 3, adj = 0, padj = -0.3, font = 0, cex = .9)
}
mtext("Weeks",    1, outer = TRUE, font = 2)
mtext("Subjects", 2, outer = TRUE, las = 0, font = 2)

do.call(par, oldPar)
