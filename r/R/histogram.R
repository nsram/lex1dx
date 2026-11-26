.sig5 <- function(v) {
  formatC(signif(as.numeric(v), 5), format = "g", digits = 5)
}

.fmt_pct <- function(pct) {
  s <- sub("\\.?0+$", "", sprintf("%.5f", pct))
  paste0(s, "%")
}

.fmt_n <- function(n) {
  if (isTRUE(all.equal(n, round(n), tolerance = 0))) as.character(as.integer(round(n)))
  else formatC(n, format = "g", digits = 8)
}

bin_table <- function(result, x, weights = NULL) {
  if (!is.numeric(x) || (!is.null(dim(x)) && length(dim(x)) > 1L))
    stop("x must be a numeric vector")
  if (length(x) == 0L) stop("x is empty")
  if (is.null(result$bounds) || is.null(result$k) && is.null(result$sizes))
    stop("result must be a lex1dx cluster() result")
  k <- length(result$sizes)
  b <- as.integer(result$bounds)
  if (b[1L] != 0L || b[length(b)] != length(x))
    stop("result$bounds must cover x")
  o <- order(x, method = "radix")
  sx <- as.numeric(x)[o]
  if (is.null(weights)) {
    sw <- NULL
    total <- length(x)
  } else {
    if (!is.numeric(weights) || length(weights) != length(x))
      stop("weights size must match the input length")
    sw <- as.numeric(weights)[o]
    total <- sum(sw)
  }
  lo <- hi <- n <- pct <- mean <- numeric(k)
  for (i in seq_len(k)) {
    a <- b[i] + 1L
    z <- b[i + 1L]
    lo[i] <- sx[a]
    hi[i] <- sx[z]
    if (is.null(sw)) {
      n[i] <- z - a + 1L
      mean[i] <- mean(sx[a:z])
    } else {
      n[i] <- sum(sw[a:z])
      mean[i] <- if (n[i] > 0) sum(sw[a:z] * sx[a:z]) / n[i] else 0
    }
    pct[i] <- if (total > 0) 100 * n[i] / total else 0
  }
  data.frame(lo = lo, hi = hi, n = n, pct = pct, mean = mean)
}

.r2_mse <- function(result, x, weights = NULL) {
  xv <- as.numeric(x)
  if (is.null(weights)) {
    mu <- mean(xv)
    sst <- sum((xv - mu)^2)
    denom <- length(xv)
  } else {
    wv <- as.numeric(weights)
    wsum <- sum(wv)
    mu <- if (wsum > 0) sum(wv * xv) / wsum else 0
    sst <- sum(wv * (xv - mu)^2)
    denom <- wsum
  }
  sse <- as.numeric(result$sse)
  r2 <- if (sst == 0) 1 else 1 - sse / sst
  mse <- if (denom == 0) 0 else sse / denom
  list(r2 = r2, mse = mse)
}

.histogram_rows <- function(result, x, weights = NULL, title = "") {
  tab <- bin_table(result, x, weights = weights)
  st <- .r2_mse(result, x, weights = weights)
  k <- nrow(tab)
  header <- character()
  if (!is.null(title) && nzchar(title)) header <- c(header, as.character(title))
  header <- c(
    header,
    sprintf("Squared correlation r\u00b2 of %d-bin model with data = %.8f", k, st$r2),
    sprintf("MSE of %d-bin model with data = %.8f", k, st$mse)
  )
  n_max <- max(tab$n, 0)
  label <- sprintf("%-26s", paste0("    ", .sig5(tab$lo), " - ", .sig5(tab$hi)))
  frac <- if (n_max > 0) tab$n / n_max else rep(0, k)
  annot <- sprintf("[%s, %s, %s]", vapply(tab$n, .fmt_n, ""),
                   vapply(tab$pct, .fmt_pct, ""), .sig5(tab$mean))
  list(header = header, label = label, frac = frac, annot = annot)
}

# Bars are drawn with rect(), not text glyphs, so the plot needs only ASCII
# text and works on any png() device (quartz, cairo, Xlib).
histogram <- function(result, x, file = NULL,
                      weights = NULL, title = "", theme = "default") {
  if (!theme %in% c("default", "light"))
    stop("theme must be 'default' or 'light'")
  if (theme == "default") {
    bg <- "#005fff"
    fg <- "#ffff00"
    bar <- "#ffff00"
  } else {
    bg <- "#ffffff"
    fg <- "#000000"
    bar <- "#555555"
  }
  rows <- .histogram_rows(result, x, weights = weights, title = title)
  k <- length(rows$label)
  nh <- length(rows$header)
  max_bar <- 60
  width_ch <- 26 + max_bar + 2 + max(nchar(rows$annot), 10) + 2
  nlines <- nh + 1 + 2 * k - 1
  fig_w <- max(10, min(14, 0.095 * width_ch))
  fig_h <- 0.40 + 0.22 * nlines
  opened <- !is.null(file)
  if (opened) {
    grDevices::png(file, width = fig_w, height = fig_h, units = "in",
                   res = 120, bg = bg)
  }
  old <- graphics::par(mar = c(0.3, 0.3, 0.3, 0.3), oma = c(0, 0, 0, 0),
                       xaxs = "i", yaxs = "i", bg = bg, family = "mono")
  # Restore par before dev.off(); the reverse order would open Rplots.pdf
  # in a session with no other device (e.g. Rscript).
  on.exit({
    graphics::par(old)
    if (opened) grDevices::dev.off()
  }, add = TRUE)
  graphics::plot.new()
  graphics::plot.window(xlim = c(0, 1), ylim = c(nlines + 0.3, -0.3))
  usr <- graphics::par("usr")
  graphics::rect(usr[1], usr[3], usr[2], usr[4], col = bg, border = NA)
  cex <- 0.85
  ch <- graphics::strwidth("0", cex = cex)
  # Shrink the bar scale if the widest row would not fit.
  need <- (26 + 2 + max(nchar(rows$annot))) * ch
  scale <- min(max_bar * ch, max(0.1, 0.98 - need))
  for (i in seq_len(nh))
    graphics::text(0, i - 1, rows$header[i], adj = c(0, 1), col = fg, cex = cex)
  for (i in seq_len(k)) {
    y <- nh + 1 + 2 * (i - 1)
    graphics::text(0, y, rows$label[i], adj = c(0, 1), col = fg, cex = cex)
    x0 <- 26 * ch
    w <- if (rows$frac[i] > 0) max(ch / 2, scale * rows$frac[i]) else 0
    if (w > 0)
      graphics::rect(x0, y + 0.05, x0 + w, y + 0.7, col = bar, border = NA)
    graphics::text(x0 + w + ch, y, rows$annot[i], adj = c(0, 1), col = fg, cex = cex)
  }
  if (is.null(file)) invisible(NULL) else invisible(file)
}
