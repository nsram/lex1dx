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

.histogram_lines <- function(result, x, weights = NULL, title = "") {
  tab <- bin_table(result, x, weights = weights)
  st <- .r2_mse(result, x, weights = weights)
  k <- nrow(tab)
  lines <- character()
  if (!is.null(title) && nzchar(title)) lines <- c(lines, as.character(title))
  lines <- c(
    lines,
    sprintf("Squared correlation r\u00b2 of %d-bin model with data = %.8f", k, st$r2),
    sprintf("MSE of %d-bin model with data = %.8f", k, st$mse),
    ""
  )
  n_max <- max(tab$n, 0)
  label_width <- 26L
  max_blocks <- 60L
  for (i in seq_len(k)) {
    label <- paste0(.sig5(tab$lo[i]), " - ", .sig5(tab$hi[i]))
    left <- sprintf("%-*s", label_width, paste0("    ", label))
    if (n_max > 0 && tab$n[i] > 0)
      blocks <- max(1L, as.integer(round(max_blocks * tab$n[i] / n_max)))
    else
      blocks <- 0L
    bar <- if (blocks > 0) paste(rep("\u2588", blocks), collapse = "") else ""
    annot <- sprintf("[%s, %s, %s]", .fmt_n(tab$n[i]), .fmt_pct(tab$pct[i]), .sig5(tab$mean[i]))
    sep <- if (blocks > 0) " " else ""
    lines <- c(lines, paste0(left, bar, sep, annot), "")
  }
  if (length(lines) && lines[length(lines)] == "")
    lines <- lines[-length(lines)]
  lines
}

histogram <- function(result, x, file = "lex1dx_histogram.png",
                      weights = NULL, title = "", theme = "default") {
  if (!theme %in% c("default", "light"))
    stop("theme must be 'default' or 'light'")
  if (theme == "default") {
    bg <- "#005fff"
    fg <- "#ffff00"
  } else {
    bg <- "#ffffff"
    fg <- "#000000"
  }
  lines <- .histogram_lines(result, x, weights = weights, title = title)
  nlines <- max(length(lines), 1L)
  width_ch <- max(c(40L, nchar(lines, type = "width"))) + 4L
  fig_w <- max(10, min(14, 0.095 * width_ch))
  fig_h <- 0.40 + 0.22 * nlines
  if (!is.null(file)) {
    opened <- FALSE
    tryCatch({
      grDevices::png(file, width = fig_w, height = fig_h, units = "in",
                     res = 120, type = "cairo", bg = bg)
      opened <- TRUE
    }, error = function(e) {
      grDevices::png(file, width = fig_w, height = fig_h, units = "in",
                     res = 120, bg = bg)
      opened <<- TRUE
    })
    on.exit(if (opened) grDevices::dev.off(), add = TRUE)
  }
  old <- graphics::par(mar = c(0.3, 0.3, 0.3, 0.3), oma = c(0, 0, 0, 0),
                       xaxs = "i", yaxs = "i", bg = bg, family = "mono")
  on.exit(graphics::par(old), add = TRUE)
  graphics::plot.new()
  graphics::plot.window(xlim = c(0, width_ch), ylim = c(nlines + 0.3, -0.3))
  usr_bg <- graphics::par("usr")
  graphics::rect(usr_bg[1], usr_bg[3], usr_bg[2], usr_bg[4], col = bg, border = NA)
  for (i in seq_along(lines)) {
    graphics::text(0, i - 1L, lines[i], adj = c(0, 1), col = fg, cex = 0.85)
  }
  if (is.null(file)) invisible(NULL) else invisible(file)
}
