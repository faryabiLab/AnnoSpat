#! /usr/bin/Rscript
my_packages <- c('spatstat','caret','optparse')                                        # Specify your packages
not_installed <- my_packages[!(my_packages %in% installed.packages()[ , "Package"])]    # Extract not installed packages
if(length(not_installed)) install.packages(not_installed)        

#install.packages(c('spatstat','caret','optparse'))

suppressMessages(library(spatstat))
suppressMessages(library(caret))
suppressMessages(library(reshape2))
suppressMessages(library(plyr))
suppressMessages(library(optparse))

# Get the valid items for a column (col) quantile (i) of n quantiles from a data
# frame (df).
filterQuantile = function(df, col, n, i) {
  quant = quantile(df[,col], probs = seq(0, 1, 1 / n))[i]
  valid = df[df[,col] >= quant,"item"]

  if(quant == 0) {
    return(c())
  } else {
    return(valid)
  }
}

# Get the subset of valid items for all columns (col) quantile (i) of n
# quantiles from a data frame (df).
filterDfOnQuantile = function(df, cols, n, i) {
  valid = sapply(cols, function(col) filterQuantile(df, col, n, i))
  validDf = droplevels(df[df$item %in% valid,])

  return(validDf)
}

# Attach and filter labels (item,label) based on a label reference. Everything
# not marked valid is labeled "Other".
labelAndFilterDf = function(df, labelDf, valids) {
  mergedDf = merge(df, labelDf, by=c("item"))
  # Remove if not keeping ALL.
  if(all(toupper(valids) == "ALL")) {
    validDf = mergedDf
  } else {
    levels(mergedDf$label) = c(levels(mergedDf$label), "Other")
    mergedDf[!(mergedDf$label %in% valids),c("label")] = as.factor("Other")
    validDf = mergedDf
  }

  return(droplevels(validDf))
}

# Find index of first switch from above to below or equal value.
findNegSwap = function(x, xs) {
  if(all(xs$obs > x)) {
    return(xs$r[length(xs$r)])
  } else if(all(xs$obs <= x)) {
    return(xs$r[1])
  } else {
    for (i in c(1:length(xs$r))) {
      if (i == length(xs$r)) {
        return(xs$r[i])
      } else if ((xs$obs[i] > x) && (xs$obs[i + 1] <= x)) {
        return(xs$r[i + 1])
      }
    }
  }
}

# Find index of first switch from below or equal value to above.
findPosSwap <- function(x, xs) {
  if (all(xs$obs <= x)) {
    return(xs$r[length(xs$r)])
  } else if (all(xs$obs > x)) {
    return(xs$r[1])
  } else {
    for (i in c(1:length(xs$r))) {
      if (i == length(xs$r)) {
        return(xs$r[i])
      } else if (xs$obs[i] <= x && xs$obs[i + 1] > x) {
        return(xs$r[i + 1])
      }
    }
  }
}

# Find the longest length of a continuously positive (or negative / equal to) stretch.
findLongestLength <- function(posFlag, x, xs) {
  if(posFlag) {
    lengthInfo = rle(xs$obs > x)
  } else {
    lengthInfo = rle(xs$obs <= x)
  }
  stretches = lengthInfo$lengths[which(lengthInfo$values)]
  if(length(stretches) == 0) {
    return(0)
  } else {
    return(xs$r[max(stretches)] - xs$r[1])
  }
}

# Find the first r value of the maximum value.
findMaxPos <- function(xs) {
  return(xs$r[min(which(xs$obs == max(xs$obs)))])
}

# Find the first r value of the minimum value.
findMinPos <- function(xs) {
  return(xs$r[min(which(xs$obs == min(xs$obs)))])
}

main = function(marks, labelFile, outFolder, inFile) {
  inDf = read.csv(inFile)

  dir.create(outFolder, showWarnings = FALSE, recursive = TRUE)

  if (file.exists(labelFile)) {
    print("Found label file")
    labelDf = read.csv(labelFile)
    df = labelAndFilterDf(inDf, labelDf, marks)
  } else {
    df = filterDfOnQuantile(inDf, marks, 10, 10)
  }

  maxWindowSize = max(max(df$CenterX) - min(df$CenterX), max(df$CenterY) - min(df$CenterY))
  # pat = ppp(df$CenterX, df$CenterY, c(-500,500), c(-500,500))
  pat = ppp(df$CenterX, df$CenterY, c(min(df$CenterX), max(df$CenterX)), c(min(df$CenterY),max(df$CenterY)))

  if (file.exists(labelFile)) {
    marks(pat) = df$label
  } else {
    # Add min-max normalized data as marks.
    marksDf = df[, marks]
    marks(pat) = predict(preProcess(marksDf, method = c("range")), newdata = marksDf)
  }

  pdf(file = file.path(outFolder, "basic_plot.pdf"))
  plot(pat)
  dev.off()

  # Marked basic plot.
  pdf(file = file.path(outFolder, "marked_plot.pdf"))
  plot(pat)
  dev.off()

  # Mark correlation function.
  pdf(file = file.path(outFolder, "mark_correlation_function.pdf"))
  plot(markcorr(pat))
  dev.off()

  # Mark cross correlation function (several columns of marks).
  pdf(file = file.path(outFolder, "mark_cross_correlation_function.pdf"))
  ## crossFn = markcrosscorr(pat)
  ## crossFn = envelopeArray(pat, markcrosscorr, nsim=9, savefuns=TRUE)
  if(all(toupper(marks) == "ALL")) {
    ## crossFn = markcrosscorr(pat, method = "loess")
    ## for(i in c(1:length(crossFn$fns))) {
    ##   crossFn$fns[[i]]$obs = crossFn$fns[[i]]$iso
    ## }
    crossFn = envelopeArray(pat, markcrosscorr, r = seq(0, 250, 250 / 5000), nsim = 9, savefuns = TRUE, method = "density")
  } else {
    ## crossFn = envelopeArray(pat, markcrosscorr, r=seq(0,250), nsim=9, savefuns=TRUE, method = "density")
    crossFn = markcrosscorr(pat)
    for(i in c(1:length(crossFn$fns))) {
      crossFn$fns[[i]]$obs = crossFn$fns[[i]]$iso
    }
  }
  plot(crossFn)
  dev.off()

  # Mark variogram (lower value, more similar mark values at a distance).
  if (!file.exists(labelFile)) {
    pdf(file = file.path(outFolder, "mark_variogram.pdf"))
    plot(markvario(pat))
    dev.off()
  }

  # Envelope.
  pdf(file = file.path(outFolder, "envelope.pdf"))
  plot(envelope(pat))
  dev.off()

  # Object
  outObject = crossFn
  outObject$label = basename(outFolder)
  saveRDS(outObject, file = file.path(outFolder,"crosscorr.rds"))

  # Statistics
  statsDf = melt(crossFn$which)
  statsDf$meanCorr = sapply(statsDf$value, function(x) mean(crossFn$fns[[x]]$obs))
  statsDf$maxCorr = sapply(statsDf$value, function(x) max(crossFn$fns[[x]]$obs))
  statsDf$minCorr = sapply(statsDf$value, function(x) min(crossFn$fns[[x]]$obs))
  statsDf$topMaxCorr = sapply(statsDf$value, function(x) max(crossFn$fns[[x]]$obs[1:round(length(crossFn$fns[[x]]$obs) / 4)]))
  statsDf$topMeanCorr = sapply(statsDf$value, function(x) mean(crossFn$fns[[x]]$obs[1:round(length(crossFn$fns[[x]]$obs) / 4)]))
  statsDf$negSwap = sapply(statsDf$value, function(x) findNegSwap(1, crossFn$fns[[x]]))
  statsDf$posSwap = sapply(statsDf$value, function(x) findPosSwap(1, crossFn$fns[[x]]))
  statsDf$longestPosLength = sapply(statsDf$value, function(x) findLongestLength(TRUE, 1, crossFn$fns[[x]]))
  statsDf$longestNegLength = sapply(statsDf$value, function(x) findLongestLength(FALSE, 1, crossFn$fns[[x]]))
  statsDf$maxPos = sapply(statsDf$value, function(x) findMaxPos(crossFn$fns[[x]]))
  statsDf$minPos = sapply(statsDf$value, function(x) findMinPos(crossFn$fns[[x]]))
  statsDf$label = basename(outFolder)

  # To get the sample size, different for quantitative vs. nominal marks
  markVars = marks(pat)
  if (is.data.frame(markVars)) {
    statsDf$n = nrow(markVars)
  } else {
    statsDf$n = length(markVars[markVars != "Other"])
  }

  # Output statistics
  write.csv(statsDf, file.path(outFolder, "stats.csv"), quote = FALSE, row.names = FALSE)

  # Curve
  curveDf <- adply(melt(crossFn$which), 1, function(xs) {
    data.frame(x = crossFn$fns[[xs$value]]$r, y = crossFn$fns[[xs$value]]$obs, label = basename(outFolder))
  })
  write.csv(curveDf, file.path(outFolder, "curve.csv"), quote = FALSE, row.names = FALSE)

  }

option_list = list(
  make_option(c("-i", "--input"), type="character",
              help="Input file for spatial data", metavar="character"),
  make_option(c("-o", "--out"), type="character", default="out.txt",
              help="[%default] Output path", metavar="character"),
  make_option(c("-l", "--labels"), type="character",
              help="Optional input file for label data", metavar="character"),
  make_option (c("-m","--marks"), help="Comma separated list of marks")
)

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
print(opt)

marks = unlist(strsplit(opt$marks, ","))

main(marks, opt$labels, opt$out, opt$input)
