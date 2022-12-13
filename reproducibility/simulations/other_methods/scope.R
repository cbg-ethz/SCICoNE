suppressWarnings(library(SCOPE))
suppressWarnings(library(GenomicRanges))
suppressWarnings(library(DNAcopy))

initialize_ploidy <- function(Y, Yhat, ref, maxPloidy = 6,
                              minPloidy = 1.5, minBinWidth = 5, SoS.plot = FALSE) {
  ploidy.SoS <- rep(NA, ncol(Y))
  
  breaks <- matrix(0, nrow(Y), ncol(Y))
  RCNP <- matrix(0, nrow(Y), ncol(Y))
  final <- matrix(0, nrow(Y), ncol(Y))
  X <- seq(minPloidy, maxPloidy, by = 0.05)
  n_ploidy <- length(X)
  SoS <- matrix(0, n_ploidy, ncol(Y))
  
  normal <- (Y + 1)/(Yhat + 1)
  
  for (k in seq_len(ncol(Y))) {
    if (k%%5 == 1) {
      cat(k, "\t")
    }
    
    lr <- log(normal[, k])
    loc <- data.frame(seq = as.character(seqnames(ref)),
                      start = start(ref), end = end(ref))
    CNA.object <- CNA(genomdat = lr, chrom = loc[, 1],
                      maploc = as.numeric(loc[, 2]), data.type = "logratio")
    CNA.smoothed <- smooth.CNA(CNA.object)
    segs <- segment(CNA.smoothed, verbose = 0,
                    min.width = minBinWidth, alpha=0.01)
    frag <- segs$output[, 2:3]
    len <- dim(frag)[1]
    bps <- array(0, len)
    for (j in seq_len(len)) {
      bps[j] <- which((loc[, 1] == frag[j, 1]) &
                        (as.numeric(loc[, 2]) == frag[j, 2]))
    }
    bps <- sort(bps)
    bps[(len = len + 1)] <- nrow(Y)
    breaks[bps, k] <- 1
    RCNP[, k][seq_len(bps[2])] <- median(normal[,
                                                k][seq_len(bps[2])])
    if (len > 2) {
      for (i in 2:(len - 1)) {
        RCNP[, k][bps[i]:(bps[i + 1] - 1)] <- median(normal[,
                                                            k][bps[i]:(bps[i + 1] - 1)])
      }
    }
    RCNP[, k] <- RCNP[, k]/mean(RCNP[, k])
    
    SCNP <- RCNP[, k] %o% X
    FSCP <- round(SCNP)
    Diff2 <- (SCNP - FSCP)^2
    SoS[, k] <- colSums(Diff2, na.rm = FALSE, dims = 1)
    ploidy.SoS[k] <- X[which.min(SoS[, k])]
    
    
    if(SoS.plot){
      par(mfrow = c(1,2))
      par(mar = c(5,4,4,2))
      hist(Y[,k], 100, main = 'Read depth distribution', 
           xlab = 'Coverage per bin')
      plot(X, SoS[,k], xlab = "ploidy", ylab = "Sum of squared errors", 
           main = "First-pass estimation of ploidy", pch = 16)
      abline(v = X[which.min(SoS[,k])], lty = 2)
    }
  }
  return(ploidy.SoS)
}

segment_data <- function(x, size) {
  n_bins <- nrow(x)
  rowsum(x, rep(1:(n_bins/size),each=size))
}

explode_cnvs <- function(x, size) {
  x[,do.call(c, lapply(1:ncol(x), rep, size))]
}


args = commandArgs(trailingOnly=TRUE)
file_name <- args[1]

# Read in data
sim_data <- as.matrix(read.table(file_name,header=FALSE, sep=','))
sim_data <- t(sim_data)
n_bins <- nrow(sim_data)
n_cells <- ncol(sim_data)
colnames(sim_data) <- paste0("cell", c(1:n_cells))
seg_size <- 10
sim_data <- segment_data(sim_data, seg_size)

# Create reference
starts <- seq(1,n_bins,seg_size)
ends <- seq(seg_size,n_bins,seg_size)
df <- data.frame(chr="chr1", start=starts, end=ends, gc=runif(length(starts)) + 40)
sim_ref <- makeGRangesFromDataFrame(df)  # strand value "." is replaced with "*"
values(sim_ref) <- cbind(values(sim_ref), gc=df$gc)

# Normalize
Gini <- get_gini(sim_data)
diploid_cells_idx <- which(Gini < quantile(Gini, probs=c(.1)))

# first-pass CODEX2 run with no latent factors
normObj.sim <- normalize_codex2_ns_noK(Y_qc = sim_data,
                                       gc_qc = sim_ref$gc,
                                       norm_index = diploid_cells_idx)

Yhat.noK.sim <- normObj.sim$Yhat
beta.hat.noK.sim <- normObj.sim$beta.hat
fGC.hat.noK.sim <- normObj.sim$fGC.hat
N.sim <- normObj.sim$N

# Ploidy initialization
#ploidy.sim <- initialize_ploidy(Y = sim_data, Yhat = Yhat.noK.sim, ref = sim_ref, minBinWidth = 5)
# Assume ploidy is 2 for all cells because estimates were way off
ploidy.sim <- rep(2.0, n_cells)

normObj.scope.sim <- normalize_scope_foreach(Y_qc = sim_data, gc_qc = sim_ref$gc,
                                             K = 1, ploidyInt = ploidy.sim,
                                             norm_index = diploid_cells_idx, T = 1:20,
                                             beta0 = beta.hat.noK.sim, nCores = 10)
Yhat.sim <- normObj.scope.sim$Yhat[[which.max(normObj.scope.sim$BIC)]]
fGC.hat.sim <- normObj.scope.sim$fGC.hat[[which.max(normObj.scope.sim$BIC)]]

# Infer CNVs
inf_cnvs <- segment_CBScs(Y = sim_data,
			    Yhat = Yhat.sim,
			    sampname = colnames(sim_data),
			    ref = sim_ref,
			    chr = "chr1",
			    mode = "integer", max.ns = 1)


init_cnv_mat <- t(inf_cnvs$iCN)
cnv_mat <- explode_cnvs(init_cnv_mat, seg_size)

# Write to output
write.table(cnv_mat, args[2], col.names=FALSE, row.names=FALSE)
