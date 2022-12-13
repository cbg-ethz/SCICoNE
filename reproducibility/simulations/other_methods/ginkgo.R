suppressWarnings(library(DNAcopy))

args = commandArgs(trailingOnly=TRUE)
#file_name <- "10nodes_10regions_100000reads_sim1"
print(args)
file_name <- args[1]
output_file_name <- args[2]
raw <- read.table(file_name, header=FALSE, sep=',')

raw <- t(raw)

# raw <- matrix(rpois(200, lambda=10), nrow=20, ncol=10)
l <- dim(raw)[1] # number of bins
w <- dim(raw)[2] # number of cells

minPloidy   = 0.
maxPloidy   = 4
minBinWidth = 5

breaks       = matrix(0,l,w)
fixed        = matrix(0,l,w)
final        = matrix(0,l,w)
stats        = matrix(0,w,10)

CNgrid       = seq(minPloidy, maxPloidy, by=0.05)
n_ploidy     = length(CNgrid)  # Number of ploidy tests during CN inference
CNmult       = matrix(0,n_ploidy,w)
CNerror      = matrix(0,n_ploidy,w)
outerColsums = matrix(0,n_ploidy,w)

normal  = sweep(raw+1, 2, colMeans(raw+1), '/')
normal2 = normal
lab     = colnames(normal)

rownames(stats) = lab
colnames(stats) = c("Reads", "Bins", "Mean", "Var", "Disp", "Min", "25th", "Median", "75th", "Max")

chrom = rep(1,l)
# loc = seq(5,l*minBinWidth,minBinWidth)
loc = seq(1, l)
F = normal[,which.min(apply(normal, 2, sd)/apply(normal,2,mean))[1]]
for (k in 1:w) { # For each cell
  # Generate basic statistics
  stats[k,1]  = sum(raw[,k])
  stats[k,2]  = l
  stats[k,3]  = round(mean(raw[,k]), digits=2)
  stats[k,4]  = round(sd(raw[,k]), digits=2)
  stats[k,5]  = round(stats[k,4]/stats[k,3], digits=2)
  stats[k,6]  = min(raw[,k])
  stats[k,7]  = quantile(raw[,k], c(.25))[[1]]
  stats[k,8]  = median(raw[,k])
  stats[k,9]  = quantile(raw[,k], c(.75))[[1]]
  stats[k,10] = max(raw[,k])

  lr = log2((normal[,k])/(F))

  CNA.object   = CNA(genomdat = lr, chrom=chrom, maploc = as.numeric(loc), data.type = 'logratio')
  CNA.smoothed = smooth.CNA(CNA.object)
  segs         = segment(CNA.smoothed, verbose=0)
  frag         = segs$output[,2:3]

  # Map breakpoints to kth sample
  len = dim(frag)[1]
  bps = array(0, len)
  for (j in 1:len)
    bps[j]=which((chrom==frag[j,1]) & (as.numeric(loc)==frag[j,2]))
  bps = sort(bps)
  bps[(len=len+1)] = l

  # Track global breakpoint locations
  breaks[bps,k] = 1

  # Modify bins to contain median read count/bin within each segment
  fixed[,k][1:bps[2]] = median(normal[,k][1:bps[2]])
  if (len > 2) {
    for(i in 2:(len-1))
      fixed[,k][bps[i]:(bps[i+1]-1)] = median(normal[,k][bps[i]:(bps[i+1]-1)])
    fixed[,k] = fixed[,k]/mean(fixed[,k])
  }
  # Determine Copy Number
  outerRaw         = fixed[,k] %o% CNgrid
  outerRound       = round(outerRaw)
  outerDiff        = (outerRaw - outerRound) ^ 2
  outerColsums[,k] = colSums(outerDiff, na.rm = FALSE, dims = 1)
  CNmult[,k]       = CNgrid[order(outerColsums[,k])]
  CNerror[,k]      = round(sort(outerColsums[,k]), digits=2)

  CN = CNmult[which(abs(CNmult[,k] - 2)<.4),k][1]

  final[,k] = round(fixed[,k]*CN)
  final[l,k] = final[l-1,k] # copy last state
}

write.table(t(final), output_file_name, col.names=FALSE, row.names=FALSE)
