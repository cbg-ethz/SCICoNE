
args = commandArgs(trailingOnly=TRUE)

#file_name <- "10nodes_10regions_100000reads_sim1"

print(args)
file_name <- args[1]


# HMMCopy script
suppressWarnings(library(HMMcopy))

# build fake data frame for HMMcopy

rfile <- system.file("extdata", "normal.wig", package = "HMMcopy")

mfile <- system.file("extdata", "map.wig", package = "HMMcopy")
gfile <- system.file("extdata", "gc.wig", package = "HMMcopy")

normal_reads <- wigsToRangedData(rfile, gfile, mfile)

tfile <- system.file("extdata", "tumour.wig", package = "HMMcopy")
tumour_copy <- correctReadcount(wigsToRangedData(tfile, gfile, mfile))

# read in the simulated data

sim_data <- read.table(file_name,header=FALSE)

sim_data <- t(t(sim_data)) # make it numeric!

# size of the data

n_cells <- nrow(sim_data)
n_bins <- ncol(sim_data)

inferred_states <- matrix(NA, nrow=n_cells, ncol=n_bins) # to store the output

# use the HMM copy example to make the fake data frame

tumour_copy_small <- tumour_copy[1:n_bins,] 

tumour_copy_small$gc <- 1
tumour_copy_small$map <- 1
tumour_copy_small$valid <- TRUE
tumour_copy_small$ideal <- FALSE
tumour_copy_small$cor.gc <- 1
tumour_copy_small$cor.map <- 1

for(ii in 1:n_cells){ # run over all cells

  print(paste("Cell number",ii, "being processed"))
  
  tumour_copy_small$copy <- sim_data[ii,] # put the ii cell in there instead

  tumour_segments <- HMMsegment(tumour_copy_small, verbose=FALSE) # get the parameters

  inferred_states[ii,] <- tumour_segments$state
}

# write the output inferred states

file_name = substr(file_name, 1, nchar(file_name)-10)
print(file_name)

write.table(inferred_states, paste0(file_name, "_HMMcopy_inferred.txt"), col.names=FALSE, row.names=FALSE)

# check the distance

#sim_truth <- read.table(paste0(file_name, "_ground_truth.txt"),header=FALSE)

#sim_truth <- t(t(sim_truth)) # make it numeric

#HMM_result <- read.table(paste0(file_name, "_HMMcopy_inferred.txt"),header=FALSE)

#HMM_result <- t(t(HMM_result)) # make it numeric

#mean((sim_truth-HMM_result)^2)
