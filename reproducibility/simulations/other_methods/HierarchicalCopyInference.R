args = commandArgs(trailingOnly=TRUE)

#file_name <- "10nodes_10regions_100000reads_sim1_d_mat.txt"

# Functions for evaluating clustering

# test example data for all the functions, remove row.names or else
# will produce warnings when using it on bootstrapping (complaining about
# duplicated row.names )

# mtcars_scaled <- scale(mtcars)
# row.names(mtcars_scaled) <- NULL

# ------------------------------------------------------------------------------------
#### Choosing the right k for clustering
library(dplyr)
library(tidyr)

## method 1. WSS :compute the total within sum square error, this measures how close
#  are the points in a cluster to each other

# [Distance] : calculates the sum squared distance of a given cluster of points,
#              note that "sum squared distance" is used here for measuring variance
Distance <- function(cluster)
{
	# the center of the cluster, mean of all the points
	center <- colMeans(cluster)

	# calculate the summed squared error between every point and
	# the center of that cluster
	distance <- apply( cluster, 1, function(row)
	{
		sum( ( row - center )^2 )
	}) %>% sum()

	return(distance)
}

# calculate the within sum squared error manually for hierarchical clustering
# [WSS] : pass in the dataset, and the resulting groups(cluster)
WSS <- function( data, groups )
{
	k <- max(groups)

	# loop through each groups (clusters) and obtain its
	# within sum squared error
	total <- lapply( 1:k, function(k)
	{
		# extract the data point within the cluster
		cluster <- subset( data, groups == k )

		distance <- Distance(cluster)
		return(distance)
	}) %>% unlist()

	return( sum(total) )
}

# testing
# sum_squared_error <- WSS( data = mtcars_scaled, groups =  groups )

# this value will will decrease as the number of clusters increases,
# because each cluster will be smaller and tighter.
# And the rate of the decrease will slow down after the optimal cluster number


## method 2 : Calinski-Harabasz index, ratio of the between cluster variance
#			  to the total within cluster variance
# http://www.mathworks.com/help/stats/clustering.evaluation.calinskiharabaszevaluation-class.html

# TSS (total sum of square) : the squared distance of all the data points from
# the dataset's centroid

# BSS (between sum of square) = TSS - WSS, measures how far apart are the clusters
# from each other
# !! a good clustering has a small WSS and a high BSS

# CHIndex = B / W, the ratio should be maximized at the optimal k
# B = BSS(k) / (k-1) ; k = # of cluster
# W = WSS(k) / (n-k) ; n = # of data points

# [CHCriterion] : calculates both Calinski-Harabasz index and within sum squared error
# @kmax          = maximum cluster number, caculates the CH index from 2 cluster to kmax

CHCriterion <- function( data, kmax, ...  )
{
	# total sum squared error (independent with the number of cluster k)
	tss <- Distance( cluster = data )

	# initialize a numeric vector storing the score
	wss <- numeric(kmax)

	# k starts from 2, cluster 1 is meaningless
	d <- dist( data, method = "euclidean" )
	clustering <- hclust( d, ... )

		for( k in 2:kmax )
		{
			groups <- cutree( clustering, k )
			wss[k] <- WSS( data = data, groups =  groups )
		}


	# between sum of square
	bss <- tss - wss[-1]

	# cluster count start from 2!
	numerator <- bss / ( 1:(kmax-1) )
	denominator <- wss[-1] / ( nrow(data) - 2:kmax )

	criteria <- data.frame( k = 2:kmax,
	                        CHIndex = numerator / denominator,
							wss = wss[-1] )

	k_selected <- criteria$k[which.max(criteria$CHIndex)[1]]

	groups <- cutree( clustering, k_selected )



	return(list( k=k_selected, groups=groups, criteria=criteria))
}


# ----------------------------------------------------------------------------------------------
# [ClusterMethod] : supports heirarchical clustering

# @data          = data frame type data, matrix also works
# @k             = specify the number of clusters
# @noise.cut     = if specified, the points of the resulting cluster whose number is smaller
#                  than it will be considered as noise, and all of these noise cluster will be
#                  grouped together as one whole cluster
# @...           = pass in other parameters for hclust or kmeans++ (same as kmeans)

ClusterMethod <- function( data, k, noise.cut = 0, ... )
{
		cluster   <- hclust( dist(data), ... )
		partition <- cutree( cluster, k )

	# equivalent to k
	cluster_num <- max(partition)

	# calculate each cluster's size
	cluster_size <- numeric(cluster_num)
	for( i in 1:cluster_num )
		cluster_size[i] <- sum( partition == i )

	# if there're cluster size smaller than the specified noise.cut, do :
	not_noise_num <- sum( cluster_size > noise.cut )

	if( cluster_num > not_noise_num )
	{
		# extract the cluster whose size is larger than noise.cut
		cluster_new <- (1:cluster_num)[ cluster_size > noise.cut ]

		# all the data points whose original cluster is smaller than the noise.cut
		# will be assigned to the same new cluster
		cluster_num <- not_noise_num + 1

		# new clustering number, assign the noise cluster's number first
		# then adjust the original cluster's number
		new <- rep( cluster_num, nrow(data) )

		for( i in 1:not_noise_num )
			new[ ( partition == cluster_new[i] ) ] <- i

		partition <- new
	}

	# boolean vector indicating which data point belongs to which cluster
	cluster_list <- lapply( 1:cluster_num, function(x)
	{
		return( partition == x )
	})

	cluster_result <- list( result      = cluster,
	                        partition   = partition,
	                        clusternum  = cluster_num,
	                        clusterlist = cluster_list )
	return(cluster_result)
}

# cluster_result <- ClusterMethod( data = mtcars_scaled, k = 5, clustermethod = "hclust" )



print(args)
file_name <- args[1]
output_file_name <- args[2]

# read in the simulated data

sim_data <- read.table(file_name,header=FALSE, sep = ",")

sim_data <- t(t(sim_data)) # make it numeric!

# size of the data

n_cells <- nrow(sim_data)
n_bins <- ncol(sim_data)

norm_data <- sim_data/rowSums(sim_data)

hclust_result <- CHCriterion(data = norm_data, kmax = 10, method = "ward.D")

inferred_states <- matrix(NA, ncol=n_cells, nrow=n_bins) # to store the output
# here is it transposed temporarily

for(kk in 1:hclust_result$k){
  selected_cells <- which(hclust_result$groups==kk)
  cluster_counts <- colSums(norm_data[selected_cells, ])
  cluster_profile <- 2*cluster_counts/median(cluster_counts) # set median to copy number state 2
  inferred_states[, selected_cells] <- cluster_profile
}

# round the inferred states to get integers states

inferred_states <- round(inferred_states)

# write the output inferred states
# where we transpose back

write.table(t(inferred_states), output_file_name, col.names=FALSE, row.names=FALSE)

# check the distance

#sim_truth <- read.table(paste0(file_name, "_ground_truth.txt"),header=FALSE)

#sim_truth <- t(t(sim_truth)) # make it numeric

#HMM_result <- read.table(paste0(file_name, "_HMMcopy_inferred.txt"),header=FALSE)

#HMM_result <- t(t(HMM_result)) # make it numeric

#mean((sim_truth-HMM_result)^2)
