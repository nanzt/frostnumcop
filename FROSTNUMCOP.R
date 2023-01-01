##############################################################
# sample code of Cao et al., 2022. ESSD
# estimate permafrost distribution based on survey-based maps
# inputs: thawing/freezing indices, 
#         survey-based permafrost maps, 
#         and soil cluster distributions
# output: permfrost distribution maps
# Zetao Cao, 2022-10-26
##############################################################

library(pso)
library(sp)
library(raster)
source("functions.R")

permaf <- 1 # value to represent permafrost
season <- 2 # value to represent seasonally frozen ground
n_size <- 5 # window size for neighboring cells 

# random weights
set.seed(123)
weight.kappa <- runif(1,5,8)
weight.boundary <- 10 - weight.kappa

# max iterations of PSO 
n_iters <- 100

# set data path and read raster
work_dir <- "sample_data/"

soil_cluster_file <- paste0(work_dir, "soil_cluster.tif")
inves_map_file <- paste0(work_dir,"investigation_map.tif")
DDT_file <- paste0(work_dir,"DDT_mean_annual.tif")
DDF_file <- paste0(work_dir,"DDF_mean_annual.tif")

out_map_file <- paste0(work_dir, "map.tif")

DDT_ras <- raster(DDT_file)
DDF_ras <- raster(DDF_file)

# avoid zero divisors
DDT_ras[DDT_ras==0] <- 1
DDF_ras[DDF_ras==0] <- 1

soil_cluster <- raster(soil_cluster_file)
inves_map <- raster(inves_map_file)
sqrt_DD <- sqrt(DDT_ras/DDF_ras)

save_boundary <- F # save boundary cells
boundary_file <- paste0(work_dir, "boundary.tif")
boundary_ras <- inves_map

out_ras <- soil_cluster

# the raster extents should be the same
if(length(soil_cluster) != length(inves_map) | length(soil_cluster) != length(sqrt_DD)){
  print("error! rasters' extent different!!!")
}

# get col and row
ncols <- soil_cluster@ncols
nrows <- soil_cluster@nrows

# form data frame
index <- 1:length(soil_cluster)
soil_cluster <- as.vector(soil_cluster)
inves_map <- as.vector(inves_map)
sqrt_DD <- as.vector(sqrt_DD)
df.all <- data.frame(index, soil_cluster, inves_map, sqrt_DD)
df.inves <- na.omit(df.all)

# get boundary raster 
#(grid with a neighborhood which contain different types of soil clusters and frozen soil)
df.inves$is.boundary <- 0
boundary_info_names <- c("index","soil_cluster", "sqrt_DD")
boundary_neighbors_permaf <- as.data.frame(matrix(nrow = 0, ncol = 3))
boundary_neighbors_season <- as.data.frame(matrix(nrow = 0, ncol = 3))
colnames(boundary_neighbors_permaf) <- boundary_info_names
colnames(boundary_neighbors_season) <- boundary_info_names

print("get boundary cells.")
for ( i in 1:length(df.inves$index)){
  grid_index <- df.inves$index[i]
  # get neighborhood
  neighbor_index <- get_neighbor(grid_index, n_size = n_size, ras_ncols = ncols, ras_nrows = nrows)
  df.neighbor <- df.all[neighbor_index,]
  df.neighbor <- na.omit(df.neighbor)
  # determine whether boundary grid satisfy the requirement
  p_l <- length(df.neighbor$inves_map[df.neighbor$inves_map == permaf])
  if ( (p_l != length(df.neighbor$inves_map) & p_l > 0)){
    df.permaf <- subset(df.neighbor, df.neighbor$inves_map == permaf )
    df.season <- subset(df.neighbor, df.neighbor$inves_map == season )
    # if met, copy neighboring cell's information into data.frame
    if ( (length(unique(df.permaf$soil_cluster)) > 1) | (length(unique(df.season$soil_cluster)) > 1) ){
      #print(grid_index)
      df.inves$is.boundary[i] <- 1
      df.permaf$index <- grid_index
      df.season$index <- grid_index
      df.permaf <- df.permaf[boundary_info_names]
      df.season <- df.season[boundary_info_names]
      boundary_neighbors_permaf <- rbind(boundary_neighbors_permaf, df.permaf)
      boundary_neighbors_season <- rbind(boundary_neighbors_season, df.season)
    }
  }
}
# amount of boundary grids
boundary_length <- length(unique(boundary_neighbors_permaf$index))

# save boundary raster 
if (save_boundary == T){
  boundary_ras[df.inves$index] <- df.inves$is.boundary
  writeRaster(boundary_ras, boundary_file, overwrite=TRUE)
}


soil_cluster_index <- sort(unique(df.inves$soil_cluster))
# number of soil c
soil_num <- length(soil_cluster_index)
e_pars <- rep(1,soil_num)
cluster_names <- c("inves_map","esti_type")


obj_fun <- function(e_pars){
  ########################################################
  # calculate kappa coefficient under a set of E parameter
  ########################################################
  for ( i in 1:soil_num){
    df.inves$E[df.inves$soil_cluster == soil_cluster_index[i]] <- e_pars[i]
  }
  # estimate frozen soil type by Frost number and E parameter
  df.inves$F2 <- df.inves$sqrt_DD*df.inves$E
  df.inves$esti_type[df.inves$F2 >= 1] <- season #2
  df.inves$esti_type[df.inves$F2 <  1] <- permaf #1
  kappa_v <- cal_kappa_2v_(df.inves, cluster_names)
  
  ########################################################
  # calculate boundary accuracy under a set of E parameter
  ########################################################
  for (i in 1:soil_num ){
    boundary_neighbors_permaf$E[boundary_neighbors_permaf$soil_cluster == soil_cluster_index[i]] <- e_pars[i]
    boundary_neighbors_season$E[boundary_neighbors_season$soil_cluster == soil_cluster_index[i]] <- e_pars[i]
  }
  
  # calculate F2
  boundary_neighbors_permaf$F2 <- boundary_neighbors_permaf$E*boundary_neighbors_permaf$sqrt_DD
  boundary_neighbors_season$F2 <- boundary_neighbors_season$E*boundary_neighbors_season$sqrt_DD
  
  f_F2 <- aggregate(boundary_neighbors_permaf$F2, by = list(boundary_neighbors_permaf$index), FUN = mean)
  s_F2 <- aggregate(boundary_neighbors_season$F2, by = list(boundary_neighbors_season$index), FUN = mean)
  
  colnames(f_F2) <- c("index", "f_F2")
  colnames(s_F2) <- c("index", "s_F2")
  df.F2 <- merge(f_F2, s_F2, by = c("index"))

  correct_boundary_length <- length(which(df.F2$f_F2 < df.F2$s_F2))
  
  beta <- correct_boundary_length/boundary_length
  
  ob.value <- weight.kappa*(1-kappa_v) + weight.boundary*(1-beta) #PSO searches the minimum value in default.
  
  return (ob.value)
}

print("PSO...")
control = list(maxit = n_iters)
p <- psoptim(e_pars, obj_fun, lower = 0.5, upper = 1.5, control = control)
### write estimate permafrost map
e_pars_new <- p$par

# store E values
for (i in 1:soil_num){
  df.all$E[df.all$soil_cluster == soil_cluster_index[i]] <- e_pars_new[i]
}

df.all$F2 <- df.all$sqrt_DD*df.all$E
df.all$esti_type[df.all$F2 >= 1] <- season #2
df.all$esti_type[df.all$F2 <  1] <- permaf #1

# permafrost map
out_ras[df.all$index] <- df.all$esti_type

writeRaster(out_ras, out_map_file, overwrite=TRUE)

