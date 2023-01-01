###########################################
# functions used permafrost mapping
# Zetao Cao, 2021-11-25
###########################################

# calculate kappa coefficient
cal_kappa_2v <- function(df, cluster_names){
  df <- df[cluster_names]
  conf_mat <- table(df)
  i1e1 <- conf_mat[1,1]
  i2e2 <- conf_mat[2,2]
  i1e2 <- conf_mat[1,2]
  i2e1 <- conf_mat[2,1]
  sum_n <- sum(conf_mat)
  po <- (i1e1+i2e2)/sum_n
  pe <- (sum(conf_mat[1,])/sum_n) * (sum(conf_mat[,1])/sum_n)  + (sum(conf_mat[2,])/sum_n) * (sum(conf_mat[,2])/sum_n)
  k <- (po-pe)/(1-pe)
  return (k)
}

# calculate kappa coefficient (divide first to avoid large number)
cal_kappa_2v_ <- function(df, cluster_names){
  df <- df[cluster_names]
  conf_mat <- table(df)
  sum_n <- sum(conf_mat)
  conf_mat <- conf_mat/sum_n
  po <- conf_mat[1,1] + conf_mat[2,2]
  pe <- sum(conf_mat[1,] * sum(conf_mat[,1])  + sum(conf_mat[2,]) * conf_mat[,2])
  k <- (po-pe)/(1-pe)
  return (k)
}

# get a boundary cell's neighboring cells
get_neighbor <- function(index_v, n_size, ras_ncols, ras_nrows){
  index_list <- c()
  radius <- (n_size-1)/2
  loc_y <- (index_v-1) %/% ras_ncols + 1
  loc_x <- (index_v-1) %% ras_ncols + 1
  for ( y in max(1,loc_y-radius) : min(ras_nrows, loc_y+radius ) ){
    loc <- (y-1)*ras_ncols + max(1, loc_x-radius):min(ras_ncols, loc_x+radius)
    index_list <- append(index_list, loc)
  }
  return (index_list)
}