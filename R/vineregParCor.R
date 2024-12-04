#' internal function
#' @noRd

vineregParCor <- function(df, cores_vine){
  #Initialize 
  var_indx <- 1:(ncol(df)-1) #Will keep track of which explanatory variables have not been added to the vine
  cond <- TRUE
  single_vars <- vector() #Will keep track of which explanatory variables have been added to the vine
  count <- 1
  prev_mdl <- list()
  total_vars <- length(var_indx)
  thr_caic <- 1000000
  u_vars <- matrix(0, nrow(df), ncol(df))
  
  #For every variable, do rank-based conversion to rank-associated uniform values
  for(j in 1:ncol(df)){
    u_vars[,j] <- u_scale(df[,j])
  }

  #Convert the rank-associated uniform values to rank-associated standard normal values. Q: why not go straight to rank-associated normal values?
  z_vars <- norm_score(u_vars)
  #Compute Pearson's correlation coefficient for every pair of variables
  cor_mat <- stats::cor(z_vars)

  while(cond==TRUE){
    thr_single <- -2
    df_cor <- vector()
    df_cor_test <- vector()
    #For every explanatory variable, find the (absolute value of) partial correlation of the response variable, explanatory variable given all explanatory variables already added to the D-vine
    for(i in 1:length(var_indx)){
      mdl_single <- abs(ParCor(cor_mat, single_vars+1, 1,  (1+var_indx[i]))) #Q: what does single_vars+1 do?
      if(mdl_single > thr_single){
        thr_single <- mdl_single
        var_single <- var_indx[i]
      }
    }
    single_vars <- union(single_vars, var_single)
    var_indx <- setdiff(var_indx, var_single)
    cand_vars <- df[,(1+single_vars)]
    name_cnd <- colnames(df)[(1+single_vars)]
    mdl <- update_mdl(df[,1], cand_vars, name_cnd, cores_vine)
    prev_mdl[[count]] <- mdl
    mdl_caic <- mdl$fit$stats$caic
    if(mdl_caic >= thr_caic){
      single_vars <- single_vars[-length(single_vars)]
      final_mdl <- prev_mdl[[count-1]]$fit
      cond <- FALSE
      break
    }
    if(length(single_vars) == total_vars){
      final_mdl <- prev_mdl[[count]]$fit
      cond <- FALSE
      break
    }
    thr_caic <- mdl_caic
    count <- count + 1
  }
  res_list <- list("method"="ParCor", "vinegraph"="Dvine","vars_indx"=single_vars,
                   "vinereg_fit"=final_mdl)
}






