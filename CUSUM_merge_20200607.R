#' Merging data-frames
#' Merging data-frames from separate years to form a complete merged data frame
#' @param data_ls a list of data-frames from separate years
#' @param tau0_vec a vector of beginning calendar/chronological dates of each data-frame in data_ls
#' @param name_ls a list of names for \code{id}, \code{chro_time}, \code{delta_vec}, \code{enrl_t}, and \code{x}. 
#'@return
#' A complete data list with records merged by \code{id}. 

merge.time <- function(data_ls, tau0_vec,
                       name_ls = list(
                         id_name = 'id',
                         chro_time_name = "chro_time",
                         delta_vec_name = "delta_vec",
                         enrl_t_name = "enrl_t",
                         covariates_name = 'x'
                       )){
  
  len_data <- length(data_ls)
  if(len_data > 1){
    for(ld in 1:len_data){
      
      if(ld == 1){
        data_tmp <- data_ls[[ld]]
        order_tmp <- order(data_tmp[,name_ls$id_name])
        data_out <- data_tmp[order_tmp,]
        id <- data_out[,name_ls$id_name]
        n <- nrow(data_out)
        rm(data_tmp)
      } else {
        data_tmp <- data_ls[[ld]]
        n_tmp <- nrow(data_tmp)
        for(i in 1:n_tmp){ # case of merging
          if(data_tmp[i,name_ls$id_name] %in% id){
            j <- max(which(id == data_tmp[i,name_ls$id_name]))
            if(abs(data_out[j,name_ls$chro_time_name]- 
                data_tmp[i,name_ls$enrl_t_name] - tau0_vec[ld] +  1) <= 1){
                     data_out[j,name_ls$chro_time_name] <-
                       data_tmp[i,name_ls$chro_time_name] + tau0_vec[ld] - 1
                     
                     data_out[j,name_ls$delta_vec_name] <-
                       max(data_tmp[i,name_ls$delta_vec_name],
                           data_out[j,name_ls$delta_vec_name])
                     
            } else {
              n <- n + 1
              data_out[n,] <- data_tmp[i,]
              data_out[n,name_ls$enrl_t_name] <- 
                data_out[n,name_ls$enrl_t_name] + tau0_vec[ld] - 1
              data_out[n,name_ls$chro_time_name] <-
                data_out[n,name_ls$chro_time_name] + tau0_vec[ld] - 1
              }

          } else {
            # a separate record (might from the existing subject)
            n <- n + 1
            data_out[n,] <- data_tmp[i,]
            data_out[n,name_ls$enrl_t_name] <- 
              data_out[n,name_ls$enrl_t_name] + tau0_vec[ld] - 1
            data_out[n,name_ls$chro_time_name] <-
              data_out[n,name_ls$chro_time_name] + tau0_vec[ld] - 1
            
          }
        }
        rm(data_tmp)
        id <- data_out[,name_ls$id_name]
        order_id <- order(id)
        data_out <- data_out[order_id,]
        id <- data_out$id
       # n <- nrow(data_out)
      }
    }
  }
  return(list(id = data_out[,name_ls$id_name],
              chro_time = data_out[,name_ls$chro_time_name],
              enrl_t = data_out[,name_ls$enrl_t_name],
              delta_vec = data_out[,name_ls$delta_vec_name],
              x = data_out[,name_ls$covariates]))
}

if(F){
  options(digits = 4)
  data <- data.gen(n=1000,
                   mu=0,
                   p_enroll=0.8,
                   lambda0=c(-log(1-0.15)/365),
                   beta=c(0.3,-0.3),
                   x_def = list(mean=c(0,0),var=diag(2)),
                   cen_rt=0.1,
                   tau0=1,
                   tau=365*2)
  data_df <- data.frame( id = data$id,
                         chro_time = data$chro_time,
                         enrl_t = data$enrl_t,
                         delta_vec = data$delta_vec,
                         x = data$x)
  !is.unsorted(data_df$id)
  data_split1 <- split.time(chro_time = data$chro_time,
                           delta_vec = data$delta_vec,
                           enrl_t = data$enrl_t,
                           tau0=1,
                           tau=365,
                           id=data$id,
                           x = data$x,
                           Lambda0 = data$Lambda0)
  data_split1_df <- data.frame( id = data_split1$id2,
                                chro_time = data_split1$chro_time2,
                                enrl_t = data_split1$enrl_t2,
                                delta_vec = data_split1$delta_vec2,
                                x = data_split1$x2)
  data_split2 <- split.time(chro_time = data$chro_time,
                            delta_vec = data$delta_vec,
                            enrl_t = data$enrl_t,
                            tau0=1+365,
                            tau=365*2,
                            id=data$id,
                            x = data$x,
                            Lambda0 = data$Lambda0)
  data_split2_df <- data.frame( id = data_split2$id2,
                                chro_time = data_split2$chro_time2,
                                enrl_t = data_split2$enrl_t2,
                                delta_vec = data_split2$delta_vec2,
                                x = data_split2$x2) 
  
  data_merged <- merge.time(data_ls = list(data_split1_df, data_split2_df),
                            tau0_vec = c(1, 365+1),
                            name_ls = list(
                              id_name = 'id',
                              chro_time_name = "chro_time",
                              delta_vec_name = "delta_vec",
                              enrl_t_name = "enrl_t",
                              covariates_name = c('x.1','x.2')
                            ))
  
  !is.unsorted(data_merged$id)
  
  identical(data_df$id,data_merged$id)
  
  range(data$chro_time-data_merged$chro_time)
  range(data$enrl_t-data_merged$enrl_t)
  range(data$delta_vec-data_merged$delta_vec)
  range(data$x-data_merged$x)
  
  
  data <- data.gen(n=1000,
                   mu=0,
                   p_enroll=0.8,
                   lambda0=c(-log(1-0.15)/365),
                   beta=c(0.3,-0.3),
                   x_def = list(mean=c(0,0),var=diag(2)),
                   cen_rt=0.1,
                   tau0=1,
                   tau=365*4)
  data_df <- data.frame( id = data$id,
                         chro_time = data$chro_time,
                         enrl_t = data$enrl_t,
                         delta_vec = data$delta_vec,
                         x = data$x)
  !is.unsorted(data_df$id)
  data_split1 <- split.time(chro_time = data$chro_time,
                            delta_vec = data$delta_vec,
                            enrl_t = data$enrl_t,
                            tau0=1,
                            tau=365,
                            id=data$id,
                            x = data$x,
                            Lambda0 = data$Lambda0)
  data_split1_df <- data.frame( id = data_split1$id2,
                                chro_time = data_split1$chro_time2,
                                enrl_t = data_split1$enrl_t2,
                                delta_vec = data_split1$delta_vec2,
                                x = data_split1$x2)
  data_split2 <- split.time(chro_time = data$chro_time,
                            delta_vec = data$delta_vec,
                            enrl_t = data$enrl_t,
                            tau0=1+365,
                            tau=365*2,
                            id=data$id,
                            x = data$x,
                            Lambda0 = data$Lambda0)
  data_split2_df <- data.frame( id = data_split2$id2,
                                chro_time = data_split2$chro_time2,
                                enrl_t = data_split2$enrl_t2,
                                delta_vec = data_split2$delta_vec2,
                                x = data_split2$x2) 
  
  data_split3 <- split.time(chro_time = data$chro_time,
                            delta_vec = data$delta_vec,
                            enrl_t = data$enrl_t,
                            tau0=1+365*2,
                            tau=365*3,
                            id=data$id,
                            x = data$x,
                            Lambda0 = data$Lambda0)
  data_split3_df <- data.frame( id = data_split3$id2,
                                chro_time = data_split3$chro_time2,
                                enrl_t = data_split3$enrl_t2,
                                delta_vec = data_split3$delta_vec2,
                                x = data_split3$x2) 
  
  data_split4 <- split.time(chro_time = data$chro_time,
                            delta_vec = data$delta_vec,
                            enrl_t = data$enrl_t,
                            tau0=1+365*3,
                            tau=365*4,
                            id=data$id,
                            x = data$x,
                            Lambda0 = data$Lambda0)
  data_split4_df <- data.frame( id = data_split4$id2,
                                chro_time = data_split4$chro_time2,
                                enrl_t = data_split4$enrl_t2,
                                delta_vec = data_split4$delta_vec2,
                                x = data_split4$x2) 
  
  data_merged <- merge.time(data_ls = list(data_split1_df, 
                                           data_split2_df,
                                           data_split3_df, 
                                           data_split4_df),
                            tau0_vec = c(1, 365+1,365*2+1,365*3+1),
                            name_ls = list(
                              id_name = 'id',
                              chro_time_name = "chro_time",
                              delta_vec_name = "delta_vec",
                              enrl_t_name = "enrl_t",
                              covariates_name = c('x.1','x.2')
                            ))
  
  !is.unsorted(data_merged$id)
  length(data_merged$id)
  identical(data_df$id,data_merged$id)
  
  range(data$chro_time-data_merged$chro_time)
  range(data$enrl_t-data_merged$enrl_t)
  range(data$delta_vec-data_merged$delta_vec)
  range(data$x-data_merged$x)
}