# The package for an integrative framework for O-E CUSUM
library(ggplot2)
library(grid)
library(MASS)
#library(reshape2)
# Lili Wang (lilywang@umich.edu) edited on April 22, 2020
#
#' The funciton to Generate data for CUSUM
#'
#' @param n total sample size that are at-risk during the monitoring window.
#' @param mu the log hazard ratio between the observed facilty and the national average. default is 0.
#' @param p_enroll the probablity of enrollment at tau0, default is 0.8
#' @param lambda0 a numeric vector for the baseline hazard starting from the time 0 of this window. If its length is 1, it will be a constant baseline hazard; otherwise it will be repeated by treating its length as the period; if its length is \code{trunc((tau-tau0))}, it will be exactly a complete baseline hazard for the complete follow-up time. The default value is 15% of chance to observe event per follow-up year. 
#' @param beta a vector of parameters, the default is 0;
#' @param x_def a list composed mean vector, and variance matrix (var) for covariates, will be generated following a multivariate normal distribution. 
#' @param cen_rt yearly censoring rate, the proportion of censoring within a year. The default is 10% of yealry censoring rate.
#' @param tau0 the starting day of the monitoring window, or in general, the starting point of the follow-up/observation/investigation/study. The default is day 1. 
#' @param tau the ending day of the monitoring window, or in general, the starting point of the follow-up/observation/investigation/study. The default is 4 years later. 
#' 
#' @return  
#' \item{chro_time}{a vector of chorological day of event/censored.}
#' \item{delta_vec}{a vector of event indicator, 1 for event, 0 for censoring.}
#' \item{enrl_t}{a vector of chorological day of enrollment.}
#' \item{xbeta}{a numeric vector of product of covariates and effect estimates.}
#' \item{beta}{regression parameters.}
#' \item{x}{covariate matrix.}
#' \item{Lambda0}{a numeric vector of cumulative baseline, should be corresponding to the chro_time and enrl_t.}
#' \item{tau0}{the chorologial starting day of the monitoring window.}
#' \item{tau}{the chorological ending day of the monitoring window.}
#' \item{n}{sample size.}
#' \item{cen_rt,mu,p_enroll}{same as input.}
#'  
#' @author Lili Wang
#' @export
data.gen <- function(n, mu = 0,
                   p_enroll = 0.8,
                   lambda0 = -log(1-0.15),
                   beta = 0,
                   x_def = list(mean = NULL, var = NULL),
                   cen_rt = 0.1,
                   tau0 = 1,
                   tau = 365.25*4){
  if(tau <= tau0){
    stop("tau has to be larger than tau0")
    }
  ndays <- trunc(tau - tau0 + 1)
  if(length(lambda0) == 1 && lambda0 > 0){
    Lambda0 <- (1:ndays) * lambda0
    stationary <- T
  }else if(length(lambda0) > 1 && all(lambda0 >= 0)){
    Lambda0 <- cumsum(rep_len(lambda0, ndays))
    stationary <- F
  }else{
    stop("Wrong input of lambda0, which should be a numeric vector 
         with non-negative entries and a length of at least 1")
  }
  if(beta != 0 && length(x_def$mean) > 0){
    if(is.null(x_def$var)) stop("variance of x covariates is not defined.")
    if((length(x_def$mean) != dim(x_def$var)[1]) | 
       (dim(x_def$var)[1] != dim(x_def$var)[2])) {
      stop("wrong dimension of the mean and the variance matrix.")
    }
    x <- mvrnorm(n, x_def$mean, x_def$var)
    xbeta <- x %*% beta
  }else{
    x <- NULL
    beta <- 0
    xbeta <- rep(0, n)
  }
  
  n_early_enroll <- round(p_enroll * n)
  n_late_enroll <- n - n_early_enroll
  enrl_t <- c(rep(tau0, n_early_enroll),
              ceiling(runif(n_late_enroll, min = tau0, max = tau))) # assume a uniform enrollment after the initial enrollment
  #tau0+c(rep(0,n_early_enroll),ceiling((ndays-1)*(1:(n_late_enroll))/n_late_enroll))
  e <- rexp(n)
  #chronological time below
  if(stationary == T){
    event_time <- ceiling(e / lambda0 / exp(xbeta) / exp(mu) + enrl_t)
  }else{
    # closest_t<-function(v){
    #     ret = which.min(abs(v))
    # if(ret == 1) ret=0 else if (ret == ndays) ret = ndays+1
    # return(ret)
    # }
    event_time <- ceiling(apply(as.vector(exp(xbeta)) %o% Lambda0 * exp(mu) - e,
                                1, closest_t_cpp, tau0, tau, ndays) + enrl_t)
  }
  
  cen_day_rt <- -log(1 - cen_rt) / 365.25
  cen_time <- ceiling(pmin(enrl_t + rexp(n, cen_day_rt), tau))
  
  delta_vec <- as.numeric(cen_time > event_time)
  chro_time <- pmin(event_time, cen_time)
  
  return(list(id = 1:n, chro_time = chro_time,
              delta_vec = delta_vec, enrl_t = enrl_t,
              x = x, beta = beta, xbeta = xbeta,
              Lambda0 = Lambda0, tau0 = tau0,
              tau = tau, n = n, cen_rt = cen_rt, 
              mu = mu, p_enroll = p_enroll))
}
if(F){
  options(digits = 4)
  data=data.gen(n = 1000,
                mu = 0,
                p_enroll = 0.8,
                lambda0 = -log(1 - 0.15) / 365,
                cen_rt = 0.1,
                tau0 = 1,tau = 365 * 4)
}

#' The function to calculate the cumulative observed and expected events
#' 
#' @param chro_time a vector of chorological day of event/censored.
#' @param delta_vec a vector of event indicator, 1 for event, 0 for censoring.
#' @param enrl_t a vector of chorological day of enrollment.
#' @param xbeta a numeric vector of product of covariates and effect estimates.
#' @param Lambda0 a numeric vector of cumulative baseline, should be corresponding to the chro_time and enrl_t.
#' @param tau0 the chorologial starting day of the monitoring window.  Note that the length of \code{Lambda0} could be larger than \code{tau0} and \code{tau}.
#' @param tau the chorological ending day of the monitoring window.
#' 
#' @return 
#' \item{O_t}{cumulative observed events.}
#' \item{E_t}{cumulative expected events.}
#' \item{O_E_t}{difference between O_t and E_t.}
#' \item{tau0}{as defined in the input.}
#' \item{tau}{as defined in the input.}
#' 
#' @author Lili Wang 
#' @export
OE.compute <- function(chro_time, 
                       delta_vec, 
                       enrl_t, 
                       xbeta, 
                       Lambda0, 
                       tau0,
                       tau){
  len <- length(chro_time)
  delta_vec_len <- length(delta_vec)
  enrl_t_len <- length(enrl_t)
  xbeta_len <- length(xbeta)
  if(any(len != delta_vec_len,
         delta_vec_len != enrl_t_len,
         enrl_t_len != xbeta_len)){
    stop("The length of input data, chro_time, delta_vec, 
                            enrl_t, and xbeta, should be identical.")
    }
  ndays <- trunc(tau - tau0 + 1)
  O_t_pre <- E_t_pre <- matrix(0, nrow = len, ncol = ndays)
  at_risk <- !((chro_time < tau0) | (enrl_t > tau))
  
  for (i in 1:len){
    if(at_risk[i]){
      true_in <- max(tau0,enrl_t[i]) 
      for(j in 1:ndays){
        if(j >= chro_time[i] - tau0 + 1 && delta_vec[i]){
          O_t_pre[i,j] <- 1
        }
        deduct <- ifelse(true_in <= 1, 0, Lambda0[true_in - 1])
        if(j >= chro_time[i] - tau0 + 1) {
          E_t_pre[i,j] <- (Lambda0[chro_time[i]] - deduct) * exp(xbeta[i])
        }else if(j >= true_in - tau0 + 1 && j < chro_time[i] - tau0 + 1){
          E_t_pre[i,j] <- (Lambda0[j + tau0-1]-deduct)*exp(xbeta[i])
        }
      }
    }
  }
  
  O_t <- colSums(O_t_pre)
  E_t <- colSums(E_t_pre)
  O_E_t <- O_t - E_t
  
  return(list(O_t = O_t, E_t = E_t, O_E_t = O_E_t, tau0 = tau0, tau = tau))
}
if(F){
  options(digits = 4)
  data <- data.gen(n = 1000,
                   mu = 0,
                   p_enroll = 0.8,
                   lambda0 = -log(1 - 0.15) / 365,
                   cen_rt = 0.1, tau0 = 1, tau = 365 * 4)
  out <- OE.compute(data$chro_time,
                    data$delta_vec,
                    data$enrl_t,
                    data$xbeta,
                    data$Lambda0,
                    data$tau0,
                    data$tau)
  out$O_E_t[365 * 4]
  out$O_t[365 * 4]
  out$E_t[365 * 4]
  out1 <- OE.compute(data$chro_time,
                     data$delta_vec,
                     data$enrl_t,
                     data$xbeta,
                     data$Lambda0,
                     tau0 = 1,
                     tau = 365 * 2)
  out2 <- OE.compute(data$chro_time,
                     data$delta_vec,
                     data$enrl_t,
                     data$xbeta,
                     data$Lambda0,
                     tau0 = 1+365 * 2,
                     tau = 65 *4 )
  
  (out1$O_E_t + out2$O_E_t)[365 * 2] # identical to out$O_E_t[365*4]
  (out1$O_t + out2$O_t)[365 * 2] # identical to out$O_t[365*4]
  (out1$E_t + out2$E_t)[365 * 2] # identical to out$E_t[365*4]
}
#' 
#' Compute the monitoring bands for the O-E CUSUM chart
#' 
#' @param O_t cumulative observed events.
#' @param E_t cumulative expected events.
#' @param tau0 chronological starting day of the monitoring window.
#' @param tau chronological ending day of the monitoring window. 
#' @param theta1,theta2 correspond to the two hypotheses, theta1 is positive, theta2 is negative, representing the log hazard ratios. 
#' @param h1,h2 the two threchold values for the two hypothesis tests, h1 is for "worse than expected" test, h2 is for "better than expected" tests. 
#' @param rho1,rho2 numeric value in \eqn{[0,1)} for head start,if 0, its default is 0, implies restarting from 0; if it is NULL, it will continue without any restarting mechanism.
#' 
#' O_t,E_t,O_E_t,M1,M2,cross1,cross2,cross1_t,cross2_t
#' tau0,tau,theta1,theta2,h1,h2,rho1,rho2
#' @return 
#' \item{O_t,E_t,tau0,tau,theta1,theta2,h1,h2,rho1,rho2}{input data.}
#' \item{O_E_t}{difference between O_t and E_t.}
#' \item{M1, M2}{monitoring bands for the two hypotheses.}
#' \item{cross1, cross2}{number of signals/rejections.}
#' \item{cross1_t,cross2_t}{signal time.}
#' @export
O_E_CUSUM <- function(O_t, 
                      E_t, 
                      theta1, 
                      theta2, 
                      h1, h2, 
                      tau0, 
                      tau, 
                      rho1 = 0, 
                      rho2 = 0){
  k1 <- (exp(theta1) - 1) / theta1 - 1
  k2 <- (exp(theta2) - 1) / theta2 - 1
  
  ndays <- trunc(tau - tau0 + 1)
  O_E_t <- O_t - E_t
  M1_pre <- O_E_t - k1 * E_t #use for M_t_restart
  M2_pre <- -O_E_t + k2 * E_t
 
  start1 <- 1
  cross1 <- 0
  cross1_t <- NULL
  M1_pt1 <- NULL
  M1_restart <- NULL
  M1_min <- 0
  
  if(!is.null(rho1)){
    for(j in 1:ndays){
      M1_min <- min(M1_min, M1_pre[j])
      M1_pt1 <- M1_min - M1_pre[j] + h1
      
      if(cross1 > 0){
        M1_pt2 <- -(M1_pre[j] - M1_pre[start1]) + (1 - rho1) * h1
        M1_restart[j] <- min(M1_pt1, M1_pt2)
      }else{
        M1_restart[j] <- M1_pt1
      }
      # Record signals
      if(M1_restart[j] < 0){
        cross1 <- cross1 + 1
        cross1_t[cross1] <- j
        start1 <- j + 1
        M1_min <- M1_pre[j + 1]
      }
    }
    M1 <- M1_restart
  }else{
    M1_no_restart <- cummin(M1_pre) - M1_pre + h1
    for(j in 1:ndays){
     if(M1_no_restart[j] < 0 & cross1 == 0){
       cross1 <- 1
       cross1_t[1] <- j
     }
    }
    M1 <- M1_no_restart
  }
  
  M2_pt1 <- NULL
  start2 <- 1
  cross2 <- 0
  cross2_t <- NULL
  M2_restart <- NULL
  M2_min <- 0
  if(!is.null(rho2)){
    for(j in 1:ndays){
      M2_min <- min(M2_min, M2_pre[j])
      M2_pt1 <- M2_min - M2_pre[j] + h2
      
      if(cross2>0){
        M2_pt2 <- -(M2_pre[j] - M2_pre[start2]) + (1-rho2) * h2
        M2_restart[j] <- min(M2_pt1, M2_pt2)
      }else{
        M2_restart[j] <- M2_pt1
      }
      # Record signals
      if(M2_restart[j] < 0){
        cross2 <- cross2 + 1
        cross2_t[cross2] <- j
        start2 <- j+1
        M2_min <- M2_pre[j + 1]
      }
    }
    M2 <- M2_restart
  }else{
    M2_no_restart <- cummin(M2_pre) - M2_pre + h2
    for(j in 1:ndays){
      if(M2_no_restart[j] < 0 & cross2 == 0){
        cross2 <- 1
        cross2_t[1] <- j
      }
    }
    M2 <- M2_no_restart
  }
  return(list(O_t = O_t,
              E_t = E_t,
              O_E_t = O_E_t,
              M1 = M1, 
              M2 = M2,
              cross1 = cross1,
              cross2 = cross2,
              cross1_t = cross1_t,
              cross2_t = cross2_t,
              tau0 = tau0,
              tau = tau,
              theta1 = theta1,
              theta2 = theta2,
              h1 = h1,
              h2 = h2,
              rho1 = rho1,
              rho2 = rho2))
}

if(F){
  options(digits = 4)
  data <- data.gen(n = 1000,
                   mu = 0, 
                   p_enroll = 0.8,
                   lambda0 = -log(1 - 0.15) / 365,
                   cen_rt = 0.1,
                   tau0 = 1, tau = 365 * 4)
  out <- OE.compute(data$chro_time,
                    data$delta_vec,
                    data$enrl_t,
                    data$xbeta,
                    data$Lambda0,
                    data$tau0,
                    data$tau)
  out_CUSUM <- O_E_CUSUM(out$O_t,
                         out$E_t,
                         theta1 = log(2),
                         theta2 = log(0.5),
                         h1 = 7,
                         h2 = 7,
                         tau0 = 1,
                         tau = 365 * 4,
                         rho1 = 0, 
                         rho2 = 0)
}

#' Standardized Tabular CUSUM
#' 
#' @param O_t cumulative observed events.
#' @param E_t cumulative expected events.
#' @param tau0 chronological starting day of the monitoring window.
#' @param tau chronological ending day of the monitoring window. 
#' @param theta log hazard ratio of the alternative hypothesis for testing. 
#' @param h  the threchold value for the hypothesis test
#' @param rho numeric value in \eqn{[0,1)} for head start,if 0, its default is 0, implies restarting from 0; if it is NULL, it will continue without any restarting mechanism.
#' 
#' @return 
#' \item{O_t,E_t,tau0,tau,theta,h,rho}{input data.}
#' \item{O_E_t}{difference between O_t and E_t.}
#' \item{G}{monitoring bands for the hypothesis test.}
#' \item{cross}{number of signals/rejections.}
#' \item{cross_t}{signal time.}
#' \item{L}{the threshold value for testing. }

#' @export
tabular_CUSUM <- function(O_t,
                          E_t,
                          theta,
                          h,
                          tau0,
                          tau,
                          rho = 0){
  L <- h * abs(theta - 0)
  ndays <- trunc(tau - tau0 + 1)
  cross <- 0
  cross_t <- NULL
  U_pre <- (theta - 0) * O_t-(exp(theta) - 1) * E_t
  G_pre <- U_pre - cummin(U_pre)
  U_min <- 0
  start <- 1
  if(!is.null(rho)){
    G_restart <- NULL
    for(j in 1:ndays){
      U_min <- min(U_min, U_pre[j])
      U_pt1 <- U_pre[j] - U_min
      if(cross > 0){
        U_pt2 <- U_pre[j] - U_pre[start] + L * rho
        G_restart[j] <- max(U_pt1, U_pt2)
      }else{
        G_restart[j] <- U_pt1
      }
      
      if(G_restart[j] > L){
        cross = cross + 1
        cross_t[cross] = j
        start = j + 1
        U_min <- U_pre[j + 1]
      }
    }
    G <- G_restart
  }else{

    for(j in 1:ndays){
      if(G_pre[j] > L & cross == 0){
        cross = 1
        cross_t[1] = j
      }
    }
    G <- G_pre
  }
  return(list(O_t = O_t,
              E_t = E_t,
              G = G,
              cross = cross,
              cross_t = cross_t,
              L = L,
              tau0 = tau0,
              tau = tau,
              theta = theta,
              h = h, rho = rho))
  
}
if(F){
  options(digits = 4)
  data <- data.gen(n = 1000, 
                   mu = 0,
                   p_enroll = 0.8,
                   lambda0 = -log(1 - 0.15)/365,
                   cen_rt = 0.1,
                   tau0 = 1, tau = 365 * 4)
  out <- OE.compute(data$chro_time,
                    data$delta_vec,
                    data$enrl_t,
                    data$xbeta,
                    data$Lambda0,
                    data$tau0,
                    data$tau)
  out_CUSUM<-O_E_CUSUM(out$O_t,
                       out$E_t,
                       theta1 = log(2),
                       theta2 = log(0.5),
                       h1 = 4, h2 = 4,
                       tau0 = 1, tau = 365 * 4,
                       rho1 = 0, rho2 = 0)
  out_CUSUM$cross1
  out_CUSUM$cross1_t
  out_CUSUM$cross2
  out_CUSUM$cross2_t
  out_tabular1 <- tabular_CUSUM(out$O_t,
                                out$E_t,
                                theta = log(2),
                                h = 4,
                                tau0 = 1,
                                tau = 365 * 4,
                                rho = 0)
  out_tabular1$cross
  out_tabular1$cross_t
  out_tabular2<-tabular_CUSUM(out$O_t,
                              out$E_t,
                              theta = -log(2),
                              h = 4,
                              tau0 = 1,
                              tau = 365 * 4,
                              rho = 0)
  out_tabular2$cross
  out_tabular2$cross_t
}
#' Plot the O-E CUSUM
#' 
#' Plot the O-E CUSUM object.
#' 
#' @param result the output object from \code{O_E_CUSUM}.
#' @param title the title of the plot, the default is \code{character(0)}.
#' 
#' @return
#' \item{result_plot}{the ggplot object.}
#' \item{result}{a copy of the input result. }
#' 
#' @export
O_E_CUSUM_plot <- function(result, title = character(0)){
  
  theta1 <- result$theta1
  theta2 <- result$theta2


  data_plot <- data.frame(time = (result$tau0):(result$tau),
                          O_E_t = result$O_E_t,
                          M1 = result$M1 + result$O_E_t,
                          M2 = -result$M2 + result$O_E_t)
  
  data_plot2 = reshape2::melt(data_plot,id="time")
  data_plot2$variable = factor(data_plot2$variable,
                             levels = c("O_E_t", "M1", "M2"),
                             labels = c("O-E", "Upper MB", "Lower MB"))
  
  plot_range <- range(data_plot2$value)
  plot_range <-  plot_range + c(-0.05,+0.05)*(plot_range[2]- plot_range[1])
  
  
  result_plot <- ggplot(data = data_plot2,
                        aes(x = time, y = value,
                             color = variable, linetype = variable))+
                geom_line(size = 0.3) + 
                ylab("Excess Deaths (O-E)") + 
                scale_color_manual(values = c(1,2,4)) + 
                scale_linetype_manual(values = c(1,2,2)) + 
                ylim(plot_range[1], plot_range[2]) +
                theme_light() + ggtitle(title) + 
                theme(plot.title = element_text(hjust = 0.5))
  
    if(result$cross1 > 0)
      result_plot = result_plot+
                    geom_vline(xintercept = result$cross1_t,
                               size = 0.5,color = 2,linetype = 2,alpha = 0.5) +
                    annotate("text", x = result$cross1_t + 15, y = plot_range[2],
                             color = 2, label = as.character(result$cross1_t))
                   
    if(result$cross2 > 0)
      result_plot = result_plot+
                    geom_vline(xintercept = result$cross2_t,
                               size = 0.5, color = 4, linetype = 2, alpha = 0.5)+
                    annotate("text", x = result$cross2_t + 15,
                             y = plot_range[1], color = 4,
                             label = as.character(result$cross2_t))
  
  return(list( result_plot = result_plot, result = result))
}
if(F){
  options(digits = 4)
  data <- data.gen(n = 500, mu = 0, p_enroll = 0.6,
                   lambda0 = -log(1 - 0.15) / 365,
                   cen_rt = 0.1, tau0 = 1, tau = 365 * 4)
  out <- OE.compute(data$chro_time,
                    data$delta_vec,
                    data$enrl_t,
                    data$xbeta,
                    data$Lambda0,
                    data$tau0,
                    data$tau)
  
  out_CUSUM <- O_E_CUSUM(out$O_t,
                         out$E_t,
                         theta1 = log(2),
                         theta2 = log(0.5),
                         h1 = 6, h2 = 6, 
                         tau0 = 1,tau=365 * 4,
                         rho1 = 0.5, rho2 = 0.5)
  plot_out <- O_E_CUSUM_plot(out_CUSUM)
  plot_out$result_plot
}
#' 
#' Plot Tabular CUSUM
#' 
#' The tool to plot standard Tabular CUSUM 
#' 
#' @param result the output object from \code{tabular_CUSUM}.
#' @param title the title of the plot, the default is \code{character(0)}.
#' 
#' @return
#' \item{result_plot}{the ggplot object.}
#' \item{result}{a copy of the input result. }
#' 
#' @export
tabular_CUSUM_plot <- function(result, title = character(0)){

  data_plot <- data.frame(time = (result$tau0):(result$tau),G = result$G)
  clr <- ifelse(result$theta > 0, 2, 4)
  ylab.txt <- ifelse(result$theta>0, "Worse than Expected", "Better than Expected")
  result_plot <- ggplot(data = data_plot, aes(x = time, y = G))+
                 geom_line(size=0.3)+
                 ggtitle(title) + 
                 ylab(ylab.txt)+
                 theme_light() + 
                 theme(plot.title = element_text(hjust = 0.5, size = 15)) +
                 geom_hline(yintercept = result$L,color = clr, linetype = "dashed")
  
  if(result$theta<0) result_plot <- result_plot+scale_y_reverse()
    
  if(result$cross>0) {
    result_plot = result_plot+ 
                  geom_vline(xintercept = result$cross_t,
                             color = clr, linetype = "dashed",
                             size = 0.5, alpha = 0.5) +
                  annotate("text", x = result$cross_t + 15,
                           y= -0.1, color = clr, 
                           label = as.character(result$cross_t))
  }
      
  return(list(result_plot = result_plot, result = result))
}
if(F){
  options(digits = 4)
  data <- data.gen(n = 1000, 
                   mu = log(1.2), 
                   p_enroll = 0.6,
                   lambda0 = -log(1 - 0.15) / 365,
                   cen_rt = 0.1,
                   tau0 = 1,
                   tau = 365 * 4)
  out <- OE.compute(data$chro_time,
                    data$delta_vec,
                    data$enrl_t,
                    data$xbeta,
                    data$Lambda0,
                    data$tau0,
                    data$tau)
  out_CUSUM <- O_E_CUSUM(out$O_t,
                         out$E_t,
                         theta1 = log(1.2),
                         theta2 = log(0.8),
                         h1 = 25.54, 
                         h2 = 22.72,
                         tau0 = 1,
                         tau = 365 * 4,
                         rho1 = 0.5,
                         rho2 = 0.5)
  plot_out <- O_E_CUSUM_plot(out_CUSUM)
  plot_out$result_plot
  out_tabular1<-tabular_CUSUM(out$O_t,
                              out$E_t,
                              theta = log(1.2),
                              h = 25.54,
                              tau0 = 1,
                              tau = 365 * 4,
                              rho = 0.5) 
  plot_out1 <- tabular_CUSUM_plot(out_tabular1)
  plot_out1$result_plot
  out_tabular2<-tabular_CUSUM(out$O_t, 
                              out$E_t,
                              theta = log(0.8),
                              h = 22.72,
                              tau0 = 1,
                              tau = 365 * 4,
                              rho = 0.5)
  plot_out2 <- tabular_CUSUM_plot(out_tabular2)
  plot_out2$result_plot
  
}
#' Split data into smaller pieces
#' 
#' Split the full follow-up into separate smaller interval of follow
#' 
#' @inheritParams OE.compute
#' @param id  identification numbers, will be added if not given
#' @param x covariates, default is null. 
#' 
#' @return 
#' The version splits data by time interval (tau0,tau]. 
#' 
split.time <- function(chro_time,
                       delta_vec,
                       enrl_t,tau0,
                       tau,
                       id = NULL,
                       x = NULL,
                       xbeta = NULL,
                       Lambda0 = NULL){
  len <-length(chro_time)
  if(is.null(id)) id <- 1:len
  delta_vec_len <- length(delta_vec)
  enrl_t_len <- length(enrl_t)
  if(any(len != delta_vec_len, delta_vec_len != enrl_t_len)){
    stop("The length of input data, chro_time, delta_vec, 
         enrl_t, and xbeta, should be identical.")}
  at_risk <- !((chro_time < tau0) | (enrl_t > tau))
  enrl_t2 <- pmax(tau0, enrl_t[at_risk]) - tau0 + 1
  chro_time2 <- pmin(tau, chro_time[at_risk]) - tau0 + 1
  delta_vec2 <- delta_vec[at_risk] * (chro_time[at_risk] <= tau)
  if(!is.null(Lambda0)) {
    if(tau0 <= 1) Lambda02 <- Lambda0[tau0:tau] else 
      Lambda02 <- Lambda0[tau0:tau] - Lambda0[tau0 - 1]
    } else Lambda02=NULL
  
  if(is.null(x)) x2 <- NULL else x2 <- x[at_risk, ]
  if(is.null(xbeta)) xbeta2 <- NULL else xbeta2 <- xbeta[at_risk]
  id2 <- id[at_risk]
  return(list(chro_time2 = chro_time2,
              delta_vec2 = delta_vec2,
              enrl_t2 = enrl_t2,
              xbeta2 = xbeta2,
              Lambda02 = Lambda02,
              tau0, tau,
              id2 = id2, x2 = x2))
}
#
if(F){
  options(digits = 4)
  data <- data.gen(n = 1000, 
                   mu = 0,
                   p_enroll = 0.8,
                   lambda0 = -log(1 - 0.15) / 365,
                   cen_rt = 0.1,
                   tau0 = 1,
                   tau = 365 * 4)
  out <- OE.compute(data$chro_time,
                    data$delta_vec,
                    data$enrl_t,
                    data$xbeta,
                    data$Lambda0,
                    data$tau0,
                    data$tau)
  out$O_E_t[365 * 4]
  out$O_t[365 * 4]
  out$E_t[365 * 4]
  out1 <- OE.compute(data$chro_time,
                     data$delta_vec,
                     data$enrl_t,
                     data$xbeta,
                     data$Lambda0,
                     tau0 = 1, tau = 365 * 2)
  out2 <- OE.compute(data$chro_time,
                     data$delta_vec,
                     data$enrl_t,
                     data$xbeta,
                     data$Lambda0,
                     tau0=1 + 365 * 2,
                     tau = 365 * 4)
  (out1$O_E_t + out2$O_E_t)[365 * 2] # identical to out$O_E_t[365*4]
  (out1$O_t + out2$O_t)[365 * 2] # identical to out$O_t[365*4]
  (out1$E_t + out2$E_t)[365 * 2] # identical to out$E_t[365*4]
  
  data_split <- split.time(data$chro_time,
                           data$delta_vec,
                           data$enrl_t,
                           data$Lambda0,
                           tau0=1,
                           tau=365 * 2,
                           xbeta = data$xbeta,
                           Lambda0 = data$Lambda0)
  out1_again <- OE.compute(data_split$chro_time2,
                           data_split$delta_vec2,
                           data_split$enrl_t2,
                           data_split$xbeta2,
                           data_split$Lambda02,
                           tau0 = 1,tau = 365 * 2)
  identical(out1_again, out1) # TRUE!
  
  data_split2 <- split.time(data$chro_time,
                            data$delta_vec,
                            data$enrl_t,
                            tau0 = 1 + 365 * 2,
                            tau = 365 * 4,
                            xbeta = data$xbeta,
                            Lambda0 = data$Lambda0)
  out2_again <- OE.compute(data_split2$chro_time2,
                           data_split2$delta_vec2,
                           data_split2$enrl_t2,
                           data_split2$xbeta2,
                           data_split2$Lambda02,
                           tau0 = 1,
                           tau = 365 * 2)
  identical(out2_again$O_t, out2$O_t) # TRUE!
  identical(out2_again$E_t, out2$E_t) # TRUE!
  # Note that their tau0 and tau are different
}

#' h values
#' 
#' Selet the control limit h to control the type I error
#' 
#' @inheritParams data.gen
#' 

O_E_CUSUM_hval <- function(nloop,
                           n, 
                           theta1, 
                           theta2,
                           mu = 0,
                           p_enroll = 0.8,
                           lambda0 = -log(1 - 0.15) / 365,
                           beta = 0,
                           x_def = list(mean = NULL, var = NULL),
                           cen_rt = 0.1, 
                           tau0 = 1,
                           tau = 365 * 4, 
                           alpha = 0.95){
  k1 <- (exp(theta1) - 1) / theta1 - 1
  k2 <- (exp(theta2) - 1) / theta2 - 1
  
  fun <- function(n, mu, p_enroll, lambda0, xbeta, cen_rt, tau0, tau){
    dt_tmp <- data.gen(n = n,
                       mu = mu,
                       p_enroll = p_enroll,
                       lambda0 = lambda0,
                       beta = beta,
                       x_def = x_def,
                       cen_rt = cen_rt,
                       tau0 = tau0, 
                       tau = tau)
    OE_tmp <- OE.compute(dt_tmp$chro_time,
                         dt_tmp$delta_vec,
                         dt_tmp$enrl_t,
                         dt_tmp$xbeta,
                         dt_tmp$Lambda0,
                         dt_tmp$tau0,
                         dt_tmp$tau)
    #print(OE_tmp$O_E_t[365*4])
    M1_pre <- OE_tmp$O_E_t - k1 * OE_tmp$E_t 
    M2_pre <- -OE_tmp$O_E_t + k2 * OE_tmp$E_t
    h1<- max(M1_pre - cummin(M1_pre))
    h2<- max(M2_pre - cummin(M2_pre)) 
    return(c(h1, h2))
  }
  h <- matrix(nrow = nloop, ncol = 2)
  pb <- txtProgressBar(min = 0, max = 1, style = 3)
  for(l in 1:nloop){
    h[l, ]<- fun(n, mu, p_enroll, lambda0, xbeta, cen_rt, tau0, tau)
    setTxtProgressBar(pb, l / nloop)
  }
  close(pb)
  h_select <- apply(h, 2, quantile, probs = 0.95)
  return(h_select)
}
if(F){
 h <- O_E_CUSUM_hval(100, 
                     n = 1000, 
                     theta1 = log(2), 
                     theta2 = -log(2), 
                     mu = 0, 
                     p_enroll = 0.8,
                     lambda0 = -log(1 - 0.15) / 365,
                     cen_rt = 0.1,
                     tau0 = 1,
                     tau = 365 * 4, 
                     alpha = 0.95)
 # 10.92 10.07 if nloop=100 #10.3 9.945  if nloop=1000 # 10.29 9.811 if nloop=10000
 h
 h <- O_E_CUSUM_hval_cpp(10000, 1000, c(log(2), -log(2)))
 O_E_CUSUM_sim_cpp(10000, 1000, c(log(2),-log(2)), c(0), 0, h)
 
}



### Utility functions

combine_rec <- function(rec_list){
  len <-length(rec_list)
  if (len > 0) out <- rec_list[[1]] else stop("empty list")
  if(len > 2) {
    for(l in 2:len){
    out <-c(out, out[length(out)] + rec_list[[l]])
     }
  }
  return(out)
}

if(F){
  setwd("/Users/liliwang/Box/Cusum_simulation/restart3/For_paper/results/20200422/")
  ### Merged data
  set.seed(2020)
  library(survival)
  options(digits = 4)
  tau0 <- 1
  tau <- 365 * 4
  theta1 <- log(1.2)
  theta2 <- log(0.8)
  rho1 <- 0
  rho2 <- 0
 lambda0 = c(0, diff(exp(0.0004*(1:480)))) # -log(1-0.15)/365 c(0,diff(exp(0.0005*(1:365))))# c(0,diff(exp(0.0004*(1:480))))
  job_type <- "aperiodic.nonstat" # "stat" "periodic.nonstat" "aperiodic.nonstat"
  
  data_control_merged <- data.gen(n=100000, mu=0, p_enroll=0.8,lambda0=lambda0,beta=c(0.3,-0.3),x_def = list(mean=c(0,0),var=diag(2)),cen_rt=0.1,tau0=tau0,tau=tau)
  
  (model_control_merged <- coxph(Surv(time=data_control_merged$enrl_t,time2=data_control_merged$chro_time,event=data_control_merged$delta_vec)~data_control_merged$x))
  
  Baseline_hat_merged <- basehaz(model_control_merged,centered=F)
  Lambda0_hat_merged <- with(Baseline_hat_merged,rep(c(0,hazard),diff(c(tau0,time,tau+1))))
  lambda0_hat <- c(0,diff(Lambda0_hat_merged))
  
  xbeta_hat_control <- data_control_merged$x%*%as.vector(coef(model_control_merged))
  
  
  # Semiparametric
  #h_cpp <- c(24.78,22.31)
  h_cpp <- O_E_CUSUM_hval_cpp(10000,n=1000,c(theta1,theta2),
                              p_enroll=0.8,lambda0=lambda0_hat,
                              xbeta=xbeta_hat_control)
  
  # h_cpp <- O_E_CUSUM_hval_resample_cpp(10000,1000, c(theta1,theta2),
  #                                               data_control_merged$enrl_t,
  #                                               data_control_merged$chro_time,
  #                                               data_control_merged$delta_vec,
  #                                              xbeta=xbeta_hat_control, 
  #                                              lambda0=lambda0_hat)
  # 
  # O_E_CUSUM_sim_cpp(10000,1000,c(theta1,theta2),
  #                   p_enroll=0.8,
  #                   lambda0=lambda0_hat,
  #                   xbeta=xbeta_hat_control,rho_vec = c(0),
  #                   mu=0,h=h_cpp,
  #                   tau0=tau0,
  #                   tau=tau)
  # 
  O_E_CUSUM_sim_resample2_cpp(10000,1000, c(theta1,theta2),
                             rho_vec = c(0,0.25,0.5),
                             h=h_cpp,
                             data_control_merged$enrl_t,
                             data_control_merged$chro_time,
                             data_control_merged$delta_vec,
                             xbeta=xbeta_hat_control,
                             lambda0=lambda0_hat,
                             tau0=tau0,
                             tau=tau)

  
  mu_type <- "H0" #"H0" "H1", "H2"
  mu=log(1) #log(1) log(1.2) log(0.8)
  
  data_sample_merged <- data.gen(n=1000,mu=mu,p_enroll=0.8,lambda0=lambda0,beta=c(0.3,-0.3),x_def = list(mean=c(0,0),var=diag(2)),cen_rt=0.1,tau0=tau0,tau=tau)
  xbeta_hat_merged <- data_sample_merged$x%*%as.vector(coef(model_control_merged))
  
  OE_merged <- OE.compute(data_sample_merged$chro_time,data_sample_merged$delta_vec,data_sample_merged$enrl_t,xbeta_hat_merged,Lambda0_hat_merged,tau0,tau)
  #h <- O_E_CUSUM_hval(100,1000,theta1=theta1,theta2=theta2)
  CUSUM_merged<-O_E_CUSUM(O_t=OE_merged$O_t,E_t=OE_merged$E_t,theta1=theta1,theta2=theta2,h1=h_cpp[1],h2=h_cpp[2],tau0=tau0,tau=tau,rho1=rho1,rho2=rho2)
  plot_merged <- O_E_CUSUM_plot(CUSUM_merged,title="merged")
  #tabular_merged<-tabular_CUSUM(O_t=OE_merged$O_t,E_t=OE_merged$E_t,theta=theta1,h=h_cpp[1],tau0=tau0,tau=tau,rho=rho1)
  #tabular_CUSUM_plot(tabular_merged)
  res=plot_merged$result_plot
  ggsave(filename=paste(job_type,mu_type,"merged.pdf",sep="_"),width = 6, height=4.5)
  save(res,file = paste(job_type,mu_type,"merged.RData",sep="_"))
}


if(F){
  
  ### Splitted ####
  data_control_splitted1 <- split.time( chro_time=data_control_merged$chro_time, delta_vec=data_control_merged$delta_vec, enrl_t=data_control_merged$enrl_t,tau0=1,tau=365,x=data_control_merged$x)
  data_control_splitted2 <- split.time( chro_time=data_control_merged$chro_time, delta_vec=data_control_merged$delta_vec, enrl_t=data_control_merged$enrl_t,tau0=1+365,tau=365*2,x=data_control_merged$x) 
  
  data_control_splitted3 <- split.time( chro_time=data_control_merged$chro_time, delta_vec=data_control_merged$delta_vec, enrl_t=data_control_merged$enrl_t,tau0=1+365*2,tau=365*3,x=data_control_merged$x)
  data_control_splitted4 <- split.time( chro_time=data_control_merged$chro_time, delta_vec=data_control_merged$delta_vec, enrl_t=data_control_merged$enrl_t,tau0=1+365*3,tau=365*4,x=data_control_merged$x) 
  data_control_stacked <- list(
    id=c(data_control_splitted1$id2,data_control_splitted2$id2,data_control_splitted3$id2,data_control_splitted4$id2),
    chro_time=c(data_control_splitted1$chro_time2,data_control_splitted2$chro_time2,data_control_splitted3$chro_time2,data_control_splitted4$chro_time2),
    enrl_t=c(data_control_splitted1$enrl_t2,data_control_splitted2$enrl_t2,data_control_splitted3$enrl_t2,data_control_splitted4$enrl_t2),
    delta_vec=c(data_control_splitted1$delta_vec2,data_control_splitted2$delta_vec2,data_control_splitted3$delta_vec2,data_control_splitted4$delta_vec2),
    x=rbind(data_control_splitted1$x2,data_control_splitted2$x2,data_control_splitted3$x2,data_control_splitted4$x2))
  (model_control_stacked <- coxph(Surv(time=data_control_stacked$enrl_t,time2=data_control_stacked$chro_time,event=data_control_stacked$delta_vec)~data_control_stacked$x))
  Baseline_hat_stacked <- basehaz(model_control_stacked,centered=F)
  Lambda0_hat_stacked <- with(Baseline_hat_stacked,rep(c(0,hazard),diff(c(tau0,time,365+1))))
   ### Splitted sample
  
  data_sample_splitted1 <- split.time( chro_time=data_sample_merged$chro_time, delta_vec=data_sample_merged$delta_vec, enrl_t=data_sample_merged$enrl_t,tau0=1,tau=365,x=data_sample_merged$x)
  
  data_sample_splitted2 <- split.time( chro_time=data_sample_merged$chro_time, delta_vec=data_sample_merged$delta_vec, enrl_t=data_sample_merged$enrl_t,tau0=1+365,tau=365*2,x=data_sample_merged$x) 
  
  data_sample_splitted3 <- split.time( chro_time=data_sample_merged$chro_time, delta_vec=data_sample_merged$delta_vec, enrl_t=data_sample_merged$enrl_t,tau0=1+365*2,tau=365*3,x=data_sample_merged$x)
  data_sample_splitted4 <- split.time( chro_time=data_sample_merged$chro_time, delta_vec=data_sample_merged$delta_vec, enrl_t=data_sample_merged$enrl_t,tau0=1+365*3,tau=365*4,x=data_sample_merged$x) 

  
  xbeta_hat_stacked1 <- data_sample_splitted1$x2%*%as.vector(coef(model_control_stacked))
  xbeta_hat_stacked2 <- data_sample_splitted2$x2%*%as.vector(coef(model_control_stacked))
  xbeta_hat_stacked3 <- data_sample_splitted3$x2%*%as.vector(coef(model_control_stacked))
  xbeta_hat_stacked4 <- data_sample_splitted4$x2%*%as.vector(coef(model_control_stacked))
  
  
  OE_stacked1 <- OE.compute(data_sample_splitted1$chro_time,data_sample_splitted1$delta_vec,data_sample_splitted1$enrl_t,xbeta_hat_stacked1,Lambda0_hat_stacked,1,365)
  OE_stacked2 <- OE.compute(data_sample_splitted2$chro_time,data_sample_splitted2$delta_vec,data_sample_splitted2$enrl_t,xbeta_hat_stacked2,Lambda0_hat_stacked,1,365)
  OE_stacked3 <- OE.compute(data_sample_splitted3$chro_time,data_sample_splitted3$delta_vec,data_sample_splitted3$enrl_t,xbeta_hat_stacked3,Lambda0_hat_stacked,1,365)
  OE_stacked4 <- OE.compute(data_sample_splitted4$chro_time,data_sample_splitted4$delta_vec,data_sample_splitted4$enrl_t,xbeta_hat_stacked4,Lambda0_hat_stacked,1,365)
  
  OE_stacked <- list(
    O_t= combine_rec(list(OE_stacked1$O_t,OE_stacked2$O_t,
                          OE_stacked3$O_t,OE_stacked4$O_t)),
    E_t= combine_rec(list(OE_stacked1$E_t,OE_stacked2$E_t,
                          OE_stacked3$E_t,OE_stacked4$E_t)),
    O_E_t= combine_rec(list(OE_stacked1$O_E_t,OE_stacked2$O_E_t,
                            OE_stacked3$O_E_t,OE_stacked4$O_E_t)),
    tau0=tau0,
    tau=tau
  )
  
  CUSUM_stacked<-O_E_CUSUM(O_t=OE_stacked$O_t,E_t=OE_stacked$E_t,theta1=theta1,theta2=theta2,h1=h_cpp[1],h2=h_cpp[2],tau0=tau0,tau=tau,rho1=rho1,rho2=rho2)
  
  plot_stacked <- O_E_CUSUM_plot(CUSUM_stacked,title = "stacked")
  res=plot_stacked$result_plot
  ggsave(filename=paste(job_type,mu_type,"stacked.pdf",sep="_"),width = 6, height=4.5)
  save(res,file = paste(job_type,mu_type,"stacked.RData",sep="_"))
  
}

if(F){
  ##### Separate ####
  (model_control_separate1 <- coxph(Surv(time=data_control_splitted1$enrl_t,time2=data_control_splitted1$chro_time,event=data_control_splitted1$delta_vec)~data_control_splitted1$x2))
  Baseline_hat_separate1 <- basehaz(model_control_separate1,centered=F)
  Lambda0_hat_separate1 <- with(Baseline_hat_separate1,rep(c(0,hazard),diff(c(tau0,time,365+1))))
  
  (model_control_separate2 <- coxph(Surv(time=data_control_splitted2$enrl_t,time2=data_control_splitted2$chro_time,event=data_control_splitted2$delta_vec)~data_control_splitted2$x2))
  Baseline_hat_separate2 <- basehaz(model_control_separate2,centered=F)
  Lambda0_hat_separate2 <- with(Baseline_hat_separate2,rep(c(0,hazard),diff(c(tau0,time,365+1))))
  
  (model_control_separate3 <- coxph(Surv(time=data_control_splitted3$enrl_t,time2=data_control_splitted3$chro_time,event=data_control_splitted3$delta_vec)~data_control_splitted3$x2))
  Baseline_hat_separate3 <- basehaz(model_control_separate3,centered=F)
  Lambda0_hat_separate3 <- with(Baseline_hat_separate3,rep(c(0,hazard),diff(c(tau0,time,365+1))))
  
  (model_control_separate4 <- coxph(Surv(time=data_control_splitted4$enrl_t,time2=data_control_splitted4$chro_time,event=data_control_splitted4$delta_vec)~data_control_splitted4$x2))
  Baseline_hat_separate4 <- basehaz(model_control_separate4,centered=F)
  Lambda0_hat_separate4 <- with(Baseline_hat_separate4,rep(c(0,hazard),diff(c(tau0,time,365+1))))
  
  
  xbeta_hat_separate1 <- data_sample_splitted1$x2%*%as.vector(coef(model_control_separate1))
  xbeta_hat_separate2 <- data_sample_splitted2$x2%*%as.vector(coef(model_control_separate2))
  xbeta_hat_separate3 <- data_sample_splitted3$x2%*%as.vector(coef(model_control_separate3))
  xbeta_hat_separate4 <- data_sample_splitted4$x2%*%as.vector(coef(model_control_separate4))
  
  OE_separate1 <- OE.compute(data_sample_splitted1$chro_time,data_sample_splitted1$delta_vec,data_sample_splitted1$enrl_t,xbeta_hat_separate1,Lambda0_hat_separate1,1,365)
  OE_separate2 <- OE.compute(data_sample_splitted2$chro_time,data_sample_splitted2$delta_vec,data_sample_splitted2$enrl_t,xbeta_hat_separate2,Lambda0_hat_separate2,1,365)
  OE_separate3 <- OE.compute(data_sample_splitted3$chro_time,data_sample_splitted3$delta_vec,data_sample_splitted3$enrl_t,xbeta_hat_separate3,Lambda0_hat_separate3,1,365)
  OE_separate4 <- OE.compute(data_sample_splitted4$chro_time,data_sample_splitted4$delta_vec,data_sample_splitted4$enrl_t,xbeta_hat_separate4,Lambda0_hat_separate4,1,365)
  
  OE_separate <- list(
    O_t= combine_rec(list(OE_separate1$O_t,OE_separate2$O_t,
                          OE_separate3$O_t,OE_separate4$O_t)),
    E_t= combine_rec(list(OE_separate1$E_t,OE_separate2$E_t,
                          OE_separate3$E_t,OE_separate4$E_t)),
    O_E_t= combine_rec(list(OE_separate1$O_E_t,OE_separate2$O_E_t,
                            OE_separate3$O_E_t,OE_separate4$O_E_t)),
    tau0=tau0,
    tau=tau
  )
  
  CUSUM_separate<-O_E_CUSUM(O_t=OE_separate$O_t,E_t=OE_separate$E_t,theta1=theta1,theta2=theta2,h1=h_cpp[1],h2=h_cpp[2],tau0=tau0,tau=tau,rho1=rho1,rho2=rho2)
  
  plot_separate <- O_E_CUSUM_plot(CUSUM_separate,title = "separate")
  res=plot_separate$result_plot
  ggsave(filename=paste(job_type,mu_type,"separate.pdf",sep="_"),width = 6, height=4.5)
  save(res,file = paste(job_type,mu_type,"separate.RData",sep="_"))
  
  
}


