################################################################################################
######### Bosch et al. 2021. Functional reorganization - herbivory processes data analyses #####

## (1) Study region - trends in environment and habitat
## (2) Trait spaces analyses
## (3) Taxonomic and Functional composition analyses
## (4) Generalized Linear Mixed Models
## (5) Supplementary tables
## (6) Turf Bite Models controlling for abundance (MaxN)
## (7) Kelp Bite Models controlling for abundance (MaxN)

## Clean working directory

rm(list=ls())

## Code developed by Matthew McLean & Nestor E. Bosch on January 2021
## Please report any bug to nbosch1989@gmail.com or Mat e-mail to be added

######################
## LIBRARY PACKAGES ##
######################

library(dplyr)
library(rlang)
library(forcats)
if(!require(FD)){install.packages("FD"); library(FD)}
if(!require(pals)){install.packages("pals"); library(pals)}
if(!require(emdbook)){install.packages("emdbook"); library(emdbook)}
if(!require(MASS)){install.packages("MASS"); library(MAWSS)}
if(!require(ggplot2)){install.packages("ggplot2"); library(ggplot2)}
if(!require(ggExtra)){install.packages("ggExtra"); library(ggExtra)}
if(!require(cowplot)){install.packages("cowplot"); library(cowplot)}
library(lubridate)
library(tidync)
library(doParallel)
library(rerddap)
library(dendextend)

# Set work directory----

work.dir=("C:/Users/22373243/Dropbox/Projects/Analysis/Analysis_Tropicalization_Functional_Changes") ### Desktop

#work.dir=("~/workspace/Body_Size_Nestor") ### Ecocloud platform

# Set sub directories----

data.dir=paste(work.dir,"Data",sep="/")
plots.dir=paste(work.dir,"Plots",sep="/")
model.out=paste(work.dir,"ModelOut",sep="/")

### Set plotting defaults ----

Theme1 <-
  theme( # use theme_get() to see available options
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.background = element_rect(fill="white"),
    legend.key = element_blank(), # switch off the rectangle around symbols in the legend
    legend.text = element_text(size=10),
    legend.title = element_text(size=14, face="bold"),
    legend.position = "right",
    legend.direction="vertical",
    text=element_text(size=14),
    strip.text.y = element_text(size = 14,angle = 0),
    axis.title.x=element_text(vjust=0.3, size=12),
    axis.title.y=element_text(vjust=0.6, angle=90, size=12),
    axis.text.x=element_text(size=10,angle = 0),
    axis.text.y=element_text(size=10),
    axis.line.x=element_line(colour="black", size=0.5,linetype='solid'),
    axis.line.y=element_line(colour="black", size=0.5,linetype='solid'),
    plot.title = element_text(size=14, face="bold"),
    strip.background = element_blank())

### Function to calculate standard errors ----

se <- function(x) sd(x) / sqrt(length(x))
se.min <- function(x) (mean(x)) - se(x)
se.max <- function(x) (mean(x)) + se(x)

## Function to find the mode

getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

### Distance to df function

dist.to.df<-function(list_dist) {
  
  # checking input is a list (even with a single element)
  if ( is.list(list_dist)==FALSE ) {
    stop("Error: input 'list_dist' should be a list (even when only one dist object is provided)")
  }
  
  
  # names and number of dist objects
  dist_nm<-names(list_dist)
  dist_nb<-length(dist_nm)
  
  # labels of the fist dist object
  dist1_labels<-labels(list_dist[[1]])
  
  ## checking inputs #####
  
  # checking list contains only dist objects
  if ( any(unlist(lapply(list_dist,class))!="dist") ) {
    stop("Error: input 'list_dist' should contain only 'dist' object
              Correct using 'as.dist()' if necessary")
  }
  
  # checking all dist objects have names
  if ( is.null(dist_nm) | any( nchar(dist_nm)==0 ) ) {
    stop("Error: some of dist objects in 'list_dist' do not all have a name.
       Name all dist objects within the list (e.g. with distance metric).")
  }
  
  # checking 1st dist objects has labels
  if ( is.null(dist1_labels) | any( nchar(dist1_labels)==0 ) ) {
    stop("Error: first dist object in 'list_dist' does not have labels.
       Provide row names as character strings for the sets*variables matrix
       before computing distance")
  }
  
  
  ## reference for pairs of sets is first object of list ####
  
  # applying dist_long
  df_dist<-dendextend::dist_long(list_dist[[1]])
  
  # reversing and renaming order of columns
  df_dist<-data.frame(df_dist$cols, df_dist$rows, df_dist$distance)
  names(df_dist)<-c("x1","x2",dist_nm[1])
  
  
  ## if other dist object, binding values with first ones as new column(s) ####
  if (dist_nb >1) {
    
    # loop on dist objects
    for (k in 2:length(list_dist) )
    {
      # name of dist objects
      k_nm<-dist_nm[k]
      
      # names of its labels
      k_labels<-labels(list_dist[[k_nm]])
      
      # checking same labels than first dist object
      if ( is.null(k_labels) |
           any( dist1_labels != k_labels ) ) {
        stop(paste0("Error:  element ",k_nm, " does not have the same labels
      than first dist object ", dist_nm[1], ".Check labels of inputs.") )
      }
      
      # applying dist_long
      df_k<-dendextend::dist_long(list_dist[[k_nm]])
      
      # merging with first dataframe as variable with dist object name
      df_dist[[k_nm]]<-df_k$distance
      
    }# end of k
    
    
  } # end of if at least 2 dist objects
  
  
  # output
  return(df_dist)
  
  
} # end of function

## Functional alpha diversity for a set of N assemblages from Chao et al. 2019 ----

FD_MLE <- function(data, dij, tau, q){
  dij <- as.matrix(dij)
  dij[which(dij>tau,arr.ind = T)] <- tau
  a <- as.vector((1 - dij/tau) %*% data )  
  data <- data[a!=0]
  a <- a[a!=0]
  v <- data/a
  nplus <- sum(data)
  if(q==1){
    exp(sum(-v*a/nplus*log(a/nplus)))
  }else{
    (sum(v*(a/nplus)^q))^(1 / (1-q))
  }
}

FunD <- function(data, dij, tau, q, boot, datatype){
  EstiBootComm.Func = function(data, distance, datatype){
    if (datatype=="abundance") {
      n = sum(data)
    } else if (datatype=="incidence_raw") {
      n <- data[1]
      data <- data[-1]
      u=sum(data)
    }
    distance = as.matrix(distance)
    dij = distance[data!=0, data!=0]
    
    X = data[data>0]
    f1 <- sum(X == 1) ; f2 <- sum(X == 2) 	
    f0.hat <- ceiling(ifelse(f2>0, ((n-1)/n)*f1^2/2/f2, ((n-1)/n)*f1*(f1-1)/2))
    if (datatype=="abundance") {
      C1 = ifelse(f2>0, 1-f1*(n-1)*f1/n/((n-1)*f1+2*f2), 1-f1*(n-1)*(f1-1)/n/((n-1)*(f1-1)+2))
      W <- (1 - C1)/sum(X/n*(1-X/n)^n) 
      Prob.hat.Unse <- rep((1-C1)/f0.hat, f0.hat)
    } else if (datatype=="incidence_raw") {
      C1 = ifelse(f2>0, 1-f1/u*(n-1)*f1/((n-1)*f1+2*f2), 1-f1*(n-1)*(f1-1)/u/((n-1)*(f1-1)+2))
      W <- (1 - C1)/sum(X/u*(1-X/n)^n) 
      Prob.hat.Unse <- rep(u/n*(1-C1)/f0.hat, f0.hat)	
    }
    
    Prob.hat <- X/n*(1-W*(1-X/n)^n)
    Prob <- c(Prob.hat, Prob.hat.Unse)
    
    F.1 <- sum(dij[, X==1]) ; F.2 <- sum(dij[, X==2])
    F11 <- sum(dij[X==1, X==1]) ; F22 <- sum(dij[X==2, X==2])
    #
    if (datatype=="abundance") {
      F.0hat <- ifelse(F.2 > 0, ((n-1)/n) * (F.1^2/(2 * F.2)), ((n-1)/n)*(F.1*(F.1-0.01)/(2)))
      F00hat <- ifelse(F22 > 0, ((n-2)* (n-3)* (F11^2)/(4* n* (n-1)* F22)), ((n-2)* (n-3)* (F11*(F11-0.01))/(4 *n * (n-1))) )
    } else if (datatype=="incidence_raw") {
      F.0hat <- ifelse(F.2 > 0, ((n-1)/n) * (F.1^2/(2 * F.2)), ((n-1)/n)*(F.1*(F.1-0.01)/(2)))
      F00hat <- ifelse(F22 > 0, ((n-1)^2 * (F11^2)/(4* n* n* F22)), ((n-1)* (n-1)* (F11*(F11-0.01))/(4 *n * n)) )
    }
    if (f0.hat==0) {
      d=dij
    } else if (f0.hat==1) {
      random_dij = as.vector(rmultinom(1, 1000, rep(1/(length(X)*f0.hat), length(X)*f0.hat) ) )/1000
      d.0bar <- matrix(random_dij*F.0hat, length(X), f0.hat, byrow = T)
      d00 = matrix(0, f0.hat, f0.hat)
      d <- cbind(dij, d.0bar )
      aa <- cbind(t(d.0bar), d00 )
      d <- rbind(d, aa)
      diag(d) = 0
    } else {
      random_dij = as.vector(rmultinom(1, 1000, rep(1/(length(X)*f0.hat), length(X)*f0.hat) ) )/1000
      d.0bar <- matrix(random_dij*F.0hat, length(X), f0.hat, byrow = T)
      
      fo.num = (f0.hat * (f0.hat-1) )/2
      random_d00 = as.vector(rmultinom(1, 1000, rep(1/fo.num, fo.num) ) )/1000
      d00 = matrix(0, f0.hat, f0.hat)
      d00[upper.tri(d00)] = (F00hat/2)*random_d00
      d00 <- pmax(d00, t(d00))###signmatrix
      d <- cbind(dij, d.0bar )
      aa <- cbind(t(d.0bar), d00 )
      d <- rbind(d, aa)
      diag(d) = 0
    }
    return(list("pi" = Prob,"dij" = d))
  }
  dij <-  as.matrix(dij)
  out <- as.vector(dij)
  out <- out[out!=0]
  dmin <- min(out)
  dmax <- max(out)
  if (datatype=="incidence_raw") {
    data <- lapply(data, function(i) {
      c(ncol(i), rowSums(i))
    })
  }
  if (datatype=="abundance") {
    if(length(data)!=1){
      tmp <- apply(do.call(cbind,lapply(data, FUN = function(x) x/sum(x))), 1, mean)
      dmean <-  sum ( (tmp %*% t(tmp) ) * dij) 
    }else{
      tmp <- data[[1]]/sum(data[[1]])
      dmean <-  sum ( (tmp %*% t(tmp) ) * dij)   
    }
  } else {
    if(length(data)!=1){
      tmp <- apply(do.call(cbind,lapply(data, FUN = function(x) x[-1]/sum(x[-1]))), 1, mean)
      dmean <-  sum ( (tmp %*% t(tmp) ) * dij) 
    }else{
      tmp <- data[[1]][-1]/sum(data[[1]][-1])
      dmean <-  sum ( (tmp %*% t(tmp) ) * dij)   
    }
  }
  FD.CI = function(data, dij, tau, q, datatype){
    if (datatype == "abundance") {
      qFun = FD_MLE(data, dij, tau, q)
    } else {
      qFun = FD_MLE(data[-1], dij, tau, q)
    }
    if(boot!=0){
      BT = EstiBootComm.Func(data, dij, datatype)
      p_hat = BT[[1]]
      dij_boot = BT[[2]]
      dij_boot <-  as.matrix(dij_boot)
      dij_boot <- replace(dij_boot, dij_boot==0, 10^(-10))
      for (i in seq_len(nrow(dij_boot))) {
        dij_boot[i, i] <- 0
      }
      if (datatype=="abundance") {
        n=sum(data)
        Boot.X = rmultinom(boot, n, p_hat)
      } else {
        n=data[1]
        Boot.X = t(sapply(p_hat,function(i) rbinom(boot, n, i)))
      }
      qFun_sd = sd(sapply(seq_len(ncol(Boot.X)), function(i) {
        FD_MLE(Boot.X[, i], dij_boot, tau, q)
      }))
    }else{
      qFun_sd = 0
    }
    LCL = max(0, qFun - qnorm(0.975) * qFun_sd)
    UCL = qFun + qnorm(0.975) * qFun_sd
    a = round(c(qFun, qFun_sd, LCL, UCL), 4)
    a
  }
  
  Funq <- function(data, datatype){
    dminFDforq <- t(sapply(q, FUN = function(q) FD.CI(data, dij, dmin, q, datatype) ))
    dmaxFDforq <- t(sapply(q, FUN = function(q) FD.CI(data, dij, dmax, q, datatype) ))
    dmeanFDforq <-t(sapply(q, FUN = function(q) FD.CI(data, dij, dmean, q, datatype) ))
    out <- data.frame(rep(q,3), rbind(dminFDforq,dmaxFDforq,dmeanFDforq),rep(c("dmin","dmax","dmean"),each=length(q)))
  }
  Funtau <- function(data, datatype){
    q0FDfortau <- t(sapply(tau, FUN = function(tau) FD.CI(data, dij, tau, 0, datatype) ))
    q1FDfortau <- t(sapply(tau, FUN = function(tau) FD.CI(data, dij, tau, 1, datatype) ))
    q2FDfortau <- t(sapply(tau, FUN = function(tau) FD.CI(data, dij, tau, 2, datatype) ))
    out <- data.frame(rep(tau,3), rbind(q0FDfortau,q1FDfortau,q2FDfortau),rep(c("0","1","2"),each=length(tau)))
  }
  
  if(length(data)!=1){
    name = names(data)
    Outputforq <- data.frame(do.call(rbind,lapply(data, Funq, datatype=datatype)), rep(name, each=3*length(q)), row.names = NULL)
    Outputfortau <- data.frame(do.call(rbind,lapply(data, Funtau, datatype=datatype)), rep(name, each=3*length(tau)), row.names = NULL)
  }else{
    name = names(data)
    Outputforq <- data.frame(Funq(data[[1]], datatype), name, row.names = NULL)
    Outputfortau <- data.frame(Funtau(data[[1]], datatype), name, row.names = NULL)
  }
  colnames(Outputforq) <- c("q","estimate", "s.e.", "LCL", "UCL", "tau","site")
  colnames(Outputfortau) <- c("tau","estimate", "s.e.", "LCL", "UCL", "q","site")
  
  Output <- list(forq = Outputforq, fortau = Outputfortau)
  return(Output)
}

## Functional beta-diversity function from Chao et al. 2019 - Ecological Monograph ----

beta.fd.hill <- function(asb_sp_w,
                         sp_dist,
                         q=c(0,1,2),
                         tau="mean",
                         beta_type="Jaccard",
                         check.input = TRUE,
                         store.details = TRUE) {
  
  #  distance between species stored in a matrix  ####
  sp_sp_dist<-sp_dist
  
  if (is.matrix(sp_sp_dist)==FALSE) {
    sp_sp_dist<-as.matrix(sp_sp_dist)
  }
  
  
  ## check inputs if required #####
  if (check.input == TRUE) {
    
    if (any(is.na(sp_dist))) {
      stop("Error: The species*coordinates matrix contains NA. Please check.")
    }
    if (is.null(rownames(sp_sp_dist))) {
      stop("Error: No row names provided in species*coordinates matrix.
             Please add species names as row names.")
    }
    if (any(is.na(asb_sp_w))) {
      stop("Error: The species*weights matrix contains NA. Please check.")
    }
    if (is.null(rownames(asb_sp_w))) {
      stop("Error: No row names provided in species*weights dataframe.
             Please add assemblages names as row names.")
    }
    if (any(! (colnames(asb_sp_w) %in% rownames(sp_sp_dist) ) ) ) {
      stop(paste("Error: Mismatch between names in species*weight and
                   species*coordinates matrix. Please check."))
    }
    
    if(any(! q %in% c(0,1,2) ) ) {
      stop(paste("Error: q should be 0, 1 and/or 2.
                  Please check."))
    }
    
    if(any(! tau %in% c("min", "mean", "max") ) ) {
      stop(paste("Error: tau should be 'mean' or 'max'. Please check."))
    }
    
    if(any(! beta_type %in% c("Jaccard", "Sorensen") ) ) {
      stop(paste("Error: beta_type should be 'Jaccard' or 'Sorensen'. Please check."))
    }
    
  }# end of checking inputs
  
  #  preliminary operations ####
  
  # ensuring species are in the same order in both matrices
  sp_sp_dist<-sp_sp_dist[colnames(asb_sp_w),colnames(asb_sp_w)]
  
  # names and number of assemblages
  asb_nm<-row.names(asb_sp_w)
  asb_nb<-length(asb_nm)
  
  # computing total weight per assemblage  ----
  asb_totw<-apply(asb_sp_w,1,sum)
  if(any(asb_totw ==0) ) {
    stop(paste("Error: all assemblages should contain at least one species.
               Please check."))
  }
  
  
  # computing tau as mean or max on distances ----
  tau_dist<-NULL
  
  if( tau=="min") {
    tau_dist<-min(sp_dist)
    
    # special case of null distance outside diagonal
    if (tau_dist==0) {
      tau_dist<-min(sp_dist[sp_dist!=0])
      cat("Warning: some species has null functional distance,
          'tau' was set to the minimum non-null distance")
    }
  }
  
  if( tau=="mean") {
    tau_dist<-mean(sp_dist)
  }
  
  if( tau=="max") {
    tau_dist<-max(sp_dist)
  }
  
  # applying tau threshold to distance matrix
  dij_tau<-sp_sp_dist
  dij_tau[which(dij_tau>tau_dist,arr.ind = T)] <- tau_dist
  
  
  
  # dissimilarity between assemblages ####
  
  # list to store diversity values
  beta_fd_q<-list()
  malpha_fd_q<-list()
  gamma_fd_q<-list()
  
  # matrices to store diversity values of order q
  mat_res<-matrix(NA, asb_nb, asb_nb, dimnames = list(asb_nm, asb_nm) )
  if (0 %in% q)
  {
    beta_fd_q$q0<-mat_res
    gamma_fd_q$q0<-mat_res
    malpha_fd_q$q0<-mat_res
  }
  if (1 %in% q)
  {
    beta_fd_q$q1<-mat_res
    gamma_fd_q$q1<-mat_res
    malpha_fd_q$q1<-mat_res
  }
  if (2 %in% q)
  {
    beta_fd_q$q2<-mat_res
    gamma_fd_q$q2<-mat_res
    malpha_fd_q$q2<-mat_res
  }
  
  
  # combinations of assemblages
  asb_pairs<-t (combn(asb_nm,2))
  colnames(asb_pairs)<-paste0("asb.",1:2)
  asb_pairs_nb<-nrow(asb_pairs)
  
  # loop on pairs of assemblages
  for (x in 1:asb_pairs_nb )
  {
    
    # names of assemblages in the pair x
    asb_nm_x<-asb_pairs[x,]
    
    # computing core variables for the pair of assemblages ----
    # notations as in Chao et al 2019, page 16, bottom right (with p for +)
    
    # weights of species (rows) in the 2 assemblages (columns)
    # (nik, bottom right p16)
    x_nik<-t(asb_sp_w[asb_nm_x,])
    
    # total weight of species in the 2 assemblages
    x_npp<-sum(x_nik)
    
    # total weight of each species among the 2 assemblages
    x_nip<-apply(x_nik,1,sum)
    
    # keeping only weight and distance of species present in pair of assemblages
    x_sp<-names(which(x_nip>0))
    x_nip<-x_nip[x_sp]
    x_nik<-x_nik[x_sp,]
    x_sp_dist<-dij_tau[x_sp,x_sp]
    
    # weight of functionally distinct group of species (aik, ai+ and vi+)
    x_sp_aik<-(1-x_sp_dist/tau_dist) %*% x_nik
    x_sp_aip<-apply(x_sp_aik,1,sum)
    x_sp_vip<-x_nip/x_sp_aip
    
    # species occurrences
    x_sp_01<-x_sp_aik
    x_sp_01[which(x_sp_01>0)]<-1
    
    # computing alpha, gamma and beta diversity according to levels of q ----
    
    # q=0 ----
    if (0 %in% q)
    {
      # alpha diversity (eq 7a) with special case of 0^0=0
      # hence sum of species attribute contribution depends on their occurrence
      x_malpha_q0<-sum(x_sp_vip*x_sp_01)/2
      
      # gamma diversity (eq 6a)
      x_gamma_q0<-sum(x_sp_vip)
      
      # beta Jaccard or Sorensen
      if (beta_type=="Sorensen") {
        x_beta_q0<-(x_gamma_q0/x_malpha_q0)-1
      }
      if (beta_type=="Jaccard") {
        x_beta_q0<-(1-(x_malpha_q0/x_gamma_q0))/(0.5)
      }
      
      # storing values
      malpha_fd_q$q0[asb_nm_x[2],asb_nm_x[1]]<-x_malpha_q0
      gamma_fd_q$q0[asb_nm_x[2],asb_nm_x[1]]<-x_gamma_q0
      beta_fd_q$q0[asb_nm_x[2],asb_nm_x[1]]<-x_beta_q0
      
    } # end of q=0
    
    # q=1 -----
    if (1 %in% q)
    {
      # alpha diversity (eq 7b) with special case of 0^0=0
      # hence sum of species attribute contribution depends on their occurrence
      x_malpha_q1<-0.5 * exp( (-1)* sum(x_sp_vip*(x_sp_aik/x_npp)*
                                          log(x_sp_aik/x_npp),na.rm=T ) )
      
      # gamma diversity (eq 6b)
      x_gamma_q1<- exp( (-1)* sum(x_sp_vip*(x_sp_aip/x_npp)*
                                    log(x_sp_aip/x_npp) ) )
      
      # beta Jaccard or Sorensen are identical
      x_beta_q1<-log(x_gamma_q1/x_malpha_q1)/log(2)
      
      # storing values
      malpha_fd_q$q1[asb_nm_x[2],asb_nm_x[1]]<-x_malpha_q1
      gamma_fd_q$q1[asb_nm_x[2],asb_nm_x[1]]<-x_gamma_q1
      beta_fd_q$q1[asb_nm_x[2],asb_nm_x[1]]<-x_beta_q1
      
    } # end of q=1
    
    
    # q=2 ----
    if (2 %in% q)
    {
      # alpha diversity (eq 7a) with special case of 0^0=0
      # hence sum of species attribute contribution depends on their occurrence
      x_malpha_q2<- 0.5 / ( sum(x_sp_vip*((x_sp_aik/x_npp)^2)) )
      
      # gamma diversity (eq 6a)
      x_gamma_q2<-1 / ( sum(x_sp_vip*((x_sp_aip/x_npp)^2)) )
      
      # beta Jaccard or Sorensen
      if (beta_type=="Sorensen") {
        x_beta_q2<-(1-(x_malpha_q2/x_gamma_q2))/(0.5)
      }
      if (beta_type=="Jaccard") {
        x_beta_q2<-(x_gamma_q2/x_malpha_q2)-1
      }
      
      # storing values
      malpha_fd_q$q2[asb_nm_x[2],asb_nm_x[1]]<-x_malpha_q2
      gamma_fd_q$q2[asb_nm_x[2],asb_nm_x[1]]<-x_gamma_q2
      beta_fd_q$q2[asb_nm_x[2],asb_nm_x[1]]<-x_beta_q2
    } # end of q=2
    
    
  } # end of loop on pairs
  
  
  # matrix with indices values as dist objects
  malpha_fd_q<-lapply(malpha_fd_q, as.dist)
  gamma_fd_q<-lapply(gamma_fd_q, as.dist)
  beta_fd_q<-lapply(beta_fd_q, as.dist)
  
  # returning outputs
  res<-beta_fd_q
  
  if (store.details==TRUE)
  {
    res<-list(beta_fd_q=beta_fd_q,
              details=list( malpha_fd_q=malpha_fd_q, gamma_fd_q=gamma_fd_q)
    )
  }
  return(res)
  
  
} # end of function

##################################################################################################################
############################## (1) Study region - Environmental Trends ###########################################

## Load libraries for creating map of WA with inset of Australia

library(rnaturalearth)
library(sp)
library(tidyverse)
library(ggspatial)
library(cowplot)

# TO install high resolution (so that Porto Santo and Desertas show up) we need the highres maps:

library(devtools)
install_github("https://github.com/ropensci/rnaturalearthhires")

world <- ne_countries(scale = "large", returnclass = "sf") # create a world map

## Plot Sout-Western Australia - Kalbarri to Esperance

# Add rectanlge around PG

d <- data.frame(x1=c(113.5), x2=c(115), y1=c(-29), y2=c(-28), r=c("PG"))

Map1<-ggplot(world) + geom_sf(fill="darkgrey") + coord_sf(xlim = c(113,123), ylim = c(-36,-27), expand = FALSE) + 
  theme_bw() +
  ggtitle(('(a)'))+
  annotation_scale(location="tl", style="ticks")+ # plot base map with the coordinates we need
  geom_rect(data=d, aes(xmin=x1,xmax=x2,ymin=y1,ymax=y2), fill="transparent", color="red")
Map1

# Create a plot for Australia

d.2 <- data.frame(x1=c(113), x2=c(123), y1=c(-36), y2=c(-27), r=c("SWA"))

A <- ggplot(world) + geom_sf(fill="white",colour="black") + coord_sf(xlim = c(113,160), ylim = c(-45,-9), expand = FALSE) + 
  geom_rect(data=d.2, aes(xmin=x1,xmax=x2,ymin=y1,ymax=y2), fill="transparent", color="red") +
  theme_void()
A

## Join Maps with inset Australia
library(cowplot)

gg_inset_map1 = ggdraw() +
  draw_plot(A,x=0.5,width = 0.5) +
  draw_plot(Map1, x = 0, width = 0.5,y=0.25,height = 0.5)
gg_inset_map1

## Export as png and pdf

setwd(plots.dir)
dir()

name<-'Study_Region_Broad'

ggsave(paste(name,".png",sep="."),width = 14, height = 12.5,units = "cm",dpi=600)
ggsave(paste(name,".pdf",sep="."),width = 14, height = 12.5,units = "cm",dpi=600,useDingbats=FALSE)


### Download specifc data for PG from Reynolds OISST data ----

# The information for the NOAA OISST data
rerddap::info(datasetid = "ncdcOisst21Agg_LonPM180", url = "https://coastwatch.pfeg.noaa.gov/erddap/")

# This function downloads and prepares data based on user provided start and end dates

citation('rerddap')
citation('heatwaveR')

OISST_sub_dl <- function(time_df){
  OISST_dat <- griddap(x = "ncdcOisst21Agg_LonPM180", 
                       url = "https://coastwatch.pfeg.noaa.gov/erddap/", 
                       time = c(time_df$start, time_df$end), 
                       zlev = c(0, 0),
                       latitude = c(-29, -28),
                       longitude = c(114, 115),
                       fields = "sst")$data %>% 
    mutate(time = as.Date(stringr::str_remove(time, "T00:00:00Z"))) %>% 
    dplyr::rename(t = time, temp = sst) %>% 
    dplyr::select(lon, lat, t, temp) %>% 
    na.omit()
}

# Date download range by start and end dates per year
dl_years <- data.frame(date_index = 1:5,
                       start = as.Date(c("1982-01-01", "1990-01-01", 
                                         "1998-01-01", "2006-01-01", "2014-01-01")),
                       end = as.Date(c("1989-12-31", "1997-12-31", 
                                       "2005-12-31", "2013-12-31", "2019-12-31")))

# Download all of the data with one nested request
# The time this takes will vary greatly based on connection speed
system.time(
  OISST_data <- dl_years %>% 
    dplyr::group_by(date_index) %>% 
    dplyr::group_modify(~OISST_sub_dl(.x)) %>% 
    dplyr::ungroup() 
) # 636 seconds, ~127 seconds per batch

OISST_data<-OISST_data%>%
  dplyr::select(lon, lat, t, temp)%>%
  glimpse()

## Group OISST temperature data

OISST_data<-OISST_data%>%
  dplyr::group_by(t)%>%
  summarise(temp=mean(temp))%>%
  glimpse()


## Export OISST data 

setwd(data.dir)
dir()

write.csv(OISST_data,"PG_Long_Term_SST.csv")

#### Bring in already downloaded OISST data

setwd(data.dir)
dir()

OISST_data<-read.csv('PG_Long_Term_SST.csv')%>%
  dplyr::select(-X)%>%
  glimpse()
OISST_data$t<-as.Date(OISST_data$t)

#### Calculate MHWs for PG ----

library(heatwaveR)

# Detect the events in a time series

ts <- ts2clm(OISST_data, climatologyPeriod = c("1982-01-01", "2019-12-31"))
mhw <- detect_event(ts)
MHW_cat <- category(mhw, S = TRUE, name = "WA") ## For categorizing heatwaves based on Hobay et al. 2018

# Look at the top few events

# View just a few metrics
table.MWHs<-as.data.frame(mhw$event %>% 
  dplyr::ungroup() %>%
  dplyr::select(event_no, duration, date_start, date_peak, intensity_max, intensity_cumulative) %>% 
  dplyr::arrange(-intensity_max)) %>% 
  head(5)

write.csv(table.MWHs,"table.MWHs")

tail(MHW_cat)

# View just a few metrics
mhw$event %>% 
  dplyr::ungroup() %>%
  dplyr::select(event_no, duration, date_start, date_peak, intensity_max, intensity_cumulative) %>% 
  dplyr::arrange(-intensity_max) %>% 
  head(5)

## Visualize MHWs for PG ----

## Plots via ggplot 2

## Subset for dates you want to plot

myfunc <- function(x,y){mhw$climatology[mhw$climatology$t >= x & mhw$climatology$t <= y,]}

DATE1 <- as.Date("2011-01-01")
DATE2 <- as.Date("2013-12-31")

Test <- myfunc(DATE1,DATE2)  

clim_cat <- Test%>%
  dplyr::mutate(diff = thresh - seas,
                thresh_2x = thresh + diff,
                thresh_3x = thresh_2x + diff,
                thresh_4x = thresh_3x + diff) %>%
  glimpse()


# Set line colours
lineColCat <- c(
  "Temperature" = "black",
  "Climatology" = "gray20",
  "Threshold" = "darkgreen",
  "2x Threshold" = "darkgreen",
  "3x Threshold" = "darkgreen",
  "4x Threshold" = "darkgreen"
)

# Set category fill colours
fillColCat <- c(
  "Moderate" = "#ffc866",
  "Strong" = "#ff6900",
  "Severe" = "#9e0000",
  "Extreme" = "#2d0000"
)

## Plot via ggplot2

## Create rectangles for stereo-DOVs and RUV samples

DATE1 <- as.Date("2013-05-01")
DATE2 <- as.Date("2013-07-01")

DATE3 <- as.Date("2013-11-01")
DATE4 <- as.Date("2013-12-01")

mhw.PG<-ggplot(data = clim_cat, aes(x = t, y = temp)) +
  geom_flame(aes(y2 = thresh, fill = "Moderate")) +
  geom_flame(aes(y2 = thresh_2x, fill = "Strong")) +
  geom_flame(aes(y2 = thresh_3x, fill = "Severe")) +
  geom_flame(aes(y2 = thresh_4x, fill = "Extreme")) +
  geom_line(aes(y = thresh_2x, col = "2x Threshold"), size = 0.3, linetype = "dashed",show.legend = F) +
  geom_line(aes(y = thresh_3x, col = "3x Threshold"), size = 0.3, linetype = "dotdash",show.legend = F) +
  geom_line(aes(y = thresh_4x, col = "4x Threshold"), size = 0.3, linetype = "dotted",show.legend = F) +
  geom_line(aes(y = seas, col = "Climatology"), size = 0.3,show.legend = F) +
  geom_line(aes(y = thresh, col = "Threshold"), size = 0.3,show.legend = F) +
  geom_line(aes(y = temp, col = "Temperature"), size = 0.3,show.legend = F) +
  scale_colour_manual(name = NULL, values = lineColCat,
                      breaks = c("Temperature", "Climatology", "Threshold",
                                 "2x Threshold", "3x Threshold", "4x Threshold")) +
  scale_fill_manual(name = NULL, values = fillColCat, guide = FALSE) +
  scale_x_date(date_labels = "%b %Y") +
  guides(colour = guide_legend(override.aes = list(linetype = c("solid", "solid", "solid",
                                                                "dashed", "dotdash", "dotted"),
                                                   size = c(0.6, 0.7, 0.7, 0.7, 0.7, 0.7)))) +
  labs(y = "Temperature (°C)", x = NULL)+
  # geom_rect(aes(xmin = DATE1, xmax = DATE2, 
  #                             ymin = -Inf, ymax = Inf), alpha = 0.01) + 
  # geom_rect(aes(xmin = DATE3, xmax = DATE4, 
  #               ymin = -Inf, ymax = Inf), alpha = 0.01) + 
  ggtitle('(b)')+
  theme_classic()+
  Theme1
mhw.PG

## Export as png and pdf

setwd(plots.dir)
dir()

name<-'MHWs_PG'

ggsave(paste(name,".png",sep="."),width = 21, height = 6,units = "cm",dpi=600)
ggsave(paste(name,".pdf",sep="."),width = 21, height = 6,units = "cm",dpi=600,useDingbats=FALSE)

## Calculate MCSs for WA ----

ts_MCS <- ts2clm(OISST_data, climatologyPeriod = c("1982-01-01", "2019-12-31"), pctile = 10)
MCS <- detect_event(ts_MCS, coldSpells = T)
MCS_cat <- category(MCS, S = TRUE, name = "WA")

# Look at the top few events

table.MCSs<-as.data.frame(MCS$event %>% 
                            dplyr::ungroup() %>%
                            dplyr::select(event_no, duration, date_start, date_peak, intensity_max, intensity_cumulative) %>% 
                            dplyr::arrange(intensity_max)) %>% 
  head(20)

write.csv(table.MCSs,"table.MCSs.csv")

tail(MCS_cat)

## Customize plots via ggplot2

## Subset for dates you want to plot

myfunc <- function(x,y){MCS$climatology[MCS$climatology$t >= x & MCS$climatology$t <= y,]}

DATE5 <- as.Date("2016-01-01")
DATE6 <- as.Date("2019-12-31")

Test <- myfunc(DATE5,DATE6)  

# Create category breaks and select slice of data.frame

MCS_clim_cat <-Test%>%
  dplyr::mutate(diff = thresh - seas,
                thresh_2x = thresh + diff,
                thresh_3x = thresh_2x + diff,
                thresh_4x = thresh_3x + diff) %>%
  glimpse()

# Set line colours
lineColCat <- c(
  "Temperature" = "black",
  "Climatology" = "gray20",
  "Threshold" = "darkgreen",
  "2x Threshold" = "darkgreen",
  "3x Threshold" = "darkgreen",
  "4x Threshold" = "darkgreen"
)

# Set category fill colours
fillColCat <- c(
  "Moderate" = "#A4D4E0",
  "Strong" = "#5B80A6",
  "Severe" = "#2A3C66",
  "Extreme" = "#111433"
)

## Plot in ggplot2

DATE7 <- as.Date("2019-05-01")
DATE8 <- as.Date("2019-07-01")

DATE9 <- as.Date("2019-11-01")
DATE10 <- as.Date("2019-12-01")

mcs.PG<-ggplot(data = MCS_clim_cat, aes(x = t, y = temp)) +
  geom_flame(aes(y = thresh, y2 = temp, fill = "Moderate")) +
  geom_flame(aes(y = thresh_2x, y2 = temp, fill = "Strong")) +
  geom_flame(aes(y = thresh_3x, y2 = temp, fill = "Severe")) +
  geom_flame(aes(y = thresh_4x, y2 = temp, fill = "Extreme")) +
  geom_line(aes(y = thresh_2x, col = "2x Threshold"), size = 0.3, linetype = "dashed",show.legend = F) +
  geom_line(aes(y = thresh_3x, col = "3x Threshold"), size = 0.3, linetype = "dotdash",show.legend = F) +
  geom_line(aes(y = thresh_4x, col = "4x Threshold"), size = 0.3, linetype = "dotted",show.legend = F) +
  geom_line(aes(y = seas, col = "Climatology"), size = 0.3,show.legend = F) +
  geom_line(aes(y = thresh, col = "Threshold"), size = 0.3,show.legend = F) +
  geom_line(aes(y = temp, col = "Temperature"), size = 0.3,show.legend = F) +
  scale_colour_manual(name = NULL, values = lineColCat,
                      breaks = c("Temperature", "Climatology", "Threshold",
                                 "2x Threshold", "3x Threshold", "4x Threshold")) +
  scale_fill_manual(name = NULL, values = fillColCat, guide = FALSE) +
  scale_x_date(date_labels = "%b %Y") +
  guides(colour = guide_legend(override.aes = list(linetype = c("solid", "solid", "solid",
                                                                "dashed", "dotdash", "dotted"),
                                                   size = c(0.6, 0.7, 0.7, 0.7, 0.7, 0.7)))) +
  labs(y = "Temperature (°C)", x = NULL)+
  # geom_rect(aes(xmin = DATE7, xmax = DATE8, 
  #               ymin = -Inf, ymax = Inf), alpha = 0.01) + 
  # geom_rect(aes(xmin = DATE9, xmax = DATE10, 
  #               ymin = -Inf, ymax = Inf), alpha = 0.01) + 
  ggtitle('(c)')+
  theme_classic()+
  Theme1
mcs.PG

## Export as png and pdf

setwd(plots.dir)
dir()

name<-'MCS_PG'

ggsave(paste(name,".png",sep="."),width = 21, height = 6,units = "cm",dpi=600)
ggsave(paste(name,".pdf",sep="."),width = 21, height = 6,units = "cm",dpi=600,useDingbats=FALSE)

## Create Palletes and export

# The MCS colour palette
MCS_colours <- c(
  "Moderate" = "#A4D4E0",
  "Strong" = "#5B80A6",
  "Severe" = "#2A3C66",
  "Extreme" = "#111433"
)

# The MHW colour palette
MHW_colours <- c(
  "Moderate" = "#ffc866",
  "Strong" = "#ff6900",
  "Severe" = "#9e0000",
  "Extreme" = "#2d0000"
)

# Create the colour palette for plotting by itself
colour_palette <- data.frame(category = factor(c("I Moderate", "II Strong", "III Severe", "IV Extreme"),
                                               levels = c("I Moderate", "II Strong", "III Severe", "IV Extreme")),
                             MHW = c(MHW_colours[1], MHW_colours[2], MHW_colours[3], MHW_colours[4]),
                             MCS = c(MCS_colours[1], MCS_colours[2], MCS_colours[3], MCS_colours[4])) %>% 
  pivot_longer(cols = c(MHW, MCS), names_to = "event", values_to = "colour")

# Show the palettes side-by-side
ggplot(data = colour_palette, aes(x = category, y = event)) +
  geom_tile(fill = colour_palette$colour) +
  coord_cartesian(expand = F) +
  labs(x = NULL, y = NULL)

## Export as png and pdf

setwd(plots.dir)
dir()

name<-'Palette'

ggsave(paste(name,".png",sep="."),width = 16, height = 2,units = "cm",dpi=600)
ggsave(paste(name,".pdf",sep="."),width = 16, height = 2,units = "cm",dpi=600,useDingbats=FALSE)

#### Long-term trends in habitat cover ----

setwd(data.dir)
dir()

dat<-read.csv('PG_Long_Term_Habitat.csv')%>%
  filter(habitat.type%in%c('Turf','Canopy_OT','Ecklonia'))%>%
  dplyr::select(-X)%>%
  filter(site%in%c('PGN1','PGN2','PGS2','PGS3'))%>%
  glimpse()
unique(dat$site)
n_distinct(dat$id)

## Get summaries of transects per year

test<-dat%>%
  group_by(year,site)%>%
  summarise(N=n_distinct(id))%>%
  glimpse()

## Check distribution of variables - Turf, Canopy_OT, Ecklonia

n_distinct(dat$habitat.type) ## 8 groups
unique(dat$habitat.type)

## Summarize data per year - for Turf, Canopy_OT and Kelp

dat.sum<-dat%>%
  filter(habitat.type%in%c('Turf','Canopy_OT','Ecklonia'))%>%
  group_by(year,habitat.type)%>%
  summarise(cover_mean=mean(perc.cover),
            cover_se=se(perc.cover))%>%
  glimpse()

## Relevel levels of the factor

unique(dat$habitat.type)

dat<-dat%>%
  mutate(habitat.type=fct_relevel(habitat.type,"Turf","Canopy_OT","Ecklonia"))%>%
  glimpse()
dat$year<-as.factor(dat$year)

dat.sum<-dat.sum%>%
  mutate(habitat.type=fct_relevel(habitat.type,"Turf","Canopy_OT","Ecklonia"))%>%
  glimpse()
dat.sum$year<-as.factor(dat.sum$year)

## Plot via ggplot

ggplot(data=dat.sum,aes(x=year,y=cover_mean))+
  ggtitle('(e)')+
  ylab('Reef cover (%)')+
  xlab('')+
  geom_jitter(data=dat,aes(x=year,y=perc.cover,colour=habitat.type),position=position_jitterdodge(jitter.width = 0.2, jitter.height = 0,
                                                                      dodge.width = 0.75, seed = NA), size=2,alpha=0.5,show.legend = F)+
  geom_errorbar(data = dat.sum, aes(x = year,y=cover_mean,
                                    ymin=cover_mean-cover_se, ymax=cover_mean+cover_se,fill=habitat.type),color = "black",width=.2,
                position=position_jitterdodge(jitter.width = 0.2, jitter.height = 0,dodge.width = 0.75, seed = NA), size=1)+
  geom_point(data=dat.sum,
             aes(x=year,y=cover_mean,fill=habitat.type),color="black",
             size=4.5,alpha=2,
             position=position_dodge(.75),show.legend = F)+
  geom_point(data = dat.sum, aes(x=year,y=cover_mean,colour=habitat.type),
             size=4,position=position_dodge(.75),show.legend = T)+
  geom_path(data=dat.sum,
            aes(x=year,y=cover_mean,group = habitat.type,color=habitat.type),
            position=position_dodge(.75),
            size=1,show.legend = F)+
  scale_fill_manual(labels = c("Turf","Canopy_OT","Ecklonia"),
                      values=c("brown","darkgreen","blue"))+
  theme_classic()+
  Theme1+
  theme(legend.position = "bottom",
        legend.direction = "horizontal")

## Export as png and pdf

setwd(plots.dir)
dir()

name<-'Habitat_cover_no_legend'

ggsave(paste(name,".png",sep="."),width = 15, height = 10,units = "cm",dpi=600)
ggsave(paste(name,".pdf",sep="."),width = 15, height = 10,units = "cm",dpi=600,useDingbats=FALSE)

################################################################################################
############################## (2) Trait spaces analyses #######################################

## (A) Plot of species position in trait space ----

## Load trait data and retain only species present in the data

setwd(data.dir)
dir()

traits<-read.csv('traits_herbivores.csv')%>%
  filter(!is.na(ThermalMP_5_95))%>%
  dplyr::select(scientific,herbivore.guild,ThermalMP_5_95,MaxLength)%>%
  filter(!is.na(herbivore.guild))%>% ## 27 nominal herbivores
  filter(!grepl('sp',scientific))%>%
  glimpse()

#######################################
## PLOT SPECIES BY NAME AND POSITION ##
#######################################

par(mfrow=c(1,1))
plot(traits$ThermalMP_5_95, traits$MaxLength,
     xlab="Temperature Preference", ylab="MaxLength",
     cex=0)
text(traits$ThermalMP_5_95, traits$MaxLength,traits$scientific,
     cex=0.6)
title("Species Names")

####################################
## COLOR SPECIES BY TROPHIC GUILD ##
## Convex hull on trophic guilds ##
####################################

## First general plot to get the legend

ggplot()+
  geom_point(data=traits,
             aes(x=ThermalMP_5_95,y=MaxLength,fill=herbivore.guild),color="black",
             size=2.5,alpha=2)+
  geom_point(data = traits, aes(x=ThermalMP_5_95,y=MaxLength,colour=herbivore.guild),
             size=2)+
  geom_text(data=traits,aes(x=ThermalMP_5_95,y=MaxLength,label=scientific),vjust=-1)+
  geom_vline(xintercept = 23,linetype="dashed")+
  scale_x_continuous(limits = c(17.5,30),breaks = c(20.0,22.5,25,27.5))+
  theme_classic()+
  ggtitle('(a)')+
  xlab('Thermal affinity (°C)')+
  ylab('Max Length (cm)')+
  Theme1+
  theme(legend.position = "bottom",
        legend.direction = "horizontal")

## Export plot

setwd(plots.dir)
dir()

name<-'Trait_space_legend'

ggsave(paste(name,".png",sep="."),width = 21, height = 8,units = "cm",dpi=600)
ggsave(paste(name,".pdf",sep="."),width = 21, height = 8,units = "cm",dpi=600,useDingbats=FALSE)

## (B) stereo-DOV - species densities analyses ----

setwd(data.dir)
dir()

traits<-read.csv('traits_herbivores.csv')%>%
  filter(!is.na(ThermalMP_5_95))%>%
  #dplyr::select(scientific,herbivore.guild,ThermalMP_5_95,MaxLength)%>%
  filter(!is.na(herbivore.guild))%>% ## 27 nominal herbivores
  glimpse()

## Load density data

x<-as.vector(unique(traits$scientific))

setwd(data.dir)
dir()

dat<-read.csv('dovs_clean.csv')%>%
  filter(!scientific%in%c('Ophisternon bengalense'))%>%
  filter(!Class%in%c('Elasmobranchii'))%>%
  filter(!scientific%in%c('Cirripectes hutchinsi'))%>%
  #filter(year!=2017)%>% Multi-year data
  filter(!Site%in%c('PGS1'))%>%
  filter(scientific%in%x)%>%
  glimpse()
unique(dat$Site)
unique(dat$year)
n_distinct(dat$scientific) ## Of the 27 - only 12 species observed in s-DOVs - but note some species were pooled for analyses
unique(dat$scientific)

## Filter out species not present from trait data

x<-as.vector(unique(dat$scientific))

traits<-traits%>%
  filter(scientific%in%x)%>%
  glimpse()

## Change nomenclature to match matrix

traits$scientific<-gsub(' ','.',traits$scientific)

### Check for missing trait information

sum(is.na(traits))/prod(dim(traits))*100
apply(traits,2,function(col)sum(is.na(col))/length(col))*100

## Convert first column to row names

traits<-data.frame(traits)
rownames(traits) <- traits[,1] #Assigning row names from 1st column 
traits[,1] <- NULL #Removing the first column

## See sampling effort per year and site

test<-dat%>%
  dplyr::group_by(year,Site)%>%
  dplyr::summarise(N=n_distinct(id))%>%
  glimpse()

## Create a grouping factor for each time period

dat<-dat%>%
  mutate(time.period=ifelse(year%in%c('2006'),"pre-heatwave",
                            ifelse(year%in%c('2013'),"post-heatwave","post-cooling")))%>%
  glimpse()

unique(dat$time.period)

dat<-dat%>%
  mutate(time.period=fct_relevel(time.period,"pre-heatwave","post-heatwave","post-cooling"))%>%
  glimpse()

## Resample data to select 8 transects per time period

## Firt we convert to long format for each year separately 

library(tidyr)

dataframe=data.frame()

number<-1:99
time.period<-as.vector(unique(dat$time.period))
site<-as.vector(unique(dat$Site))

for(j in number) {

for (i in time.period) {
  
  for(k in site) {

Abun<-dat%>%
  filter(time.period%in%i)%>%
  filter(Site%in%k)%>%
  dplyr::group_by(id,scientific)%>%
  summarise(number=sum(number))%>%
  ungroup()%>%
  pivot_wider(names_from = "scientific",values_from = number) %>%
  mutate_all(~replace_na(., 0))%>%
  sample_n(8,replace = T)%>%
  glimpse()

## Reconvert to long format

Abun<-Abun%>%
  pivot_longer(cols = c(2:10), names_to = "scientific", values_to = "number")%>%
  glimpse()

## Create metadata 

metadata<-dat%>%
  dplyr::select(Location,Site,time.period,id)%>%
  glimpse()

new.dat<-left_join(Abun,metadata,by="id")%>%
  distinct(id,scientific,.keep_all = TRUE)%>%
  glimpse()

dataframe<-rbind(dataframe,new.dat)

  }
}
}
  
## Group data and obtain mean 

test<-dataframe%>%
  group_by(time.period,Site,scientific)%>%
  summarise(number=mean(number))%>%
  glimpse()

## Recheck sampling effort

dat<-test

## Clean working directory

rm(Abun,metadata,new.dat,dataframe)

## Summarize data by year

test<-dat%>%
  dplyr::group_by(time.period,scientific)%>%
  dplyr::summarise(number=sum(number))%>%
  glimpse()

## Check the shape of abundance distributions after transformation of species abundances

plot.pre.heatwave<-ggplot(test%>%filter(time.period%in%c('pre-heatwave'))%>%filter(number>0),aes(x=reorder(scientific,-number),y=number,fill=scientific))+
  geom_bar(stat = "identity",colour="black",show.legend = F)+
  scale_y_continuous(expand = c(0,0),limits = c(0,45))+
  theme_classic()+
  theme(axis.text.x = element_text(angle=90))

plot.post.heatwave<-ggplot(test%>%filter(time.period%in%c('post-heatwave'))%>%filter(number>0),aes(x=reorder(scientific,-number),y=number,fill=scientific))+
  geom_bar(stat = "identity",colour="black",show.legend = F)+
  scale_y_continuous(expand = c(0,0),limits = c(0,45))+
  theme_classic()+
  theme(axis.text.x = element_text(angle=90))

plot.post.cooling<-ggplot(test%>%filter(time.period%in%c('post-cooling'))%>%filter(number>0),aes(x=reorder(scientific,-number),y=number,fill=scientific))+
  geom_bar(stat = "identity",colour="black",show.legend = F)+
  scale_y_continuous(expand = c(0,0),limits = c(0,45))+
  theme_classic()+
  theme(axis.text.x = element_text(angle=90))

## Arrange plots in a grob

library(ggpubr)

ggarrange(plot.pre.heatwave,plot.post.heatwave,plot.post.cooling,nrow=1,ncol=3,align = 'hv')

## Export plot

setwd(plots.dir)
dir()

name<-'Supplementary_Dominance'

ggsave(paste(name,".png",sep="."),width = 21, height = 8,units = "cm",dpi=600)
ggsave(paste(name,".pdf",sep="."),width = 21, height = 8,units = "cm",dpi=600,useDingbats=FALSE)

## Check normal histograms

plot.new()
par(mfrow=c(1,2))

hist(test$number)
hist(log10(test$number+1))

## Create a species x abundances matrix for each year

library(tidyr)
library(dplyr)

abundances<-test%>%
  dplyr::select(time.period,scientific,number)%>%
  pivot_wider(names_from = "scientific",values_from = number) %>%
  mutate_all(~replace_na(., 0))%>%
  glimpse()

## Convert species to column names

abundances<-data.frame(abundances)
rownames(abundances) <- abundances[,1] #Assigning row names from 1st column 
abundances[,1] <- NULL #Removing the first column

## Transform species abundances

abundances <- log10(abundances+1)

#Order species in the community matrix in same order as traits

x<-as.vector(row.names(traits))

abundances <- abundances[,x]

## Check identical species are the same in trait and species matrix

identical(rownames(traits),colnames(abundances))

#######################################
## PLOT SPECIES DENSITY BEFORE/AFTER ##
#######################################

# REPEAT EACH SPECIES ROW AS MANY TIMES AS ITS ABUNDANCE (ON LOG SCALE - MULTIPLIED BY 10 AND ROUNDED)
# TO INCORPORATE ABUNDANCE INTO DENSITY ESTIMATE

n.times <- round(as.numeric(abundances[1,])*10,0)
density_2006 <- traits[rep(seq_len(nrow(traits)), n.times),]

n.times <- round(as.numeric(abundances[2,])*10,0)
density_2013 <- traits[rep(seq_len(nrow(traits)), n.times),]

n.times <- round(as.numeric(abundances[3,])*10,0)
density_2019 <- traits[rep(seq_len(nrow(traits)), n.times),]

## Calculate the community centroid (weighthed mean trait values)

## 2006

names(traits)

trait_space <- traits[,c("ThermalMP_5_95","mode.2006")]

centroids.2006 <- functcomp(trait_space, as.matrix(abundances[1,]))

## 2013

names(traits)

trait_space <- traits[,c("ThermalMP_5_95","mode.2013")]

centroids.2013 <- functcomp(trait_space, as.matrix(abundances[2,]))

## 2019

names(traits)

trait_space <- traits[,c("ThermalMP_5_95","mode.2019")]

centroids.2019 <- functcomp(trait_space, as.matrix(abundances[3,]))

## MARGINAL DENSITY USING COWPLOT ##

## 2006 

range(traits$ThermalMP_5_95)
range(traits$MaxLength)

## Relabel factor variables

unique(density_2006$herbivore.guild)

density_2006<-density_2006%>%
  mutate(herbivore.guild=fct_relevel(herbivore.guild,"Browser","Algal farmer","Scraper"))%>%
  glimpse()

## Estimate 50%, 75% AND 95% Highest Density Estimates 

getLevel <- function(x,y,prob=c(0.5,0.75,0.95)) {
  kk <- MASS::kde2d(x,y)
  dx <- diff(kk$x[1:2])
  dy <- diff(kk$y[1:2])
  sz <- sort(kk$z)
  c1 <- cumsum(sz) * dx * dy
  approx(c1, sz, xout = 1 - prob)$y
}
L95 <- getLevel(density_2006$ThermalMP_5_95,density_2006$mode.2006)

## Plot 2D density estimates with marginal densities for each trait axis
names(centroids.2006)
d.2006 <- ggplot(density_2006, aes(ThermalMP_5_95,mode.2006)) + 
  ggtitle('(a)')+
  ylab('Body size (cm)')+
  xlab('Thermal affinity (°C)')+
  stat_density_2d(geom = "polygon", aes(alpha = ..level..),show.legend = F,colour="black",breaks=L95) +
  geom_point(data=density_2006, aes(ThermalMP_5_95,mode.2006,
                                    fill=herbivore.guild),colour="black",size=1.5,show.legend = F) +
  geom_point(data=density_2006, aes(ThermalMP_5_95,mode.2006,
                                    color=herbivore.guild),size=1,show.legend = F) +
  geom_point(data=centroids.2006,aes(x=ThermalMP_5_95,y=mode.2006),shape=23,size=3,fill="orange")+
  scale_colour_manual(labels = c("Browser","Algal farmer","Scraper"),values=c("yellow4","tomato1","steelblue1"))+
  xlim(18.0405,30) + ylim(0,75)+
  theme_classic()+
  Theme1

xdens <- axis_canvas(d.2006, axis="x") +
  geom_density(data=density_2006, aes(x=ThermalMP_5_95), fill="black",alpha=0.2)

ydens <- axis_canvas(d.2006, axis="y", coord_flip = TRUE) +
  geom_density(data=density_2006, aes(x=mode.2006), fill="black",alpha=0.2) +
  coord_flip()

p1 <- insert_xaxis_grob(d.2006, xdens,position = "top")
p2 <- insert_yaxis_grob(p1, ydens,position = "right")
ggdraw(p2)

## 2013

range(traits$ThermalMP_5_95)
range(traits$MaxLength)

## Relabel factor variables

unique(density_2013$herbivore.guild)

density_2013<-density_2013%>%
  mutate(herbivore.guild=fct_relevel(herbivore.guild,"Sediment sucker","Browser","Algal farmer","Scraper"))%>%
  glimpse()

## Estimate 50%, 75% AND 95% Highest Density Estimates 

getLevel <- function(x,y,prob=c(0.5,0.75,0.95)) {
  kk <- MASS::kde2d(x,y)
  dx <- diff(kk$x[1:2])
  dy <- diff(kk$y[1:2])
  sz <- sort(kk$z)
  c1 <- cumsum(sz) * dx * dy
  approx(c1, sz, xout = 1 - prob)$y
}
L95 <- getLevel(density_2013$ThermalMP_5_95,density_2013$mode.2013)

## Plot 2D density estimates with marginal densities for each trait axis

d.2013 <- ggplot(density_2013, aes(ThermalMP_5_95,mode.2013)) + 
  ggtitle('(b)')+
  ylab('')+
  xlab('')+
  stat_density_2d(geom = "polygon", aes(alpha = ..level..),show.legend = F,colour="black",breaks=L95) +
  geom_point(data=density_2013, aes(ThermalMP_5_95,mode.2013,
                                    fill=herbivore.guild),colour="black",size=1.5,show.legend = F) +
  geom_point(data=density_2013, aes(ThermalMP_5_95,mode.2013,
                                    color=herbivore.guild),size=1,show.legend = F) +
  geom_point(data=centroids.2013,aes(x=ThermalMP_5_95,y=mode.2013),shape=23,size=3,fill="orange")+
  scale_colour_manual(labels = c("Sediment sucker","Browser","Algal farmer","Scraper"),
                      values=c("orchid2","yellow4","tomato1","steelblue1"))+
  xlim(18.0405,30) + ylim(0,75)+
  theme_classic()+
  Theme1

xdens <- axis_canvas(d.2013, axis="x") +
  geom_density(data=density_2013, aes(x=ThermalMP_5_95), fill="black",alpha=0.2)

ydens <- axis_canvas(d.2013, axis="y", coord_flip = TRUE) +
  geom_density(data=density_2013, aes(x=mode.2013), fill="black",alpha=0.2) +
  coord_flip()

p3 <- insert_xaxis_grob(d.2013, xdens,position = "top")
p4 <- insert_yaxis_grob(p3, ydens,position = "right")
ggdraw(p4)

## 2019

range(traits$ThermalMP_5_95)
range(traits$MaxLength)

## Relabel factor variables

unique(density_2019$herbivore.guild)

density_2019<-density_2019%>%
  mutate(herbivore.guild=fct_relevel(herbivore.guild,"Sediment sucker","Browser","Algal farmer","Scraper"))%>%
  glimpse()

## Estimate 50%, 75% AND 95% Highest Density Estimates 

getLevel <- function(x,y,prob=c(0.5,0.75,0.95)) {
  kk <- MASS::kde2d(x,y)
  dx <- diff(kk$x[1:2])
  dy <- diff(kk$y[1:2])
  sz <- sort(kk$z)
  c1 <- cumsum(sz) * dx * dy
  approx(c1, sz, xout = 1 - prob)$y
}
L95 <- getLevel(density_2019$ThermalMP_5_95,density_2019$mode.2019)

## Plot 2D density estimates with marginal densities for each trait axis

d.2019 <- ggplot(density_2019, aes(ThermalMP_5_95,mode.2019)) + 
  ggtitle('(c)')+
  ylab('')+
  xlab('')+
  stat_density_2d(geom = "polygon", aes(alpha = ..level..),show.legend = F,colour="black",breaks=L95) +
  geom_point(data=density_2019, aes(ThermalMP_5_95,mode.2019,
                                    fill=herbivore.guild),colour="black",size=1.5,show.legend = F) +
  geom_point(data=density_2019, aes(ThermalMP_5_95,mode.2019,
                                    color=herbivore.guild),size=1,show.legend = F) +
  geom_point(data=centroids.2019,aes(x=ThermalMP_5_95,y=mode.2019),shape=23,size=3,fill="orange")+
  scale_colour_manual(labels = c("Sediment sucker","Browser","Algal farmer","Scraper"),
                      values=c("orchid2","yellow4","tomato1","steelblue1"))+
  xlim(18.0405,30) + ylim(0,75)+
  theme_classic()+
  Theme1

xdens <- axis_canvas(d.2019, axis="x") +
  geom_density(data=density_2019, aes(x=ThermalMP_5_95), fill="black",alpha=0.2)

ydens <- axis_canvas(d.2019, axis="y", coord_flip = TRUE) +
  geom_density(data=density_2019, aes(x=mode.2019), fill="black",alpha=0.2) +
  coord_flip()

p5 <- insert_xaxis_grob(d.2019, xdens,position = "top")
p6 <- insert_yaxis_grob(p5, ydens,position = "right")

ggdraw(p6)

## Arrange plots in a grob

library(ggpubr)

ggarrange(p2,p4,p6,ncol=3,nrow = 1)

## Export plot

setwd(plots.dir)
dir()

name<-'2D-DOVs'

ggsave(paste(name,".png",sep="."),width = 21, height = 8,units = "cm",dpi=600)
ggsave(paste(name,".pdf",sep="."),width = 21, height = 8,units = "cm",dpi=600,useDingbats=FALSE)

####### Changes in the cetroid - abundance-weighted mean trait values

####### HDI STI ----

library(HDInterval)

#### 2006

## 50%

plot.new()
par(mfrow=c(1,1))

hist(density_2006$ThermalMP_5_95,freq = FALSE)

hdi<-hdi(density_2006$ThermalMP_5_95, credMass=0.5)
hdi
segments(hdi[1], 0, hdi[2], 0, lwd=3, col='red')

dens2 <- density(density_2006$ThermalMP_5_95)
lines(dens2, lwd=2, col='blue')

(hdiD2 <- hdi(dens2,credMass=0.5, allowSplit=TRUE))
segments(hdiD2[,1], 0, hdiD2[,2], 0, lwd=3, col='blue')

## 75%

plot.new()
par(mfrow=c(1,1))

hist(density_2006$ThermalMP_5_95,freq = FALSE)

hdi<-hdi(density_2006$ThermalMP_5_95, credMass=0.75)
hdi

segments(hdi[1], 0, hdi[2], 0, lwd=3, col='red')

dens2 <- density(density_2006$ThermalMP_5_95)
lines(dens2, lwd=2, col='blue')

(hdiD2 <- hdi(dens2,credMass=0.75, allowSplit=TRUE))
segments(hdiD2[, 1], 0, hdiD2[, 2], 0, lwd=3, col='blue')

## 95%

plot.new()
par(mfrow=c(1,1))

hist(density_2006$ThermalMP_5_95,freq = FALSE)

hdi<-hdi(density_2006$ThermalMP_5_95, credMass=0.95)
hdi

segments(hdi[1], 0, hdi[2], 0, lwd=3, col='red')

dens2 <- density(density_2006$ThermalMP_5_95)
lines(dens2, lwd=2, col='blue')

(hdiD2 <- hdi(dens2,credMass=0.95, allowSplit=TRUE))
segments(hdiD2[, 1], 0, hdiD2[, 2], 0, lwd=3, col='blue')

### 2013 STI

## 50%

plot.new()
par(mfrow=c(1,1))

hist(density_2013$ThermalMP_5_95,freq = FALSE)

hdi<-hdi(density_2013$ThermalMP_5_95, credMass=0.5)
hdi

segments(hdi[1], 0, hdi[2], 0, lwd=3, col='red')

dens2 <- density(density_2013$ThermalMP_5_95)
lines(dens2, lwd=2, col='blue')

(hdiD2 <- hdi(dens2,credMass=0.5, allowSplit=TRUE))
segments(hdiD2[,1], 0, hdiD2[,2], 0, lwd=3, col='blue')

## 75%

plot.new()
par(mfrow=c(1,1))

hist(density_2013$ThermalMP_5_95,freq = FALSE)

hdi<-hdi(density_2013$ThermalMP_5_95, credMass=0.75)
hdi

segments(hdi[1], 0, hdi[2], 0, lwd=3, col='red')

dens2 <- density(density_2013$ThermalMP_5_95)
lines(dens2, lwd=2, col='blue')

(hdiD2 <- hdi(dens2,credMass=0.75, allowSplit=TRUE))
segments(hdiD2[, 1], 0, hdiD2[, 2], 0, lwd=3, col='blue')

## 95%

plot.new()
par(mfrow=c(1,1))

hist(density_2013$ThermalMP_5_95,freq = FALSE)

hdi<-hdi(density_2013$ThermalMP_5_95, credMass=0.95)
hdi

segments(hdi[1], 0, hdi[2], 0, lwd=3, col='red')

dens2 <- density(density_2013$ThermalMP_5_95)
lines(dens2, lwd=2, col='blue')

(hdiD2 <- hdi(dens2,credMass=0.95, allowSplit=TRUE))
segments(hdiD2[, 1], 0, hdiD2[, 2], 0, lwd=3, col='blue')

### 2019 STI

## 50%

plot.new()
par(mfrow=c(1,1))

hist(density_2019$ThermalMP_5_95,freq = FALSE)

hdi<-hdi(density_2019$ThermalMP_5_95, credMass=0.5)
hdi

segments(hdi[1], 0, hdi[2], 0, lwd=3, col='red')

dens2 <- density(density_2019$ThermalMP_5_95)
lines(dens2, lwd=2, col='blue')

(hdiD2 <- hdi(dens2,credMass=0.5, allowSplit=TRUE))
segments(hdiD2[,1], 0, hdiD2[,2], 0, lwd=3, col='blue')

## 75%

plot.new()
par(mfrow=c(1,1))

hist(density_2019$ThermalMP_5_95,freq = FALSE)

hdi<-hdi(density_2019$ThermalMP_5_95, credMass=0.75)
hdi

segments(hdi[1], 0, hdi[2], 0, lwd=3, col='red')

dens2 <- density(density_2019$ThermalMP_5_95)
lines(dens2, lwd=2, col='blue')

(hdiD2 <- hdi(dens2,credMass=0.75, allowSplit=TRUE))
segments(hdiD2[, 1], 0, hdiD2[, 2], 0, lwd=3, col='blue')

## 95%

plot.new()
par(mfrow=c(1,1))

hist(density_2019$ThermalMP_5_95,freq = FALSE)

hdi<-hdi(density_2019$ThermalMP_5_95, credMass=0.95)
hdi

segments(hdi[1], 0, hdi[2], 0, lwd=3, col='red')

dens2 <- density(density_2019$ThermalMP_5_95)
lines(dens2, lwd=2, col='blue')

(hdiD2 <- hdi(dens2,credMass=0.95, allowSplit=TRUE))
segments(hdiD2[, 1], 0, hdiD2[, 2], 0, lwd=3, col='blue')

####### HDI MaxLength ----

library(HDInterval)

#### 2006

## 50%

plot.new()
par(mfrow=c(1,1))

hist(density_2006$mode.2006,freq = FALSE)

hdi<-hdi(density_2006$mode.2006, credMass=0.5)
hdi

segments(hdi[1], 0, hdi[2], 0, lwd=3, col='red')

dens2 <- density(density_2006$mode.2006)
lines(dens2, lwd=2, col='blue')

(hdiD2 <- hdi(dens2,credMass=0.5, allowSplit=TRUE))
segments(hdiD2[,1], 0, hdiD2[,2], 0, lwd=3, col='blue')

## 75%

plot.new()
par(mfrow=c(1,1))

hist(density_2006$mode.2006,freq = FALSE)

hdi<-hdi(density_2006$mode.2006, credMass=0.75)
hdi

segments(hdi[1], 0, hdi[2], 0, lwd=3, col='red')

dens2 <- density(density_2006$mode.2006)
lines(dens2, lwd=2, col='blue')

(hdiD2 <- hdi(dens2,credMass=0.75, allowSplit=TRUE))
segments(hdiD2[, 1], 0, hdiD2[, 2], 0, lwd=3, col='blue')

## 95%

plot.new()
par(mfrow=c(1,1))

hist(density_2006$mode.2006,freq = FALSE)

hdi<-hdi(density_2006$mode.2006, credMass=0.95)
hdi

segments(hdi[1], 0, hdi[2], 0, lwd=3, col='red')

dens2 <- density(density_2006$mode.2006)
lines(dens2, lwd=2, col='blue')

(hdiD2 <- hdi(dens2,credMass=0.95, allowSplit=TRUE))
segments(hdiD2[, 1], 0, hdiD2[, 2], 0, lwd=3, col='blue')

### 2013 STI

## 50%

plot.new()
par(mfrow=c(1,1))

hist(density_2013$mode.2013,freq = FALSE)

hdi<-hdi(density_2013$mode.2013, credMass=0.5)
hdi

segments(hdi[1], 0, hdi[2], 0, lwd=3, col='red')

dens2 <- density(density_2013$mode.2013)
lines(dens2, lwd=2, col='blue')

(hdiD2 <- hdi(dens2,credMass=0.5, allowSplit=TRUE))
segments(hdiD2[,1], 0, hdiD2[,2], 0, lwd=3, col='blue')

## 75%

plot.new()
par(mfrow=c(1,1))

hist(density_2013$mode.2013,freq = FALSE)

hdi<-hdi(density_2013$mode.2013, credMass=0.75)
hdi

segments(hdi[1], 0, hdi[2], 0, lwd=3, col='red')

dens2 <- density(density_2013$mode.2013)
lines(dens2, lwd=2, col='blue')

(hdiD2 <- hdi(dens2,credMass=0.75, allowSplit=TRUE))
segments(hdiD2[, 1], 0, hdiD2[, 2], 0, lwd=3, col='blue')

## 95%

plot.new()
par(mfrow=c(1,1))

hist(density_2013$mode.2013,freq = FALSE)

hdi<-hdi(density_2013$mode.2013, credMass=0.95)
hdi

segments(hdi[1], 0, hdi[2], 0, lwd=3, col='red')

dens2 <- density(density_2013$mode.2013)
lines(dens2, lwd=2, col='blue')

(hdiD2 <- hdi(dens2,credMass=0.95, allowSplit=TRUE))
segments(hdiD2[, 1], 0, hdiD2[, 2], 0, lwd=3, col='blue')

### 2019 STI

## 50%

plot.new()
par(mfrow=c(1,1))

hist(density_2019$mode.2019,freq = FALSE)

hdi<-hdi(density_2019$mode.2019, credMass=0.5)
hdi

segments(hdi[1], 0, hdi[2], 0, lwd=3, col='red')

dens2 <- density(density_2019$mode.2019)
lines(dens2, lwd=2, col='blue')

(hdiD2 <- hdi(dens2,credMass=0.5, allowSplit=TRUE))
segments(hdiD2[,1], 0, hdiD2[,2], 0, lwd=3, col='blue')

## 75%

plot.new()
par(mfrow=c(1,1))

hist(density_2019$mode.2019,freq = FALSE)

hdi<-hdi(density_2019$mode.2019, credMass=0.75)
hdi

segments(hdi[1], 0, hdi[2], 0, lwd=3, col='red')

dens2 <- density(density_2019$mode.2019)
lines(dens2, lwd=2, col='blue')

(hdiD2 <- hdi(dens2,credMass=0.75, allowSplit=TRUE))
segments(hdiD2[, 1], 0, hdiD2[, 2], 0, lwd=3, col='blue')

## 95%

plot.new()
par(mfrow=c(1,1))

hist(density_2019$mode.2019,freq = FALSE)

hdi<-hdi(density_2019$mode.2019, credMass=0.95)
hdi

segments(hdi[1], 0, hdi[2], 0, lwd=3, col='red')

dens2 <- density(density_2019$mode.2019)
lines(dens2, lwd=2, col='blue')

(hdiD2 <- hdi(dens2,credMass=0.95, allowSplit=TRUE))
segments(hdiD2[, 1], 0, hdiD2[, 2], 0, lwd=3, col='blue')

######################################################################################################
#################################### (3) TD and FD compositional dissimilarities #####################

setwd(data.dir)
dir()

traits<-read.csv('traits_herbivores.csv')%>%
  filter(!is.na(ThermalMP_5_95))%>%
  #dplyr::select(scientific,herbivore.guild,ThermalMP_5_95,MaxLength)%>%
  filter(!is.na(herbivore.guild))%>% ## 27 nominal herbivores
  glimpse()

## Load density data

x<-as.vector(unique(traits$scientific))

setwd(data.dir)
dir()

dat<-read.csv('dovs_clean.csv')%>%
  filter(!scientific%in%c('Ophisternon bengalense'))%>%
  filter(!Class%in%c('Elasmobranchii'))%>%
  filter(!scientific%in%c('Cirripectes hutchinsi'))%>%
  #filter(year!=2017)%>%
  filter(!Site%in%c('PGS1'))%>%
  filter(scientific%in%x)%>%
  glimpse()
unique(dat$Site)
unique(dat$year)
n_distinct(dat$scientific) ## Of the 27 - only 12 species observed in s-DOVs - but note some species were pooled for analyses
unique(dat$scientific)

## Filter out species not present from trait data

x<-as.vector(unique(dat$scientific))

traits<-traits%>%
  filter(scientific%in%x)%>%
  glimpse()

## Change nomenclature to match matrix

traits$scientific<-gsub(' ','.',traits$scientific)

### Check for missing trait information

sum(is.na(traits))/prod(dim(traits))*100
apply(traits,2,function(col)sum(is.na(col))/length(col))*100

## Convert first column to row names

traits<-data.frame(traits)
rownames(traits) <- traits[,1] #Assigning row names from 1st column 
traits[,1] <- NULL #Removing the first column

## See sampling effort per year and site

test<-dat%>%
  dplyr::group_by(year,Site)%>%
  dplyr::summarise(N=n_distinct(id))%>%
  glimpse()

## Create a grouping factor for each time period

dat<-dat%>%
  mutate(time.period=ifelse(year%in%c('2006'),"pre-heatwave",
                            ifelse(year%in%c('2013'),"post-heatwave","post-cooling")))%>%
  glimpse()

unique(dat$time.period)

dat<-dat%>%
  mutate(time.period=fct_relevel(time.period,"pre-heatwave","post-heatwave","post-cooling"))%>%
  glimpse()

## Resample data to select 8 transects per time period

## Firt we convert to long format for each year separately 

library(tidyr)

dataframe=data.frame()

number<-1:99
time.period<-as.vector(unique(dat$time.period))
site<-as.vector(unique(dat$Site))

for(j in number) {
  
  for (i in time.period) {
    
    for(k in site) {
      
      Abun<-dat%>%
        filter(time.period%in%i)%>%
        filter(Site%in%k)%>%
        dplyr::group_by(id,scientific)%>%
        summarise(number=sum(number))%>%
        ungroup()%>%
        pivot_wider(names_from = "scientific",values_from = number) %>%
        mutate_all(~replace_na(., 0))%>%
        sample_n(8,replace = T)%>%
        glimpse()
      
      ## Reconvert to long format
      
      Abun<-Abun%>%
        pivot_longer(cols = c(2:10), names_to = "scientific", values_to = "number")%>%
        glimpse()
      
      ## Create metadata 
      
      metadata<-dat%>%
        dplyr::select(Location,Site,time.period,id)%>%
        glimpse()
      
      new.dat<-left_join(Abun,metadata,by="id")%>%
        distinct(id,scientific,.keep_all = TRUE)%>%
        glimpse()
      
      dataframe<-rbind(dataframe,new.dat)
      
    }
  }
}

## Group data and obtain mean 

test<-dataframe%>%
  group_by(time.period,Site,scientific)%>%
  summarise(number=mean(number))%>%
  glimpse()

## Recheck sampling effort

dat<-test

## Clean working directory

rm(Abun,metadata,new.dat,dataframe)

## We use the Attribute diversity framework proposed by Anne Chao ----
## This allow us to test whether dissimilarites are mainly driven by changing compositions or abundance-structure

## Calculate functional distances between species - 
## We would use a Gower distance - to also account for changes in the guild of herbivores
## !! Note - we will need to retweak the code to adjust to different modal lengths across years

names(traits)

gowmat <- as.matrix(gowdis(traits[,c(2,3,8)], ord = "podani"))
min(gowmat)
mean(gowmat)
max(gowmat)

## Export Gower distance matrix

setwd(data.dir)
dir()

write.csv(gowmat,"gowmat.csv")

## Run in a loop for each Site

dataframe=data.frame()

site<-as.list(unique(dat$Site))

for(i in site) {

## Create a species x abundance matrix - with species as columns and years as rows

Abun<-as.data.frame(dat)%>%
  filter(Site%in%i)%>%
  dplyr::select(time.period,scientific,number)%>%
  pivot_wider(names_from = "scientific",values_from = number) %>%
  mutate_all(~replace_na(., 0))%>%
  glimpse()

## Set years as row names

Abun<-data.frame(Abun)
rownames(Abun) <- Abun[,1] #Assigning row names from 1st column 
Abun[,1] <- NULL #Removing the first column

#computing relative abundances, hence all communities have same total abundance
Abun<-Abun/apply(Abun,1,sum)
round(Abun,2)

sum(Abun[1,])

## Convert abundances to a list of N sites with abundances as vectors

Abun<-as.matrix(Abun)

## Calculate Sorensen's

## TD 

TD_sor<-beta.fd.hill(asb_sp_w=Abun, sp_dist=gowmat,
                          q=c(0,1,2), tau="min", beta_type="Sorensen")

TD_sor<-dist.to.df(TD_sor$beta_fd_q)

TD_sor<-TD_sor%>%
  mutate(year=paste(x1,x2,sep = "-"))%>%
  dplyr::select(-x1,-x2)%>%
  glimpse()

TD_sor<-TD_sor%>%
  mutate(metric="TD",
         Measurements="Sorensen")

## FD

FD_sor<-beta.fd.hill(asb_sp_w=Abun, sp_dist=gowmat,
                          q=c(0,1,2), tau="mean", beta_type="Sorensen")

FD_sor<-dist.to.df(FD_sor$beta_fd_q)

FD_sor<-FD_sor%>%
  mutate(year=paste(x1,x2,sep = "-"))%>%
  dplyr::select(-x1,-x2)%>%
  glimpse()

FD_sor<-FD_sor%>%
  mutate(metric="FD",
         Measurements="Sorensen")

## Bind rows

Diversity.beta<-bind_rows(TD_sor,FD_sor)%>%
  glimpse()

## Organize columns

names(Diversity.beta)

x<-c("metric","Measurements","year","q0","q1","q2")

Diversity.beta<-Diversity.beta[,x]

## Convert to long format

Diversity.beta<-Diversity.beta%>%
  pivot_longer(cols = c(q0,q1,q2), names_to = "q", values_to = "response")%>%
  glimpse()

Diversity.beta$Site<-i

dataframe<-rbind(dataframe,Diversity.beta)

}

## Export csv

setwd(data.dir)
dir()

write.csv(dataframe,"Beta.Abundance.csv")

## Bring in dataframe

setwd(data.dir)
dir()

dat<-read.csv('Beta.Abundance.csv')%>%
  glimpse()

## Get summaries per year

bbs.summary<-dat%>%
  group_by(metric,year,q)%>%
  summarise(mean=mean(response),
            se=se(response))%>%
  glimpse()

## Plot of TD dissimilarities

glimpse(bbs.summary)

bbs.summary<-bbs.summary%>%
  mutate(time.period=fct_relevel(year,"pre-heatwave-post-heatwave","pre-heatwave-post-cooling","post-heatwave-post-cooling"))%>%
  glimpse()
label<-c("pre-heatwave-post-heatwave","pre-heatwave-post-cooling","post-heatwave-post-cooling")

TD.plot<-ggplot(bbs.summary%>%filter(metric%in%c('TD')),aes(x=q,y=mean,fill=q))+
  geom_col(colour="black")+
  geom_errorbar(aes(ymin=mean, ymax=mean+se),color = "black",width=.2,size=1)+
  ggtitle('')+
  ylab('Taxonomic dissimilarities')+
  xlab('')+
  scale_y_continuous(expand=c(0,0))+
  scale_fill_manual(labels = c("q0","q1","q2"),
                    values=c("black","darkgrey","lightgrey"))+
  facet_wrap(~time.period)+
  theme_classic()
TD.plot

## Export plot 

setwd(plots.dir)
dir()

name<-'Beta_TD_Abundance'

ggsave(paste(name,".png",sep="."),width = 21, height = 7,units = "cm",dpi=600)
ggsave(paste(name,".pdf",sep="."),width = 21, height = 7,units = "cm",dpi=600,useDingbats=FALSE)

## Plot of FD dissimilarities

glimpse(bbs.summary)

FD.plot<-ggplot(bbs.summary%>%filter(metric%in%c('FD')),aes(x=q,y=mean,fill=q))+
  geom_col(colour="black")+
  geom_errorbar(aes(ymin=mean, ymax=mean+se),color = "black",width=.2,size=1)+
  ggtitle('')+
  ylab('Functional dissimilarities')+
  xlab('')+
  scale_y_continuous(expand=c(0,0))+
  scale_fill_manual(labels = c("q0","q1","q2"),
                    values=c("black","darkgrey","lightgrey"))+
  facet_wrap(~time.period)+
  theme_classic()
FD.plot

## Export plot 

setwd(plots.dir)
dir()

name<-'Beta_FD_Abundance'

ggsave(paste(name,".png",sep="."),width = 21, height = 7,units = "cm",dpi=600)
ggsave(paste(name,".pdf",sep="."),width = 21, height = 7,units = "cm",dpi=600,useDingbats=FALSE)

## Heatmap of taxonomic and functional dissimilarities ----

# colour ramps-
re <- colorRampPalette(c("mistyrose", "red2","darkred"))(200)
bl <- colorRampPalette(c("grey62","grey14","black"))(200)  
glimpse(bbs.summary)

# Reorder levels of the factor

unique(bbs.summary$time.period)

bbs.summary<-bbs.summary%>%
  dplyr::mutate(time.period=fct_relevel(time.period,'post-heatwave-post-cooling',
                                        'pre-heatwave-post-cooling',
                                        'pre-heatwave-post-heatwave'))%>%
  glimpse()
  
## Plot heatmap-

## TD -

ggheat.TD <- ggplot(bbs.summary%>%filter(metric%in%c('TD')), aes(x=time.period,y=q,fill=mean))+ 
  ggtitle('(a) TD')+
  geom_tile(show.legend=T) + 
  scale_fill_gradient(low="white",high=bl,labels=c('0','0.25','0.5','0.75','1'))+
  xlab(NULL)+
  ylab(NULL)+
  scale_y_discrete(expand=c(0,0))+
  scale_x_discrete(expand=c(0,0))+
  coord_flip()+
  Theme1+
  theme(axis.text.x = element_text(angle = 45,hjust = 1,vjust=1,size=10,colour = "black"),
        axis.text.y = element_text(angle = 0,hjust = 1,vjust=1,size=10,colour = "black"),
        legend.spacing.x = unit(1, 'cm'),
        legend.title=element_text(size = 14),
        legend.text=element_text(size=14),
        legend.key.size = unit(1, "cm"),
        legend.position = "right",
        legend.direction = "vertical")
ggheat.TD

## Export plot

setwd(plots.dir)
dir()

name<-'Fig2_Heapmap_TD'

ggsave(paste(name,".png",sep="."),width = 10.5, height = 10,units = "cm",dpi = 600)
ggsave(paste(name,".pdf",sep="."),width = 10.5, height = 10,units = "cm",dpi = 600,useDingbats=FALSE)

## FD -

ggheat.FD <- ggplot(bbs.summary%>%filter(metric%in%c('FD')), aes(x=time.period,y=q,fill=mean))+ 
  ggtitle('(b) FD')+
  geom_tile(show.legend=T) + 
  scale_fill_gradient(low="white",high=bl,labels=c('0','0.25','0.5','0.75','1'))+
  xlab(NULL)+
  ylab(NULL)+
  scale_y_discrete(expand=c(0,0))+
  scale_x_discrete(expand=c(0,0))+
  coord_flip()+
  Theme1+
  theme(axis.text.x = element_text(angle = 45,hjust = 1,vjust=1,size=10,colour = "black"),
        axis.text.y = element_text(angle = 0,hjust = 1,vjust=1,size=10,colour = "black"),
        legend.spacing.x = unit(1, 'cm'),
        legend.title=element_text(size = 14),
        legend.text=element_text(size=14),
        legend.key.size = unit(1, "cm"),
        legend.position = "none")
ggheat.FD

## Export plot

setwd(plots.dir)
dir()

name<-'Fig2_Heapmap_FD'

ggsave(paste(name,".png",sep="."),width = 10.5, height = 10,units = "cm",dpi = 600)
ggsave(paste(name,".pdf",sep="."),width = 10.5, height = 10,units = "cm",dpi = 600,useDingbats=FALSE)


## (3.4) GLM testing for differences in dissimilarities ----

## Bring in data

setwd(data.dir)
dir()

beta.TD<-read.csv('beta_TD.csv')%>%
  glimpse()

beta.FD<-read.csv('Beta_FD_2D.csv')%>%
  glimpse()

beta.abundance<-read.csv('Beta.Abundance.csv')%>%
  glimpse()

## Split abundance-based dissimilarities in TD and FD

beta.TD.abundance<-beta.abundance%>%
  filter(metric%in%('TD'))%>%
  glimpse()

beta.FD.abundance<-beta.abundance%>%
  filter(metric%in%('FD'))%>%
  glimpse()

## Spread q values for modelling

beta.TD.abundance<-beta.TD.abundance%>%
  pivot_wider(names_from = "q",values_from = response) %>%
  mutate_all(~replace_na(., 0))%>%
  glimpse()

beta.FD.abundance<-beta.FD.abundance%>%
  pivot_wider(names_from = "q",values_from = response) %>%
  mutate_all(~replace_na(., 0))%>%
  glimpse()

## GLMs models

## (A) TD - incidence-based ----

## (A.1) Total

glimpse(beta.TD)

hist(beta.TD$beta.sor)

mod<-lm(beta.sor ~ year,
             data = beta.TD)
summary(mod)
car::Anova(mod)

## Check residual pattern

library(DHARMa)

mod.res <- simulateResiduals(mod)
plot(mod.res)

## Pairwise comparisons

library(emmeans)

emmeans(mod, pairwise ~ year)

## (A.2) Turnover

glimpse(beta.TD)

hist(beta.TD$beta.sim)

mod<-lm(beta.sim ~ year,
        data = beta.TD)
summary(mod)
car::Anova(mod)

## Check residual pattern

library(DHARMa)

mod.res <- simulateResiduals(mod)
plot(mod.res)

## Pairwise comparisons

library(emmeans)

emmeans(mod, pairwise ~ year)

## (A.3) Nestedness

glimpse(beta.TD)

hist(beta.TD$beta.sne)

mod<-lm(beta.sne ~ year,
        data = beta.TD)
summary(mod)
car::Anova(mod)

## Check residual pattern

library(DHARMa)

mod.res <- simulateResiduals(mod)
plot(mod.res)

## Pairwise comparisons

library(emmeans)

emmeans(mod, pairwise ~ year)

## (B) FD - incidence-based ----

## (A.1) Total

glimpse(beta.FD)

hist(beta.FD$beta.sor)

mod<-lm(beta.sor ~ year,
        data = beta.FD)
summary(mod)
car::Anova(mod)

## Check residual pattern

library(DHARMa)

mod.res <- simulateResiduals(mod)
plot(mod.res)

## Pairwise comparisons

library(emmeans)

emmeans(mod, pairwise ~ year)

## (A.2) Turnover

glimpse(beta.FD)

hist(beta.FD$beta.sim)

mod<-lm(beta.sim ~ year,
        data = beta.FD)
summary(mod)
car::Anova(mod)

## Check residual pattern

library(DHARMa)

mod.res <- simulateResiduals(mod)
plot(mod.res)

## Pairwise comparisons

library(emmeans)

emmeans(mod, pairwise ~ year)

## (A.3) Nestedness

glimpse(beta.FD)

hist(beta.FD$beta.sne)

mod<-lm(beta.sne ~ year,
        data = beta.FD)
summary(mod)
car::Anova(mod)

## Check residual pattern

library(DHARMa)

mod.res <- simulateResiduals(mod)
plot(mod.res)

## Pairwise comparisons

library(emmeans)

emmeans(mod, pairwise ~ year)

## (C) TD - abundance-based ----

## (A.1) "q" = 0

glimpse(beta.TD.abundance)

hist(beta.TD.abundance$q0)

mod<-lm(q0 ~ year,
        data = beta.TD.abundance)
summary(mod)
car::Anova(mod)

## Check residual pattern

library(DHARMa)

mod.res <- simulateResiduals(mod)
plot(mod.res)

## Pairwise comparisons

library(emmeans)

emmeans(mod, pairwise ~ year)

## (A.2) "q" = 1

glimpse(beta.TD.abundance)

hist(beta.TD.abundance$q1)

mod<-lm(q1 ~ year,
        data = beta.TD.abundance)
summary(mod)
car::Anova(mod)

## Check residual pattern

library(DHARMa)

mod.res <- simulateResiduals(mod)
plot(mod.res)

## Pairwise comparisons

library(emmeans)

emmeans(mod, pairwise ~ year)

## (A.3) q2

glimpse(beta.TD.abundance)

hist(beta.TD.abundance$q2)

mod<-lm(q2 ~ year,
        data = beta.TD.abundance)
summary(mod)
car::Anova(mod)

## Check residual pattern

library(DHARMa)

mod.res <- simulateResiduals(mod)
plot(mod.res)

## Pairwise comparisons

library(emmeans)

emmeans(mod, pairwise ~ year)

## (D) FD - abundance-based ----

## (D.1) "q" = 0

glimpse(beta.FD.abundance)

hist(beta.FD.abundance$q0)

mod<-lm(q0 ~ year,
        data = beta.FD.abundance)
summary(mod)
car::Anova(mod)

## Check residual pattern

library(DHARMa)

mod.res <- simulateResiduals(mod)
plot(mod.res)

## Pairwise comparisons

library(emmeans)

emmeans(mod, pairwise ~ year)

## (A.2) "q" = 1

glimpse(beta.FD.abundance)

hist(beta.FD.abundance$q1)

mod<-lm(q1 ~ year,
        data = beta.FD.abundance)
summary(mod)
car::Anova(mod)

## Check residual pattern

library(DHARMa)

mod.res <- simulateResiduals(mod)
plot(mod.res)

## Pairwise comparisons

library(emmeans)

emmeans(mod, pairwise ~ year)

## (A.3) q2

glimpse(beta.FD.abundance)

hist(beta.FD.abundance$q2)

mod<-lm(q2 ~ year,
        data = beta.FD.abundance)
summary(mod)
car::Anova(mod)

## Check residual pattern

library(DHARMa)

mod.res <- simulateResiduals(mod)
plot(mod.res)

## Pairwise comparisons

library(emmeans)

emmeans(mod, pairwise ~ year)

#########################################################################################################################
############################### (4) GLMM - Changes in abundances and consumption ########################################

## (A) DOVs ----

setwd(data.dir)
dir()

dat<-read.csv('dovs_clean.csv')%>%
  filter(!scientific%in%c('Ophisternon bengalense'))%>%
  filter(!Class%in%c('Elasmobranchii'))%>%
  #filter(year!=2017)%>%
  filter(!Site%in%c('PGS1'))%>%
  glimpse()
n_distinct(dat$scientific) ## 65 species
unique(dat$scientific)

## Quick check on samples,year, site for unbalanced design

test<-dat%>%
  dplyr::group_by(year,Site)%>%
  dplyr::summarise(N=n_distinct(id))%>%
  glimpse()

## Bring in traits and filter only herbivores

traits<-read.csv('traits_herbivores.csv')%>%
  filter(!is.na(herbivore.guild))%>%
  filter(!is.na(ThermalMP_5_95))%>%
  dplyr::select(scientific,Realised.function,herbivore.guild,ThermalMP_5_95,MaxLength)%>%
  mutate(thermal.guild=ifelse(ThermalMP_5_95>=23,"Tropical","Temperate"))%>%
  mutate(RF.TG=paste(Realised.function,thermal.guild,sep="_"),
         HG.TG=paste(herbivore.guild,thermal.guild,sep="_"))%>%
  glimpse()

## Filter only herbivores from the data

x<-as.vector(unique(traits$scientific))

dat<-dat%>%
  filter(scientific%in%x)%>%
  glimpse()
n_distinct(dat$scientific) ## 10 species
unique(dat$scientific)

## Create a grouping factor for each time period

dat<-dat%>%
  mutate(time.period=ifelse(year%in%c('2006'),"pre-heatwave",
                            ifelse(year%in%c('2013'),"post-heatwave","post-cooling")))%>%
  glimpse()

unique(dat$time.period)

dat<-dat%>%
  mutate(time.period=fct_relevel(time.period,"pre-heatwave","post-heatwave","post-cooling"))%>%
  glimpse()

## Supplementary plot on mean +- se densities

dat.summary<-dat%>%
  group_by(time.period,scientific)%>%
  summarise(mean=mean(number),
            se=se(number))%>%
  glimpse()

## Abundance-rank plot

plot.pre.heatwave<-ggplot(dat.summary%>%filter(time.period%in%c('pre-heatwave')),aes(x=scientific,y=mean))+
  geom_col(fill="lightgrey",colour="black")+
  geom_errorbar(aes(ymin=mean, ymax=mean+se),color = "black",width=.2,size=1)+
  ggtitle('(a)')+
  ylab('Densities (ind./125m2)')+
  xlab('')+
  scale_y_continuous(expand=c(0,0),limits = c(0,16))+
  theme_classic()+
  theme(axis.text.x = element_text(angle=45,hjust=1,vjust=1))
plot.pre.heatwave

plot.post.heatwave<-ggplot(dat.summary%>%filter(time.period%in%c('post-heatwave')),aes(x=scientific,y=mean))+
  geom_col(fill="lightgrey",colour="black")+
  geom_errorbar(aes(ymin=mean, ymax=mean+se),color = "black",width=.2,size=1)+
  ggtitle('(b)')+
  ylab('')+
  xlab('')+
  scale_y_continuous(expand=c(0,0),limits = c(0,16))+
  theme_classic()+
  theme(axis.text.x = element_text(angle=45,hjust=1,vjust=1))
plot.post.heatwave

plot.post.cooling<-ggplot(dat.summary%>%filter(time.period%in%c('post-cooling')),aes(x=scientific,y=mean))+
  geom_col(fill="lightgrey",colour="black")+
  geom_errorbar(aes(ymin=mean, ymax=mean+se),color = "black",width=.2,size=1)+
  ggtitle('(c)')+
  ylab('')+
  xlab('')+
  scale_y_continuous(expand=c(0,0),limits = c(0,16))+
  theme_classic()+
  theme(axis.text.x = element_text(angle=45,hjust=1,vjust=1))
plot.post.cooling

## Arrange plots in a grob

library(ggpubr)

ggarrange(plot.pre.heatwave,plot.post.heatwave,plot.post.cooling,
          nrow=1,ncol=3,align = 'hv')

## Export plot

setwd(plots.dir)
dir()

name<-'Appendix_FigureA3_Dominance_DOVs'

ggsave(paste(name,".png",sep="."),width = 21, height = 10,units = "cm")
ggsave(paste(name,".pdf",sep="."),width = 21, height = 10,units = "cm",dpi=600,useDingbats=FALSE)

## Join trait information and summarize data - by trophic guild - thermal guild combination?

dat<-left_join(dat,traits,by="scientific")%>%
  glimpse()
n_distinct(dat$RF.TG) # 5 groups
unique(dat$RF.TG)

dat.sum<-dat%>%
  filter(!herbivore.guild%in%c('Excavator'))%>%
  dplyr::group_by(Location,Site,time.period,year,id,Realised.function,RF.TG,thermal.guild)%>%
  dplyr::summarise(response=sum(number))%>%
  glimpse()

n_distinct(dat$id)

## Check distribution of densities across years

ggplot(dat.sum,aes(x=reorder(RF.TG,-response),y=response,fill=RF.TG))+
  geom_bar(stat = "identity",colour="black",show.legend = F)+
  facet_wrap(~time.period,scales = "free")+
  theme_classic()+
  theme(axis.text.x = element_text(angle=90))

ggplot(dat.sum,aes(x=reorder(RF.TG,-response),y=log10(response+1),fill=RF.TG))+
  geom_bar(stat = "identity",colour="black",show.legend = F)+
  facet_wrap(~time.period,scales = "free")+
  theme_classic()+
  theme(axis.text.x = element_text(angle=90))

## Run GLMMs on HG.TG and save coefficients =- se

## Create a loop to investigate the distribution across HG.TG groups
## All models with a negative binomial distribution

n_distinct(dat.sum$Realised.function) ## 3 groups
unique(dat.sum$Realised.function)

pred.vars.cont=c("Macroalgal removal","Turf + sediment removal","Turf removal")

plot.new()
par(mfrow=c(3,2))

for (i in pred.vars.cont) {
  x<-dat.sum%>%filter(Realised.function%in%i)
  hist(x$response,main = paste(i))#Looks best
}

## Run models in a loop

i<-c("Turf + sediment removal","Turf removal")

# Create an empty dataframe

dataframe=data.frame()

for (i in pred.vars.cont) {

dat.mod<-dat.sum%>%
  filter(Realised.function%in%i)%>%
  glimpse()
unique(dat.mod$Realised.function)
unique(dat.mod$thermal.guild)

## Recode factors for modelling

dat.mod$time.period<-as.factor(dat.mod$time.period)
dat.mod$thermal.guild<-as.factor(dat.mod$thermal.guild)

## Investigate sample size distribution

plot(dat.mod$time.period)
plot(dat.mod$thermal.guild)

## Run model

library(glmmTMB)

mod<-glmmTMB(response ~ time.period*thermal.guild + (1|Site),
             data = dat.mod,family = tweedie())
summary(mod)
car::Anova(mod)

system.time(mod_d1 <- drop1(mod,test="Chisq",all.cols=TRUE))
print(mod_d1)

## Pairwise comparisons

emmeans(mod, pairwise ~ time.period:thermal.guild)

## Calculate R2 = 1 - (Residual Deviance/Null Deviance)

library(performance)

r2_nakagawa(mod)

## Check residual pattern

library(DHARMa)

mod.res <- simulateResiduals(mod)
plot(mod.res)

## Test zero-inflation

testZeroInflation(mod.res)

## Plot the effect of gravity by management status

library(sjPlot)

plot_model(mod,type = "eff",terms = c("time.period"))

## Predict values from model and store in a csv

unique(dat.mod$thermal.guild)
unique(dat.mod$time.period)

testdata <-expand.grid(time.period=c('pre-heatwave','post-heatwave','post-cooling'),
                       thermal.guild=c('Temperate','Tropical'),
                       Site="PGN1")%>%
  glimpse()

fits <- predict(mod, newdata=testdata, type = "response",se.fit=T,re.form = NA)

## Group predictions

predictions = testdata%>%data.frame(fits)%>%
  group_by(time.period,thermal.guild)%>% #only change here
  dplyr::summarise(response=mean(fit),se.fit=mean(se.fit))%>%
  ungroup()%>%
  glimpse()
predictions$Realised.function<-i

dataframe<-rbind(dataframe,predictions)

}

## Export dataframe 

setwd(data.dir)
dir()

write.csv(dataframe,"Estimates_test.csv")

## Bring abundance estimates 

setwd(data.dir)
dir()

dataframe<-read.csv("Estimates_Abundance.csv")%>%
  glimpse()

## Plot in ggplot

glimpse(dat.sum)
glimpse(dataframe)
unique(dat.sum$Realised.function)

dat.sum$time.period<-as.factor(dat.sum$time.period)
dataframe$time.period<-as.factor(dataframe$time.period)

## Re-level factors year

unique(dat.sum$thermal.guild)

dat.sum<-dat.sum%>%
  mutate(thermal.guild=fct_relevel(thermal.guild,"Temperate","Tropical"))%>%
  glimpse()

dataframe<-dataframe%>%
  mutate(thermal.guild=fct_relevel(thermal.guild,"Temperate","Tropical"))%>%
  glimpse()

unique(dataframe$Realised.function)

Macroalgae<-ggplot(data=dataframe%>%filter(Realised.function%in%c('Macroalgal removal')),aes(x=time.period,y=response))+
  ggtitle('(a)')+
  ylab('Densities (individuals/125m2)')+
  xlab('')+
  geom_jitter(data=dat.sum%>%filter(Realised.function%in%c('Macroalgal removal')),aes(x=time.period,y=response,colour=thermal.guild),position=position_jitterdodge(jitter.width = 0.2, jitter.height = 0,
                                                                                                  dodge.width = 0.75, seed = NA), size=2,alpha=0.25,show.legend = F)+
  geom_errorbar(data = dataframe%>%filter(Realised.function%in%c('Macroalgal removal')),aes(x = time.period,y=response,
                                    ymin=response-se.fit, ymax=response+se.fit,fill=thermal.guild),color = "black",width=.2,
                position=position_jitterdodge(jitter.width = 0, jitter.height = 0,dodge.width = 0.75, seed = NA), size=1)+
  geom_point(data=dataframe%>%filter(Realised.function%in%c('Macroalgal removal')),
             aes(x=time.period,y=response,fill=thermal.guild),color="black",
             size=4.5,alpha=2,
             position=position_dodge(.75),show.legend = F)+
  geom_point(data = dataframe%>%filter(Realised.function%in%c('Macroalgal removal')), aes(x=time.period,y=response,colour=thermal.guild),
             size=4,position=position_dodge(.75),show.legend = F)+
  geom_path(data=dataframe%>%filter(Realised.function%in%c('Macroalgal removal')),
            aes(x=time.period,y=response,group = thermal.guild,color=thermal.guild),
            position=position_dodge(.75),
            size=1,show.legend = F)+
  scale_y_continuous(limits = c(0,400),trans="log1p")+
  scale_colour_manual(labels = c("Temperate","Tropical"),
                    values=c("blue","red"))+
  theme_classic()+
  Theme1+
  theme(legend.position = "bottom",
        legend.direction = "horizontal")
Macroalgae

Turf<-ggplot(data=dataframe%>%filter(Realised.function%in%c('Turf + sediment removal','Turf removal')),aes(x=time.period,y=response))+
  ggtitle('(b)')+
  ylab('')+
  xlab('')+
  geom_jitter(data=dat.sum%>%filter(Realised.function%in%c('Turf + sediment removal','Turf removal')),aes(x=time.period,y=response,colour=thermal.guild),position=position_jitterdodge(jitter.width = 0.2, jitter.height = 0,
                                                                                                                                               dodge.width = 0.75, seed = NA), size=2,alpha=0.25,show.legend = F)+
  geom_errorbar(data = dataframe%>%filter(Realised.function%in%c('Turf + sediment removal','Turf removal')), aes(x = time.period,y=response,
                                                                                ymin=response-se.fit, ymax=response+se.fit,fill=thermal.guild),color = "black",width=.2,
                position=position_jitterdodge(jitter.width = 0, jitter.height = 0,dodge.width = 0.75, seed = NA), size=1)+
  geom_point(data=dataframe%>%filter(Realised.function%in%c('Turf + sediment removal','Turf removal')),
             aes(x=time.period,y=response,fill=thermal.guild),color="black",
             size=4.5,alpha=2,
             position=position_dodge(.75),show.legend = F)+
  geom_point(data = dataframe%>%filter(Realised.function%in%c('Turf + sediment removal','Turf removal')), aes(x=time.period,y=response,colour=thermal.guild),
             size=4,position=position_dodge(.75),show.legend = F)+
  geom_path(data=dataframe%>%filter(Realised.function%in%c('Turf + sediment removal','Turf removal')),
            aes(x=time.period,y=response,group = thermal.guild,color=thermal.guild),
            position=position_dodge(.75),
            size=1,show.legend = F)+
  scale_y_continuous(limits=c(0,400),trans="log1p")+
  scale_colour_manual(labels = c("Temperate","Tropical"),
                      values=c("blue","red"))+
  theme_classic()+
  Theme1+
  theme(legend.position = "bottom",
        legend.direction = "horizontal",
        axis.title.y = element_blank())
Turf

Sediments<-ggplot(data=dataframe%>%filter(Realised.function%in%c('Sediment removal')),aes(x=time.period,y=response))+
  ggtitle('(c)')+
  ylab('')+
  xlab('')+
  geom_jitter(data=dat.sum%>%filter(Realised.function%in%c('Turf + sediment removal')),aes(x=time.period,y=response,colour=thermal.guild),position=position_jitterdodge(jitter.width = 0.2, jitter.height = 0,
                                                                                                                                                                                       dodge.width = 0.75, seed = NA), size=2,alpha=0.25,show.legend = F)+
  geom_errorbar(data = dataframe%>%filter(Realised.function%in%c('Sediment removal')), aes(x = time.period,y=response,
                                                                                                                 ymin=response-se.fit, ymax=response+se.fit,fill=thermal.guild),color = "black",width=.2,
                position=position_jitterdodge(jitter.width = 0, jitter.height = 0,dodge.width = 0.75, seed = NA), size=1)+
  geom_point(data=dataframe%>%filter(Realised.function%in%c('Sediment removal')),
             aes(x=time.period,y=response,fill=thermal.guild),color="black",
             size=4.5,alpha=2,
             position=position_dodge(.75),show.legend = F)+
  geom_point(data = dataframe%>%filter(Realised.function%in%c('Sediment removal')), aes(x=time.period,y=response,colour=thermal.guild),
             size=4,position=position_dodge(.75),show.legend = F)+
  geom_path(data=dataframe%>%filter(Realised.function%in%c('Sediment removal')),
            aes(x=time.period,y=response,group = thermal.guild,color=thermal.guild),
            position=position_dodge(.75),
            size=1,show.legend = F)+
  scale_y_continuous(limits=c(0,400),trans="log1p")+
  scale_colour_manual(labels = c("Tropical"),
                      values=c("red"))+
  theme_classic()+
  Theme1+
  theme(legend.position = "bottom",
        legend.direction = "horizontal",
        axis.title.y = element_blank())
Sediments

## Arrange plots in a grob

library(ggpubr)

ggarrange(Macroalgae,Turf,Sediments,
          nrow=1,ncol = 4,
          align = 'hv',
          common.legend = TRUE,
          legend = "bottom")

## Export plot

setwd(plots.dir)
dir()

name<-'DOVs_Herbivore_Change_17_05_2021'

ggsave(paste(name,".png",sep="."),width = 21, height = 7,units = "cm")
ggsave(paste(name,".pdf",sep="."),width = 21, height = 7,units = "cm",dpi=600,useDingbats=FALSE)

## Now we can calculate the relative proportion of each phenotype ('how' function) for each time period

prop<-dat%>%
  filter(!herbivore.guild%in%c('Excavator'))%>%
  dplyr::group_by(Location,Site,time.period,year,id,herbivore.guild,HG.TG,thermal.guild)%>%
  dplyr::summarise(response=sum(number))%>%
  glimpse()

unique(prop$HG.TG)

dataframe.1=data.frame()

period<-as.vector(unique(prop$time.period)) ## For n in 3 time periods

for (i in period) {
  
  ## Create frequency of occurrence for each latitudinal bin and frequency class
  
  test<-prop%>%
    filter(time.period%in%i)%>%
    glimpse()
  unique(test$time.period)
  n<-sum(test$response)
  
  test.1<-test%>%
    dplyr::group_by(time.period,HG.TG)%>%
    dplyr::summarise(number=sum(response))%>%
    mutate(percentage=(number/n)*100)%>%
    ungroup()%>%
    glimpse()
  
  dataframe.1<-rbind(dataframe.1,test.1)
  
}

## Plot of relative frequencies of common, frequent and rare species

library(forcats)

unique(dataframe.1$HG.TG)

dataframe.1<-dataframe.1%>%
  mutate(HG.TG=fct_relevel(HG.TG,"Browser_Temperate","Browser_Tropical","Algal farmer_Temperate","Algal farmer_Tropical",
                           "Scraper_Tropical","Sediment sucker_Tropical"))%>%
  glimpse()

proportion<-ggplot(dataframe.1,aes(fill=HG.TG, y=percentage, x=time.period)) + 
  ggtitle('(d)')+
  xlab('')+
  ylab('Contribution (%)')+
  geom_bar(position="stack", stat="identity",colour="black")+
  scale_fill_manual(labels=c("Browser_Temperate","Browser_Tropical","Algal farmer_Temperate","Algal farmer_Tropical",
                             "Scraper_Tropical","Sediment sucker_Tropical"),
                    values = c('lightsalmon4','lightsalmon1','green4','green1','darkorange1','darkorchid1'))+
  theme_classic()+
  Theme1
proportion

## Arrange plots in a grob

library(ggpubr)

ggarrange(Macroalgae,Turf,Sediments,proportion,
          nrow=1,ncol = 4,
          align = 'hv',
          common.legend = TRUE,
          legend = "none")

## Export plot

setwd(plots.dir)
dir()

name<-'DOVs_Herbivore_Change_17_05_2021'

ggsave(paste(name,".png",sep="."),width = 21, height = 7,units = "cm")
ggsave(paste(name,".pdf",sep="."),width = 21, height = 7,units = "cm",dpi=600,useDingbats=FALSE)

## (B) Turf ----

setwd(data.dir)
dir()

dat<-read.csv('turf_bites_clean.csv')%>%
  filter(!scientific%in%c('Ophisternon bengalense'))%>%
  filter(!Class%in%c('Elasmobranchii'))%>%
  glimpse()
n_distinct(dat$scientific) ## 10 species
unique(dat$scientific)

## Quick check on samples,year, site for unbalanced design

test<-dat%>%
  dplyr::group_by(year,site,date)%>%
  dplyr::summarise(N=n_distinct(id))%>%
  glimpse()
n_distinct(dat$id)

## Bring in traits and filter only herbivores

traits<-read.csv('traits_herbivores.csv')%>%
  filter(!is.na(herbivore.guild))%>%
  filter(!is.na(ThermalMP_5_95))%>%
  dplyr::select(scientific,Realised.function,herbivore.guild,ThermalMP_5_95,MaxLength)%>%
  mutate(thermal.guild=ifelse(ThermalMP_5_95>=23,"Tropical","Temperate"))%>%
  mutate(RF.TG=paste(Realised.function,thermal.guild,sep="_"),
         HG.TG=paste(herbivore.guild,thermal.guild,sep="_"))%>%
  glimpse()

## Filter only herbivores from the data

x<-as.vector(unique(traits$scientific))

dat<-dat%>%
  filter(scientific%in%x)%>%
  glimpse()
n_distinct(dat$scientific) ## 10 species
unique(dat$scientific)

## Supplementary plot on mean +- se densities

dat.summary<-dat%>%
  group_by(year,scientific)%>%
  summarise(mean=mean(mass.std.bites),
            se=se(mass.std.bites))%>%
  glimpse()

## Abundance-rank plot

plot.2013<-ggplot(dat.summary%>%filter(year%in%c('2013')),aes(x=scientific,y=mean))+
  geom_col(fill="lightgrey",colour="black")+
  geom_errorbar(aes(ymin=mean, ymax=mean+se),color = "black",width=.2,size=1)+
  ggtitle('(a)')+
  ylab('')+
  xlab('')+
  scale_y_continuous(expand=c(0,0),limits = c(0,200),trans = "log1p")+
  theme_classic()+
  theme(axis.text.x = element_text(angle=45,hjust=1,vjust=1))
plot.2013

plot.2019<-ggplot(dat.summary%>%filter(year%in%c('2019')),aes(x=scientific,y=mean))+
  geom_col(fill="lightgrey",colour="black")+
  geom_errorbar(aes(ymin=mean, ymax=mean+se),color = "black",width=.2,size=1)+
  ggtitle('(b)')+
  ylab('')+
  xlab('')+
  scale_y_continuous(expand=c(0,0),limits = c(0,200),trans = "log1p")+
  theme_classic()+
  theme(axis.text.x = element_text(angle=45,hjust=1,vjust=1))
plot.2019

## Arrange plots in a grob

library(ggpubr)

ggarrange(plot.2013,plot.2019,
          nrow=1,ncol=2,align = 'hv')

## Export plot

setwd(plots.dir)
dir()

name<-'Appendix_FigureA4_Dominance_Turf'

ggsave(paste(name,".png",sep="."),width = 21, height = 10,units = "cm")
ggsave(paste(name,".pdf",sep="."),width = 21, height = 10,units = "cm",dpi=600,useDingbats=FALSE)

## (B.1) Patterns of consumption for each realised ecosystem function 
## Join trait information and summarize data - by trophic guild - thermal guild combination?

dat<-left_join(dat,traits,by="scientific")%>%
  glimpse()
n_distinct(dat$RF.TG) # 5 groups
unique(dat$RF.TG)

dat.sum<-dat%>%
  filter(!herbivore.guild%in%c('Excavator'))%>%
  dplyr::group_by(location,site,year,date,id,Realised.function,RF.TG,thermal.guild)%>%
  dplyr::summarise(response=sum(mass.std.bites))%>%
  glimpse()

n_distinct(dat$id)

## Check distribution of densities across years

ggplot(dat.sum,aes(x=reorder(RF.TG,-response),y=response,fill=RF.TG))+
  geom_bar(stat = "identity",colour="black",show.legend = F)+
  facet_wrap(~year,scales = "free")+
  theme_classic()+
  theme(axis.text.x = element_text(angle=90))

ggplot(dat.sum,aes(x=reorder(RF.TG,-response),y=log10(response+1),fill=RF.TG))+
  geom_bar(stat = "identity",colour="black",show.legend = F)+
  facet_wrap(~year,scales = "free")+
  theme_classic()+
  theme(axis.text.x = element_text(angle=90))

## Run GLMMs on HG.TG and save coefficients =- se

## Create a loop to investigate the distribution across HG.TG groups
## All models with a negative binomial distribution

n_distinct(dat.sum$Realised.function) ## 3 groups
unique(dat.sum$Realised.function)

pred.vars.cont=c("Macroalgal removal","Turf + sediment removal","Turf removal")

plot.new()
par(mfrow=c(1,3))

for (i in pred.vars.cont) {
  x<-dat.sum%>%filter(Realised.function%in%i)
  hist(x$response,main = paste(i))#Looks best
}

## Run models in a loop

library(performance)
library(sjPlot)

i<-c("Turf + sediment removal","Turf removal")

# Create an empty dataframe

dataframe=data.frame()

for (i in pred.vars.cont) {
  
  dat.mod<-dat.sum%>%
    filter(Realised.function%in%i)%>%
    glimpse()
  unique(dat.mod$Realised.function)
  
  ## Recode factors for modelling
  
  dat.mod$year<-as.factor(dat.mod$year)
  dat.mod$thermal.guild<-as.factor(dat.mod$thermal.guild)
  dat.mod$date<-as.factor(dat.mod$date)
  
  ## Investigate sample size distribution
  
  plot(dat.mod$year)
  plot(dat.mod$thermal.guild)
  plot(dat.mod$date)
  
  ## Run model
  
  library(glmmTMB)
  library(MuMIn)
  
  mod<-glmmTMB(response ~ year*thermal.guild + (1|site) + (1|date),
               data = dat.mod,family = tweedie())
  summary(mod)
  car::Anova(mod)
  
  system.time(mod_d1 <- drop1(mod,test="Chisq",all.cols=TRUE))
  print(mod_d1)
  
  ## Pairwise comparisons
  
  emmeans(mod, pairwise ~ year:thermal.guild)
  
  ## Calculate R2 = 1 - (Residual Deviance/Null Deviance)
  
  r2_nakagawa(mod)
  
  ## Check residual pattern
  
  library(DHARMa)
  
  mod.res <- simulateResiduals(mod)
  plot(mod.res)
  
  ## Test zero-inflation
  
  testZeroInflation(mod.res)
  
  ## Plot the effect of gravity by management status
  
  plot_model(mod,type = "eff",terms = c("year"))
  
  ## Predict values from model and store in a csv
  
  unique(dat.mod$thermal.guild)
  
  testdata <-expand.grid(year=c('2013','2019'),
                         thermal.guild=c('Temperate','Tropical'),
                         site="PGN1",
                         date="1")%>%
    glimpse()
  
  fits <- predict(mod, newdata=testdata, type = "response",se.fit=T,re.form = NA)
  
  ## Group predictions
  
  predictions = testdata%>%data.frame(fits)%>%
    group_by(year,thermal.guild)%>% #only change here
    dplyr::summarise(response=mean(fit),se.fit=mean(se.fit))%>%
    ungroup()%>%
    glimpse()
  predictions$Realised.function<-i
  
  dataframe<-rbind(dataframe,predictions)
  
}

## Export dataframe 

setwd(data.dir)
dir()

write.csv(dataframe,"Estimates_Test.csv")

## Bring abundance estimates 

setwd(data.dir)
dir()

dataframe<-read.csv("Estimates_Turf.csv")%>%
  glimpse()

## Plot in ggplot

glimpse(dat.sum)
glimpse(dataframe)
unique(dat.sum$Realised.function)

dat.sum$year<-as.factor(dat.sum$year)
dataframe$year<-as.factor(dataframe$year)

## Re-level factors year

unique(dat.sum$thermal.guild)

dat.sum<-dat.sum%>%
  mutate(thermal.guild=fct_relevel(thermal.guild,"Temperate","Tropical"))%>%
  glimpse()

dataframe<-dataframe%>%
  mutate(thermal.guild=fct_relevel(thermal.guild,"Temperate","Tropical"))%>%
  glimpse()

unique(dataframe$Realised.function)

Macroalgae<-ggplot(data=dataframe%>%filter(Realised.function%in%c('Macroalgal removal')),aes(x=year,y=response))+
  ggtitle('(e)')+
  ylab('Turf bites (kg/m2h)')+
  xlab('')+
  geom_jitter(data=dat.sum%>%filter(Realised.function%in%c('Macroalgal removal')),aes(x=year,y=response,colour=thermal.guild),position=position_jitterdodge(jitter.width = 0.2, jitter.height = 0,
                                                                                                                                                                   dodge.width = 0.75, seed = NA), size=2,alpha=0.25,show.legend = F)+
  geom_errorbar(data = dataframe%>%filter(Realised.function%in%c('Macroalgal removal')),aes(x = year,y=response,
                                                                                            ymin=response-se.fit, ymax=response+se.fit,fill=thermal.guild),color = "black",width=.2,
                position=position_jitterdodge(jitter.width = 0, jitter.height = 0,dodge.width = 0.75, seed = NA), size=1)+
  geom_point(data=dataframe%>%filter(Realised.function%in%c('Macroalgal removal')),
             aes(x=year,y=response,fill=thermal.guild),color="black",
             size=4.5,alpha=2,
             position=position_dodge(.75),show.legend = F)+
  geom_point(data = dataframe%>%filter(Realised.function%in%c('Macroalgal removal')), aes(x=year,y=response,colour=thermal.guild),
             size=4,position=position_dodge(.75),show.legend = F)+
  geom_path(data=dataframe%>%filter(Realised.function%in%c('Macroalgal removal')),
            aes(x=year,y=response,group = thermal.guild,color=thermal.guild),
            position=position_dodge(.75),
            size=1,show.legend = F)+
  scale_y_continuous(limits = c(0,400),trans="log1p")+
  scale_colour_manual(labels = c("Temperate","Tropical"),
                      values=c("blue","red"))+
  theme_classic()+
  Theme1+
  theme(legend.position = "bottom",
        legend.direction = "horizontal")
Macroalgae

Turf<-ggplot(data=dataframe%>%filter(Realised.function%in%c('Turf removal')),aes(x=year,y=response))+
  ggtitle('(f)')+
  ylab('')+
  xlab('')+
  geom_jitter(data=dat.sum%>%filter(Realised.function%in%c('Turf + sediment removal','Turf removal')),aes(x=year,y=response,colour=thermal.guild),position=position_jitterdodge(jitter.width = 0.2, jitter.height = 0,
                                                                                                                                                                                       dodge.width = 0.75, seed = NA), size=2,alpha=0.25,show.legend = F)+
  geom_errorbar(data = dataframe%>%filter(Realised.function%in%c('Turf removal')), aes(x = year,y=response,
                                                                                       ymin=response-se.fit, ymax=response+se.fit,fill=thermal.guild),color = "black",width=.2,
                position=position_jitterdodge(jitter.width = 0, jitter.height = 0,dodge.width = 0.75, seed = NA), size=1)+
  geom_point(data=dataframe%>%filter(Realised.function%in%c('Turf removal')),
             aes(x=year,y=response,fill=thermal.guild),color="black",
             size=4.5,alpha=2,
             position=position_dodge(.75),show.legend = F)+
  geom_point(data = dataframe%>%filter(Realised.function%in%c('Turf removal')), aes(x=year,y=response,colour=thermal.guild),
             size=4,position=position_dodge(.75),show.legend = F)+
  geom_path(data=dataframe%>%filter(Realised.function%in%c('Turf removal')),
            aes(x=year,y=response,group = thermal.guild,color=thermal.guild),
            position=position_dodge(.75),
            size=1,show.legend = F)+
  scale_y_continuous(limits=c(0,400),trans="log1p")+
  scale_colour_manual(labels = c("Temperate","Tropical"),
                      values=c("blue","red"))+
  theme_classic()+
  Theme1+
  theme(legend.position = "bottom",
        legend.direction = "horizontal",
        axis.title.y = element_blank())
Turf

Sediments<-ggplot(data=dataframe%>%filter(Realised.function%in%c('Sediment removal')),aes(x=year,y=response))+
  ggtitle('(g)')+
  ylab('')+
  xlab('')+
  geom_jitter(data=dat.sum%>%filter(Realised.function%in%c('Turf + sediment removal')),aes(x=year,y=response,colour=thermal.guild),position=position_jitterdodge(jitter.width = 0.2, jitter.height = 0,
                                                                                                                                                                        dodge.width = 0.75, seed = NA), size=2,alpha=0.25,show.legend = F)+
  geom_errorbar(data = dataframe%>%filter(Realised.function%in%c('Sediment removal')), aes(x = year,y=response,
                                                                                           ymin=response-se.fit, ymax=response+se.fit,fill=thermal.guild),color = "black",width=.2,
                position=position_jitterdodge(jitter.width = 0, jitter.height = 0,dodge.width = 0.75, seed = NA), size=1)+
  geom_point(data=dataframe%>%filter(Realised.function%in%c('Sediment removal')),
             aes(x=year,y=response,fill=thermal.guild),color="black",
             size=4.5,alpha=2,
             position=position_dodge(.75),show.legend = F)+
  geom_point(data = dataframe%>%filter(Realised.function%in%c('Sediment removal')), aes(x=year,y=response,colour=thermal.guild),
             size=4,position=position_dodge(.75),show.legend = F)+
  geom_path(data=dataframe%>%filter(Realised.function%in%c('Sediment removal')),
            aes(x=year,y=response,group = thermal.guild,color=thermal.guild),
            position=position_dodge(.75),
            size=1,show.legend = F)+
  scale_y_continuous(limits=c(0,400),trans="log1p")+
  scale_colour_manual(labels = c("Tropical"),
                      values=c("red"))+
  theme_classic()+
  Theme1+
  theme(legend.position = "bottom",
        legend.direction = "horizontal",
        axis.title.y = element_blank())
Sediments

## Now we can calculate the relative proportion of each phenotype ('how' function) for each time period

glimpse(dat)

prop<-dat%>%
  filter(!herbivore.guild%in%c('Excavator'))%>%
  dplyr::group_by(location,site,year,date,id,herbivore.guild,HG.TG,thermal.guild)%>%
  dplyr::summarise(response=sum(mass.std.bites))%>%
  glimpse()
n_distinct(prop$id)
unique(prop$HG.TG)

dataframe.1=data.frame()

period<-as.vector(unique(prop$year)) ## For n in 3 time periods

for (i in period) {
  
  ## Create frequency of occurrence for each latitudinal bin and frequency class
  
  test<-prop%>%
    filter(year%in%i)%>%
    glimpse()
  unique(test$year)
  n<-sum(test$response)
  
  test.1<-test%>%
    dplyr::group_by(year,HG.TG)%>%
    dplyr::summarise(number=sum(response))%>%
    mutate(percentage=(number/n)*100)%>%
    ungroup()%>%
    glimpse()
  
  dataframe.1<-rbind(dataframe.1,test.1)
  
}

## Plot of relative frequencies of common, frequent and rare species

library(forcats)

unique(dataframe.1$HG.TG)

dataframe.1<-dataframe.1%>%
  mutate(HG.TG=fct_relevel(HG.TG,"Browser_Temperate","Browser_Tropical","Algal farmer_Temperate","Algal farmer_Tropical",
                           "Scraper_Tropical"))%>%
  glimpse()

proportion<-ggplot(dataframe.1,aes(fill=HG.TG, y=percentage, x=year)) + 
  ggtitle('(h)')+
  xlab('')+
  ylab('Contribution (%)')+
  geom_bar(position="stack", stat="identity",colour="black")+
  scale_fill_manual(labels=c("Browser_Temperate","Browser_Tropical","Algal farmer_Temperate","Algal farmer_Tropical",
                             "Scraper_Tropical"),
                    values = c('lightsalmon4','lightsalmon1','green4','green1','darkorange1'))+
  theme_classic()+
  Theme1
proportion

## Arrange plots in a grob

library(ggpubr)

ggarrange(Macroalgae,Turf,Sediments,proportion,
          nrow=1,ncol = 4,
          align = 'hv',
          common.legend = TRUE,
          legend = "none")

## Export plot

setwd(plots.dir)
dir()

name<-'Turf_Herbivore_Change_17_05_2021'

ggsave(paste(name,".png",sep="."),width = 21, height = 7,units = "cm")
ggsave(paste(name,".pdf",sep="."),width = 21, height = 7,units = "cm",dpi=600,useDingbats=FALSE)

## (C) Kelp ----

setwd(data.dir)
dir()

dat<-read.csv('kelp_bites_clean.csv')%>%
  filter(!scientific%in%c('Ophisternon bengalense'))%>%
  filter(!Class%in%c('Elasmobranchii'))%>%
  filter(!observer%in%c('Claude Spencer'))%>%
  glimpse()
n_distinct(dat$scientific) ## 7 species
unique(dat$scientific)
unique(dat$observer)

## Quick check on samples,year, site for unbalanced design

test<-dat%>%
  dplyr::group_by(year,site,date)%>%
  dplyr::summarise(N=n_distinct(id))%>%
  glimpse()

## Bring in traits and filter only herbivores

traits<-read.csv('traits_herbivores.csv')%>%
  filter(!is.na(herbivore.guild))%>%
  filter(!is.na(ThermalMP_5_95))%>%
  dplyr::select(scientific,Realised.function,herbivore.guild,ThermalMP_5_95,MaxLength)%>%
  mutate(thermal.guild=ifelse(ThermalMP_5_95>=23,"Tropical","Temperate"))%>%
  mutate(RF.TG=paste(Realised.function,thermal.guild,sep="_"),
         HG.TG=paste(herbivore.guild,thermal.guild,sep="_"))%>%
  glimpse()

## Filter only herbivores from the data

x<-as.vector(unique(traits$scientific))

dat<-dat%>%
  filter(scientific%in%x)%>%
  glimpse()
n_distinct(dat$scientific) ## 7 species
unique(dat$scientific)

## Supplementary plot on mean +- se densities

dat.summary<-dat%>%
  group_by(year,scientific)%>%
  summarise(mean=mean(mass.std.bites),
            se=se(mass.std.bites))%>%
  glimpse()

## Abundance-rank plot

plot.2013<-ggplot(dat.summary%>%filter(year%in%c('2013')),aes(x=reorder(scientific,-mean),y=mean))+
  geom_col(fill="lightgrey",colour="black")+
  geom_errorbar(aes(ymin=mean, ymax=mean+se),color = "black",width=.2,size=1)+
  ggtitle('(a)')+
  ylab('')+
  xlab('')+
  scale_y_continuous(expand=c(0,0),limits = c(0,25000),trans="log1p")+
  theme_classic()+
  theme(axis.text.x = element_text(angle=45,hjust=1,vjust=1))
plot.2013

plot.2019<-ggplot(dat.summary%>%filter(year%in%c('2019'))%>%filter(scientific%in%c('Kyphosus bigibbus','Kyphosus sydneyanus','Siganus fuscescens','Scarus ghobban')),aes(x=reorder(scientific,-mean),y=mean))+
  geom_col(fill="lightgrey",colour="black")+
  geom_errorbar(aes(ymin=mean, ymax=mean+se),color = "black",width=.2,size=1)+
  ggtitle('(b)')+
  ylab('')+
  xlab('')+
  scale_y_continuous(expand=c(0,0),limits = c(0,25000),trans="log1p")+
  theme_classic()+
  theme(axis.text.x = element_text(angle=45,hjust=1,vjust=1))
plot.2019

## Arrange plots in a grob

library(ggpubr)

ggarrange(plot.2013,plot.2019,
          nrow=1,ncol=2,align = 'hv')

## Export plot

setwd(plots.dir)
dir()

name<-'Appendix_FigureA4_Dominance_Kelp'

ggsave(paste(name,".png",sep="."),width = 21, height = 10,units = "cm")
ggsave(paste(name,".pdf",sep="."),width = 21, height = 10,units = "cm",dpi=600,useDingbats=FALSE)

## Join trait information and summarize data - by trophic guild - thermal guild combination?

dat<-left_join(dat,traits,by="scientific")%>%
  glimpse()
n_distinct(dat$RF.TG) # 5 groups
unique(dat$RF.TG)

dat.sum<-dat%>%
  filter(Realised.function%in%c('Macroalgal removal'))%>%
  dplyr::group_by(location,site,year,date,id,Realised.function,RF.TG,thermal.guild)%>%
  dplyr::summarise(response=sum(mass.std.bites))%>%
  glimpse()
unique(dat.sum$Realised.function)
n_distinct(dat$id)

## Check distribution of densities across years

ggplot(dat.sum,aes(x=reorder(RF.TG,-response),y=response,fill=RF.TG))+
  geom_bar(stat = "identity",colour="black",show.legend = F)+
  facet_wrap(~year,scales = "free")+
  theme_classic()+
  theme(axis.text.x = element_text(angle=90))

ggplot(dat.sum,aes(x=reorder(RF.TG,-response),y=log10(response+1),fill=RF.TG))+
  geom_bar(stat = "identity",colour="black",show.legend = F)+
  facet_wrap(~year,scales = "free")+
  theme_classic()+
  theme(axis.text.x = element_text(angle=90))

## Run GLMMs on HG.TG and save coefficients =- se

## Create a loop to investigate the distribution across HG.TG groups
## All models with a negative binomial distribution

n_distinct(dat.sum$Realised.function) ## 1 groups
unique(dat.sum$Realised.function)

pred.vars.cont=c("Macroalgal removal")

plot.new()
par(mfrow=c(1,3))

for (i in pred.vars.cont) {
  x<-dat.sum%>%filter(Realised.function%in%i)
  hist(x$response,main = paste(i))#Looks best
}

## Run models in a loop

library(performance)
library(sjPlot)

i<-c("Macroalgal removal")

# Create an empty dataframe

dataframe=data.frame()

for (i in pred.vars.cont) {
  
  dat.mod<-dat.sum%>%
    filter(Realised.function%in%i)%>%
    glimpse()
  unique(dat.mod$Realised.function)
  
  ## Recode factors for modelling
  
  dat.mod$year<-as.factor(dat.mod$year)
  dat.mod$thermal.guild<-as.factor(dat.mod$thermal.guild)
  dat.mod$date<-as.factor(dat.mod$date)
  
  ## Investigate sample size distribution
  
  plot(dat.mod$year)
  plot(dat.mod$thermal.guild)
  plot(dat.mod$date)
  
  ## Run model
  
  library(glmmTMB)
  library(MuMIn)
  
  mod<-glmmTMB(response ~ year*thermal.guild + (1|site) + (1|date),
               data = dat.mod,family = tweedie())
  summary(mod)
  car::Anova(mod)
  
  system.time(mod_d1 <- drop1(mod,test="Chisq",all.cols=TRUE))
  print(mod_d1)
  
  ## Calculate R2 = 1 - (Residual Deviance/Null Deviance)
  
  r2_nakagawa(mod)
  
  ## Pairwise comparisons
  
  emmeans(mod, pairwise ~ year:thermal.guild)
  
  ## Check residual pattern
  
  library(DHARMa)
  
  mod.res <- simulateResiduals(mod)
  plot(mod.res)
  
  ## Test zero-inflation
  
  testZeroInflation(mod.res)
  
  ## Plot the effect of gravity by management status
  
  plot_model(mod,type = "eff",terms = c("year","thermal.guild"))
  
  ## Predict values from model and store in a csv
  
  unique(dat.mod$thermal.guild)
  
  testdata <-expand.grid(year=c('2013','2019'),
                         thermal.guild=c('Temperate','Tropical'),
                         site="PGN1",
                         date="1")%>%
    glimpse()
  
  fits <- predict(mod, newdata=testdata, type = "response",se.fit=T,re.form = NA)
  
  ## Group predictions
  
  predictions = testdata%>%data.frame(fits)%>%
    group_by(year,thermal.guild)%>% #only change here
    dplyr::summarise(response=mean(fit),se.fit=mean(se.fit))%>%
    ungroup()%>%
    glimpse()
  predictions$Realised.function<-i
  
  dataframe<-rbind(dataframe,predictions)
  
}

## Export dataframe 

setwd(data.dir)
dir()

write.csv(dataframe,"Estimates_Test.csv")

## Bring abundance estimates 

setwd(data.dir)
dir()

dataframe<-read.csv("Estimates_Kelp.csv")%>%
  glimpse()

## Plot in ggplot

glimpse(dat.sum)
glimpse(dataframe)
unique(dat.sum$Realised.function)

dat.sum$year<-as.factor(dat.sum$year)
dataframe$year<-as.factor(dataframe$year)

## Re-level factors year

unique(dat.sum$thermal.guild)

dat.sum<-dat.sum%>%
  mutate(thermal.guild=fct_relevel(thermal.guild,"Temperate","Tropical"))%>%
  glimpse()

dataframe<-dataframe%>%
  mutate(thermal.guild=fct_relevel(thermal.guild,"Temperate","Tropical"))%>%
  glimpse()

unique(dataframe$Realised.function)

Macroalgae<-ggplot(data=dataframe%>%filter(Realised.function%in%c('Macroalgal removal')),aes(x=year,y=response))+
  ggtitle('(i)')+
  ylab('Kelp bites (kg bites/h)')+
  xlab('')+
  geom_jitter(data=dat.sum%>%filter(Realised.function%in%c('Macroalgal removal')),aes(x=year,y=response,colour=thermal.guild),position=position_jitterdodge(jitter.width = 0.2, jitter.height = 0,
                                                                                                                                                            dodge.width = 0.75, seed = NA), size=2,alpha=0.25,show.legend = F)+
  geom_errorbar(data = dataframe%>%filter(Realised.function%in%c('Macroalgal removal')),aes(x = year,y=response,
                                                                                            ymin=response-se.fit, ymax=response+se.fit,fill=thermal.guild),color = "black",width=.2,
                position=position_jitterdodge(jitter.width = 0, jitter.height = 0,dodge.width = 0.75, seed = NA), size=1)+
  geom_point(data=dataframe%>%filter(Realised.function%in%c('Macroalgal removal')),
             aes(x=year,y=response,fill=thermal.guild),color="black",
             size=4.5,alpha=2,
             position=position_dodge(.75),show.legend = F)+
  geom_point(data = dataframe%>%filter(Realised.function%in%c('Macroalgal removal')), aes(x=year,y=response,colour=thermal.guild),
             size=4,position=position_dodge(.75),show.legend = F)+
  geom_path(data=dataframe%>%filter(Realised.function%in%c('Macroalgal removal')),
            aes(x=year,y=response,group = thermal.guild,color=thermal.guild),
            position=position_dodge(.75),
            size=1,show.legend = F)+
  scale_y_continuous(trans="log1p")+
  scale_colour_manual(labels = c("Temperate","Tropical"),
                      values=c("blue","red"))+
  theme_classic()+
  Theme1+
  theme(legend.position = "bottom",
        legend.direction = "horizontal")
Macroalgae

## Now we can calculate the relative proportion of each phenotype ('how' function) for each time period

glimpse(dat)

prop<-dat%>%
  filter(Realised.function%in%c('Macroalgal removal'))%>%
  dplyr::group_by(location,site,year,date,id,herbivore.guild,HG.TG,thermal.guild)%>%
  dplyr::summarise(response=sum(mass.std.bites))%>%
  glimpse()
n_distinct(prop$id)
unique(prop$HG.TG)

dataframe.1=data.frame()

period<-as.vector(unique(prop$year)) ## For n in 3 time periods

for (i in period) {
  
  ## Create frequency of occurrence for each latitudinal bin and frequency class
  
  test<-prop%>%
    filter(year%in%i)%>%
    glimpse()
  unique(test$year)
  n<-sum(test$response)
  
  test.1<-test%>%
    dplyr::group_by(year,HG.TG)%>%
    dplyr::summarise(number=sum(response))%>%
    mutate(percentage=(number/n)*100)%>%
    ungroup()%>%
    glimpse()
  
  dataframe.1<-rbind(dataframe.1,test.1)
  
}

## Plot of relative frequencies of common, frequent and rare species

library(forcats)

unique(dataframe.1$HG.TG)

dataframe.1<-dataframe.1%>%
  mutate(HG.TG=fct_relevel(HG.TG,"Browser_Temperate","Browser_Tropical"))%>%
  glimpse()

proportion<-ggplot(dataframe.1,aes(fill=HG.TG, y=percentage, x=year)) + 
  ggtitle('(j)')+
  xlab('')+
  ylab('Contribution (%)')+
  geom_bar(position="stack", stat="identity",colour="black")+
  scale_fill_manual(labels=c("Browser_Temperate","Browser_Tropical"),
                    values = c('lightsalmon4','lightsalmon1'))+
  theme_classic()+
  Theme1
proportion

## Arrange plots in a grob

library(ggpubr)

ggarrange(Macroalgae,proportion,
          nrow=1,ncol = 4,
          align = 'hv',
          common.legend = TRUE,
          legend = "none")

## Export plot

setwd(plots.dir)
dir()

name<-'Kelp_Herbivore_Change_17_05_2021'

ggsave(paste(name,".png",sep="."),width = 21, height = 7,units = "cm")
ggsave(paste(name,".pdf",sep="."),width = 21, height = 7,units = "cm",dpi=600,useDingbats=FALSE)

## (D) Kelp Assays ----

setwd(data.dir)
dir()

dat<-read.csv("Kelp_Assays_clean.csv")%>%
  dplyr::filter(!is.na(Weight_in_kelp))%>%
  glimpse()

## Investigate the distribution of the response variable

n_distinct(dat$Stage) ## 2 groups - Adult and Juveniles
unique(dat$Stage)

pred.vars.cont=c("A","J") ## Tweedie

plot.new()
par(mfrow=c(1,2))

for (i in pred.vars.cont) {
  x<-dat%>%filter(Stage%in%i)
  hist(x$Weight_in_kelp,main = paste(i))#Looks best
}

## Run models in a loop

library(performance)
library(sjPlot)

i<-"J"

# Create an empty dataframe

dataframe=data.frame()

for (i in pred.vars.cont) {
  
  dat.sum<-dat%>%
    filter(Stage%in%i)%>%
    glimpse()
  unique(dat.sum$Stage)
  
  ## Recode factors for modelling
  
  dat.sum$Year<-as.factor(dat.sum$Year)
  dat.sum$Stage<-as.factor(dat.sum$Stage)
  dat.sum$Period<-as.factor(dat.sum$Period)
  
  ## Investigate sample size distribution
  
  plot(dat.sum$Year)
  plot(dat.sum$Period)
  
  ## Run model
  
  library(glmmTMB)
  library(MuMIn)
  
  mod<-glmmTMB(Weight_in_kelp ~ Year*Period + (1|Site) + (1|Period),
               data = dat.sum,family = tweedie())
  summary(mod)
  car::Anova(mod)
  
  system.time(mod_d1 <- drop1(mod,test="Chisq",all.cols=TRUE))
  print(mod_d1)
  
  ## Pairwise comparison
  
  emmeans(mod, pairwise ~ Year:Period)
  
  ## Calculate R2 = 1 - (Residual Deviance/Null Deviance)
  
  r2_nakagawa(mod)
  
  ## Check residual pattern
  
  library(DHARMa)
  
  mod.res <- simulateResiduals(mod)
  plot(mod.res)
  
  ## Test zero-inflation
  
  testZeroInflation(mod.res)
  
  ## Plot the effect of gravity by management status
  
  plot_model(mod,type = "eff",terms = c("Year","Period"))
  
  ## Predict values from model and store in a csv
  
  unique(dat$Year)
  unique(dat$Period)
  
  testdata <-expand.grid(Year=c('2013','2019'),
                         Period=c('0','24','48'),
                         Site="PGN1")%>%
    glimpse()
  
  fits <- predict(mod, newdata=testdata, type = "response",se.fit=T,re.form = NA)
  
  ## Group predictions
  
  predictions = testdata%>%data.frame(fits)%>%
    group_by(Year,Period)%>% #only change here
    dplyr::summarise(response=mean(fit),se.fit=mean(se.fit))%>%
    ungroup()%>%
    glimpse()
  predictions$Stage<-i
  
  dataframe<-rbind(dataframe,predictions)
  
}

## Export dataframe 

setwd(data.dir)
dir()

write.csv(dataframe,"Estimates_Essays.csv")

## Bring abundance estimates 

setwd(data.dir)
dir()

dataframe<-read.csv("Estimates_Essays.csv")%>%
  glimpse()

## Plot in ggplot

glimpse(dat)
glimpse(dataframe)

dat$Year<-as.factor(dat$Year)
dataframe$Year<-as.factor(dataframe$Year)

dat$Period<-as.factor(dat$Period)
dataframe$Period<-as.factor(dataframe$Period)

## Re-level factors year

unique(dat$Year)

dat<-dat%>%
  mutate(Year=fct_relevel(Year,"2013","2019"))%>%
  glimpse()

dataframe<-dataframe%>%
  mutate(Year=fct_relevel(Year,"2013","2019"))%>%
  glimpse()

Adult.plot<-ggplot(data=dataframe%>%filter(Stage%in%c('A')),aes(x=Period,y=response))+
  ggtitle('(i)')+
  ylab('Individual biomass (g) remaining')+
  xlab('')+
  geom_jitter(data=dat%>%filter(Stage%in%c('A')),aes(x=Period,y=Weight_in_kelp,colour=Year),position=position_jitterdodge(jitter.width = 0.2, jitter.height = 0,
                                                                                                                                               dodge.width = 0.75, seed = NA), size=2,alpha=0.25,show.legend = F)+
  geom_errorbar(data = dataframe%>%filter(Stage%in%c('A')), aes(x = Period,y=response,
                                                                                ymin=response-se.fit, ymax=response+se.fit,fill=Year),color = "black",width=.2,
                position=position_jitterdodge(jitter.width = 0, jitter.height = 0,dodge.width = 0.75, seed = NA), size=1)+
  geom_point(data=dataframe%>%filter(Stage%in%c('A')),
             aes(x=Period,y=response,fill=Year),color="black",
             size=4.5,alpha=2,
             position=position_dodge(.75),show.legend = F)+
  geom_point(data = dataframe%>%filter(Stage%in%c('A')), aes(x=Period,y=response,colour=Year),
             size=4,position=position_dodge(.75),show.legend = F)+
  geom_path(data=dataframe%>%filter(Stage%in%c('A')),
            aes(x=Period,y=response,group = Year,color=Year),
            position=position_dodge(.75),
            size=1,show.legend = F)+
  #scale_y_continuous(trans="log1p",breaks = c(1,10,100,500,1000))+
  scale_colour_manual(labels = c("2013","2019"),
                      values=c("orange","purple"))+
  theme_classic()+
  Theme1+
  theme(legend.position = "bottom",
        legend.direction = "horizontal")
Adult.plot

Juvenile.plot<-ggplot(data=dataframe%>%filter(Stage%in%c('J')),aes(x=Period,y=response))+
  ggtitle('(j)')+
  ylab('')+
  xlab('')+
  geom_jitter(data=dat%>%filter(Stage%in%c('J')),aes(x=Period,y=Weight_in_kelp,colour=Year),position=position_jitterdodge(jitter.width = 0.2, jitter.height = 0,
                                                                                                                          dodge.width = 0.75, seed = NA), size=2,alpha=0.25,show.legend = F)+
  geom_errorbar(data = dataframe%>%filter(Stage%in%c('J')), aes(x = Period,y=response,
                                                                ymin=response-se.fit, ymax=response+se.fit,fill=Year),color = "black",width=.2,
                position=position_jitterdodge(jitter.width = 0, jitter.height = 0,dodge.width = 0.75, seed = NA), size=1)+
  geom_point(data=dataframe%>%filter(Stage%in%c('J')),
             aes(x=Period,y=response,fill=Year),color="black",
             size=4.5,alpha=2,
             position=position_dodge(.75),show.legend = T)+
  geom_point(data = dataframe%>%filter(Stage%in%c('J')), aes(x=Period,y=response,colour=Year),
             size=4,position=position_dodge(.75),show.legend = T)+
  geom_path(data=dataframe%>%filter(Stage%in%c('J')),
            aes(x=Period,y=response,group = Year,color=Year),
            position=position_dodge(.75),
            size=1,show.legend = F)+
  #scale_y_continuous(trans="log1p",breaks = c(1,10,100,500,1000))+
  scale_colour_manual(labels = c("2013","2019"),
                      values=c("orange","purple"))+
  theme_classic()+
  Theme1+
  theme(legend.position = "bottom",
        legend.direction = "horizontal")
Juvenile.plot

## Arrange plots in a grob

library(ggpubr)

ggarrange(Adult.plot,Juvenile.plot,
          nrow=1,ncol = 2,
          align = 'hv',
          common.legend = TRUE,
          legend = "bottom")

## Export plot

setwd(plots.dir)
dir()

name<-'Kelp_essays_Change'

ggsave(paste(name,".png",sep="."),width = 15, height = 10,units = "cm",dpi=600)
ggsave(paste(name,".pdf",sep="."),width = 15, height = 10,units = "cm",dpi=600,useDingbats=FALSE)

##########################################################################################################
###################################### (4) Supplementary tables ##########################################

## (A) Species-specific densities ----

setwd(data.dir)
dir()

dat<-read.csv('dovs_clean.csv')%>%
  filter(!scientific%in%c('Ophisternon bengalense'))%>%
  filter(!Class%in%c('Elasmobranchii'))%>%
  filter(year!=2017)%>%
  filter(!Site%in%c('PGS1'))%>%
  glimpse()
unique(dat$year)
unique(dat$Site)
n_distinct(dat$scientific) ## 65 species
unique(dat$scientific)

## Bring in traits and filter only herbivores

traits<-read.csv('traits_final.csv')%>%
  filter(!is.na(herbivore.guild))%>%
  filter(!is.na(ThermalMP_5_95))%>%
  dplyr::select(scientific,herbivore.guild,ThermalMP_5_95,MaxLength)%>%
  mutate(thermal.guild=ifelse(ThermalMP_5_95>=23,"Tropical","Temperate"))%>%
  mutate(HG.TG=paste(herbivore.guild,thermal.guild,sep="_"))%>%
  glimpse()

## Filter only herbivores from the data

x<-as.vector(unique(traits$scientific))

dat<-dat%>%
  filter(scientific%in%x)%>%
  glimpse()
n_distinct(dat$scientific) ## 10 species
unique(dat$scientific)

## Summary table per year and site

summary<-dat%>%
  group_by(year,Site,scientific)%>%
  summarise(mean=mean(number,na.rm=T),
            se=se(number))%>%
  glimpse()

## Round mean and se variables

summary$mean<-round(summary$mean,1)
summary$se<-round(summary$se,1)

## Create a joined column

summary<-summary%>%
  mutate(variable=paste(mean,se,sep = " - "))%>%
  glimpse()

## Convert to Wide format 

library(tidyr)

sum.2006<-summary%>%
  filter(year==2006)%>%
  dplyr::select(scientific,Site,variable)%>%
  pivot_wider(names_from = "Site",values_from = variable) %>%
  mutate_all(~replace_na(., 0))%>%
  glimpse()

sum.2013<-summary%>%
  filter(year==2013)%>%
  dplyr::select(scientific,Site,variable)%>%
  pivot_wider(names_from = "Site",values_from = variable) %>%
  mutate_all(~replace_na(., 0))%>%
  glimpse()

sum.2019<-summary%>%
  filter(year==2019)%>%
  dplyr::select(scientific,Site,variable)%>%
  pivot_wider(names_from = "Site",values_from = variable) %>%
  mutate_all(~replace_na(., 0))%>%
  glimpse()

## Join datasets

sum.final<-cbind(sum.2006,sum.2013,sum.2019)%>%
  glimpse()

## Export dataset

setwd(data.dir)
dir()

write.csv(sum.final,"Densities.Summary.csv")

## (B) Species-specific abundances (MaxN) ----

setwd(data.dir)
dir()

dat<-read.csv('MaxN_clean.csv')%>%
  filter(!scientific%in%c('Ophisternon bengalense'))%>%
  filter(!Class%in%c('Elasmobranchii'))%>%
  filter(year!=2017)%>%
  glimpse()
n_distinct(dat$scientific) ## 69 species
unique(dat$scientific)

## Bring in traits and filter only herbivores

traits<-read.csv('traits_final.csv')%>%
  filter(!is.na(herbivore.guild))%>%
  filter(!is.na(ThermalMP_5_95))%>%
  dplyr::select(scientific,herbivore.guild,ThermalMP_5_95,MaxLength)%>%
  mutate(thermal.guild=ifelse(ThermalMP_5_95>=23,"Tropical","Temperate"))%>%
  mutate(HG.TG=paste(herbivore.guild,thermal.guild,sep="_"))%>%
  glimpse()

## Filter only herbivores from the data

x<-as.vector(unique(traits$scientific))

dat<-dat%>%
  filter(scientific%in%x)%>%
  glimpse()
n_distinct(dat$scientific) ## 18 species
unique(dat$scientific)

## Summary table per year and site

summary<-dat%>%
  group_by(year,site,scientific)%>%
  summarise(mean=mean(number,na.rm=T),
            se=se(number))%>%
  glimpse()

## Round mean and se variables

summary$mean<-round(summary$mean,1)
summary$se<-round(summary$se,1)

## Create a joined column

summary<-summary%>%
  mutate(variable=paste(mean,se,sep = " - "))%>%
  glimpse()

## Convert to Wide format 

sum.2013<-summary%>%
  filter(year==2013)%>%
  dplyr::select(scientific,site,variable)%>%
  pivot_wider(names_from = "site",values_from = variable) %>%
  mutate_all(~replace_na(., 0))%>%
  glimpse()

sum.2019<-summary%>%
  filter(year==2019)%>%
  dplyr::select(scientific,site,variable)%>%
  pivot_wider(names_from = "site",values_from = variable) %>%
  mutate_all(~replace_na(., 0))%>%
  glimpse()

## Join datasets

sum.final<-cbind(sum.2013,sum.2019)%>%
  glimpse()

## Export dataset

setwd(data.dir)
dir()

write.csv(sum.final,"Abundances.Summary.csv")

## (C) Species-specific values in turf bite rates ----

setwd(data.dir)
dir()

dat<-read.csv('turf_bites_clean.csv')%>%
  filter(!scientific%in%c('Ophisternon bengalense'))%>%
  filter(!Class%in%c('Elasmobranchii'))%>%
  glimpse()
n_distinct(dat$scientific) ## 10 species
unique(dat$scientific)

## Bring in traits and filter only herbivores

traits<-read.csv('traits_final.csv')%>%
  filter(!is.na(herbivore.guild))%>%
  filter(!is.na(ThermalMP_5_95))%>%
  dplyr::select(scientific,herbivore.guild,ThermalMP_5_95,MaxLength)%>%
  mutate(thermal.guild=ifelse(ThermalMP_5_95>=23,"Tropical","Temperate"))%>%
  mutate(HG.TG=paste(herbivore.guild,thermal.guild,sep="_"))%>%
  glimpse()

## Filter only herbivores from the data

x<-as.vector(unique(traits$scientific))

dat<-dat%>%
  filter(scientific%in%x)%>%
  glimpse()
n_distinct(dat$scientific) ## 10 species
unique(dat$scientific)

## Summary table per year and site

summary<-dat%>%
  group_by(year,site,scientific)%>%
  summarise(mean=mean(mass.std.bites,na.rm=T),
            se=se(mass.std.bites))%>%
  glimpse()

## Round mean and se variables

summary$mean<-round(summary$mean,1)
summary$se<-round(summary$se,1)

## Create a joined column

summary<-summary%>%
  mutate(variable=paste(mean,se,sep = " - "))%>%
  glimpse()

## Convert to Wide format 

sum.2013<-summary%>%
  filter(year==2013)%>%
  dplyr::select(scientific,site,variable)%>%
  pivot_wider(names_from = "site",values_from = variable) %>%
  mutate_all(~replace_na(., 0))%>%
  glimpse()

sum.2019<-summary%>%
  filter(year==2019)%>%
  dplyr::select(scientific,site,variable)%>%
  pivot_wider(names_from = "site",values_from = variable) %>%
  mutate_all(~replace_na(., 0))%>%
  glimpse()

## Join datasets

sum.final<-cbind(sum.2013,sum.2019)%>%
  glimpse()

## Export dataset

setwd(data.dir)
dir()

write.csv(sum.final,"Turf.Bites.Summary.csv")

## (D) Species-specific values in Kelp bite rates ----

setwd(data.dir)
dir()

dat<-read.csv('kelp_bites_clean.csv')%>%
  filter(!scientific%in%c('Ophisternon bengalense'))%>%
  filter(!Class%in%c('Elasmobranchii'))%>%
  filter(!observer%in%c('Claude Spencer'))%>%
  filter(scientific%in%c('Kyphosus bigibbus','Kyphosus cornelli','Kyphosus sydneyanus','Scarus ghobban',
                         'Siganus fuscescens'))%>%
  glimpse()
n_distinct(dat$scientific) ## 10 species
unique(dat$scientific)
unique(dat$observer)


## Summary table per year and site

summary<-dat%>%
  group_by(year,site,scientific)%>%
  summarise(mean=mean(mass.std.bites,na.rm=T),
            se=se(mass.std.bites))%>%
  glimpse()

## Round mean and se variables

summary$mean<-round(summary$mean,1)
summary$se<-round(summary$se,1)

## Create a joined column

summary<-summary%>%
  mutate(variable=paste(mean,se,sep = " - "))%>%
  glimpse()

## Convert to Wide format 

sum.2013<-summary%>%
  filter(year==2013)%>%
  dplyr::select(scientific,site,variable)%>%
  pivot_wider(names_from = "site",values_from = variable) %>%
  mutate_all(~replace_na(., 0))%>%
  glimpse()

sum.2019<-summary%>%
  filter(year==2019)%>%
  dplyr::select(scientific,site,variable)%>%
  pivot_wider(names_from = "site",values_from = variable) %>%
  mutate_all(~replace_na(., 0))%>%
  glimpse()

## Join datasets

sum.final<-cbind(sum.2013,sum.2019)%>%
  glimpse()

## Export dataset

setwd(data.dir)
dir()

write.csv(sum.final,"Kelp.Bites.Summary.csv")

#############################################################################################################
######################################### (5) Turf Models controlling for abundance/MaxN ####################

## Bring in turf bite data

setwd(data.dir)
dir()

dat<-read.csv('turf_bites_clean.csv')%>%
  filter(!scientific%in%c('Ophisternon bengalense'))%>%
  filter(!Class%in%c('Elasmobranchii'))%>%
  glimpse()
n_distinct(dat$scientific) ## 10 species
unique(dat$scientific)

## Quick check on samples,year, site for unbalanced design

test<-dat%>%
  dplyr::group_by(year,site,date)%>%
  dplyr::summarise(N=n_distinct(id))%>%
  glimpse()
n_distinct(dat$id)

## Bring in traits and filter only herbivores

traits<-read.csv('traits_final.csv')%>%
  filter(!is.na(herbivore.guild))%>%
  filter(!is.na(ThermalMP_5_95))%>%
  dplyr::select(scientific,herbivore.guild,ThermalMP_5_95,MaxLength)%>%
  mutate(thermal.guild=ifelse(ThermalMP_5_95>=23,"Tropical","Temperate"))%>%
  mutate(HG.TG=paste(herbivore.guild,thermal.guild,sep="_"))%>%
  filter(herbivore.guild%in%c('Browser'))%>%
  glimpse()

## Filter only herbivores from the data

x<-as.vector(unique(traits$scientific))

dat<-dat%>%
  filter(scientific%in%x)%>%
  glimpse()
n_distinct(dat$scientific) ## 4 species
unique(dat$scientific)

## Join trait information and summarize data - by trophic guild - thermal guild combination?

dat<-left_join(dat,traits,by="scientific")%>%
  glimpse()
unique(dat$herbivore.guild)
unique(dat$thermal.guild) # 2 groups

dat.sum<-dat%>%
  dplyr::group_by(location,site,year,date,id,thermal.guild)%>%
  dplyr::summarise(response=sum(mass.std.bites))%>%
  glimpse()

## Bring in MaxN data for each sample

library(tidyr)

setwd(data.dir)
dir()

maxn<-read.csv("Turf.Bite.Analyses.csv")%>%
  dplyr::select(id,browser.temperate.MaxN,browser.tropical.MaxN)%>%
  glimpse()

unique(dat.sum$thermal.guild)

## Change column names

colnames(maxn)<-c('id','Temperate','Tropical')

maxn<-maxn%>%
  pivot_longer(cols = c(2:3), names_to = "thermal.guild", values_to = "MaxN")%>%
  glimpse()

## Join with data analysis frame

dat.sum<-left_join(dat.sum,maxn,by=c("id","thermal.guild"))%>%
  glimpse()
dat.sum$year<-as.factor(dat.sum$year)

## Run GLMMs on HG.TG and save coefficients =- se

## Recode factors for modelling
  
dat.sum$year<-as.factor(dat.sum$year)
dat.sum$thermal.guild<-as.factor(dat.sum$thermal.guild)
dat.sum$date<-as.factor(dat.sum$date)
  
## Investigate sample size distribution
  
plot(dat.sum$year)
plot(dat.sum$thermal.guild)
plot(dat.sum$date)
  
## Run model
  
library(glmmTMB)
library(MuMIn)
  
mod<-glmmTMB(response ~ year*thermal.guild + MaxN*thermal.guild + (1|site) + (1|date),
               data = dat.sum,family = tweedie())
summary(mod)
car::Anova(mod)
  
system.time(mod_d1 <- drop1(mod,test="Chisq",all.cols=TRUE))
print(mod_d1)
  
#############################################################################################################
######################################### (6) Kelp Models controlling for abundance/MaxN ###################

## Bring in turf bite data

setwd(data.dir)
dir()

dat<-read.csv('kelp_bites_clean.csv')%>%
  filter(!scientific%in%c('Ophisternon bengalense'))%>%
  filter(!Class%in%c('Elasmobranchii'))%>%
  filter(!observer%in%c('Claude Spencer'))%>%
  glimpse()
n_distinct(dat$scientific) ## 21 species
unique(dat$scientific)
unique(dat$observer)

## Quick check on samples,year, site for unbalanced design

test<-dat%>%
  dplyr::group_by(year,site,date)%>%
  dplyr::summarise(N=n_distinct(id))%>%
  glimpse()
n_distinct(dat$id)

## Bring in traits and filter only herbivores

traits<-read.csv('traits_final.csv')%>%
  filter(!is.na(herbivore.guild))%>%
  filter(!is.na(ThermalMP_5_95))%>%
  dplyr::select(scientific,herbivore.guild,ThermalMP_5_95,MaxLength)%>%
  mutate(thermal.guild=ifelse(ThermalMP_5_95>=23,"Tropical","Temperate"))%>%
  mutate(HG.TG=paste(herbivore.guild,thermal.guild,sep="_"))%>%
  filter(herbivore.guild%in%c('Browser'))%>%
  glimpse()

## Filter only herbivores from the data

x<-as.vector(unique(traits$scientific))

dat<-dat%>%
  filter(scientific%in%x)%>%
  glimpse()
n_distinct(dat$scientific) ## 4 species
unique(dat$scientific)

## Join trait information and summarize data - by trophic guild - thermal guild combination?

dat<-left_join(dat,traits,by="scientific")%>%
  glimpse()
unique(dat$herbivore.guild)
unique(dat$thermal.guild) # 2 groups

dat.sum<-dat%>%
  dplyr::group_by(location,site,year,date,id,thermal.guild)%>%
  dplyr::summarise(response=sum(mass.std.bites))%>%
  glimpse()

## Bring in MaxN data for each sample

library(tidyr)

setwd(data.dir)
dir()

maxn<-read.csv("Kelp.Bite.Analyses.csv")%>%
  dplyr::select(id,browser.temperate.MaxN,browser.tropical.MaxN)%>%
  glimpse()
unique(dat.sum$thermal.guild)

## Change column names

colnames(maxn)<-c('id','Temperate','Tropical')

maxn<-maxn%>%
  pivot_longer(cols = c(2:3), names_to = "thermal.guild", values_to = "MaxN")%>%
  glimpse()

## Join with data analysis frame

dat.sum<-left_join(dat.sum,maxn,by=c("id","thermal.guild"))%>%
  glimpse()
dat.sum$year<-as.factor(dat.sum$year)

## Run GLMMs on HG.TG and save coefficients =- se

## Recode factors for modelling

dat.sum$year<-as.factor(dat.sum$year)
dat.sum$thermal.guild<-as.factor(dat.sum$thermal.guild)
dat.sum$date<-as.factor(dat.sum$date)

## Investigate sample size distribution

plot(dat.sum$year)
plot(dat.sum$thermal.guild)
plot(dat.sum$date)

## Run model

library(glmmTMB)
library(MuMIn)

mod<-glmmTMB(response ~ year*thermal.guild + MaxN*thermal.guild + (1|site) + (1|date),
             data = dat.sum,family = tweedie())

## Some warning messages 

library(DHARMa)

mod.res <- simulateResiduals(mod)
plot(mod.res)

summary(mod)
car::Anova(mod)

system.time(mod_d1 <- drop1(mod,test="Chisq",all.cols=TRUE))
print(mod_d1)

############################################################################################
################ (7) Testing for temporal variability ######################################

## (A) Trait spaces ----

# Bring trait data

setwd(data.dir)
dir()

traits<-read.csv('traits_herbivores.csv')%>%
  filter(!is.na(ThermalMP_5_95))%>%
  #dplyr::select(scientific,herbivore.guild,ThermalMP_5_95,MaxLength)%>%
  filter(!is.na(herbivore.guild))%>% ## 27 nominal herbivores
  glimpse()

## Load density data

x<-as.vector(unique(traits$scientific))

setwd(data.dir)
dir()

dat<-read.csv('sDOV_seasonal.csv')%>%
  filter(!scientific%in%c('Ophisternon bengalense'))%>%
  filter(!scientific%in%c('Cirripectes hutchinsi'))%>%
  #filter(year!=2017)%>% Multi-year data
  filter(!Site%in%c('PGS1'))%>%
  filter(scientific%in%x)%>%
  rename(id=SAMPLE_ID)%>%
  glimpse()
unique(dat$Site)
unique(dat$year)
n_distinct(dat$scientific) ## Of the 27 - only 12 species observed in s-DOVs - but note some species were pooled for analyses
unique(dat$scientific)

## Filter out species not present from trait data

x<-as.vector(unique(dat$scientific))

traits<-traits%>%
  filter(scientific%in%x)%>%
  glimpse()

## Change nomenclature to match matrix

traits$scientific<-gsub(' ','.',traits$scientific)

### Check for missing trait information

sum(is.na(traits))/prod(dim(traits))*100
apply(traits,2,function(col)sum(is.na(col))/length(col))*100

## Convert first column to row names

traits<-data.frame(traits)
rownames(traits) <- traits[,1] #Assigning row names from 1st column 
traits[,1] <- NULL #Removing the first column

## See sampling effort per year and site

test<-dat%>%
  dplyr::group_by(year,month,Site)%>%
  dplyr::summarise(N=n_distinct(id))%>%
  glimpse()

## (A.1) Inter-year variability ----
# Only data for March/April for 2017

dat.year<-dat%>%
  filter(!Site%in%c("PGN1"))%>%
  filter(month%in%c("3"))%>%
  glimpse()
unique(dat.year$month)
unique(dat.year$year)

## Resample data to select 8 transects per year

dat.year$year<-as.factor(dat.year$year)
levels(dat.year$year)

## Firt we convert to long format for each year separately 

library(tidyr)

dataframe=data.frame()

number<-1:99
time.period<-as.vector(unique(dat.year$year))
site<-as.vector(unique(dat.year$Site))

for(j in number) {
  
  for (i in time.period) {
    
    for(k in site) {
      
      Abun<-dat%>%
        filter(year%in%i)%>%
        filter(Site%in%k)%>%
        dplyr::group_by(id,scientific)%>%
        summarise(number=sum(number))%>%
        ungroup()%>%
        pivot_wider(names_from = "scientific",values_from = number) %>%
        mutate_all(~replace_na(., 0))%>%
        sample_n(8,replace = T)%>%
        glimpse()
      
      ## Reconvert to long format
      
      Abun<-Abun%>%
        pivot_longer(cols = c(-1), names_to = "scientific", values_to = "number")%>%
        glimpse()
      
      ## Create metadata 
      
      metadata<-dat%>%
        dplyr::select(Location,Site,year,id)%>%
        glimpse()
      
      new.dat<-left_join(Abun,metadata,by="id")%>%
        distinct(id,scientific,.keep_all = TRUE)%>%
        glimpse()
      
      dataframe<-rbind(dataframe,new.dat)
      
    }
  }
}

## Group data and obtain mean 

test<-dataframe%>%
  group_by(year,Site,scientific)%>%
  summarise(number=mean(number))%>%
  glimpse()

## Recheck sampling effort

dat.year<-test

## Clean working directory

rm(Abun,metadata,new.dat,dataframe)

## Summarize data by year

test<-dat.year%>%
  dplyr::group_by(year,scientific)%>%
  dplyr::summarise(number=sum(number))%>%
  glimpse()

## Check the shape of abundance distributions after transformation of species abundances

plot.2017<-ggplot(test%>%filter(year%in%c('2017'))%>%filter(number>0),aes(x=reorder(scientific,-number),y=number,fill=scientific))+
  geom_bar(stat = "identity",colour="black",show.legend = F)+
  scale_y_continuous(expand = c(0,0),limits = c(0,45))+
  theme_classic()+
  theme(axis.text.x = element_text(angle=90))

plot.2019<-ggplot(test%>%filter(year%in%c('2019'))%>%filter(number>0),aes(x=reorder(scientific,-number),y=number,fill=scientific))+
  geom_bar(stat = "identity",colour="black",show.legend = F)+
  scale_y_continuous(expand = c(0,0),limits = c(0,45))+
  theme_classic()+
  theme(axis.text.x = element_text(angle=90))

## Arrange plots in a grob

library(ggpubr)

ggarrange(plot.2017,plot.2019,nrow=1,ncol=3,align = 'hv')

## Check normal histograms

plot.new()
par(mfrow=c(1,2))

hist(test$number)
hist(log10(test$number+1))

## Create a species x abundances matrix for each year

library(tidyr)
library(dplyr)

abundances<-test%>%
  dplyr::select(year,scientific,number)%>%
  pivot_wider(names_from = "scientific",values_from = number) %>%
  mutate_all(~replace_na(., 0))%>%
  glimpse()

## Convert species to column names

abundances<-data.frame(abundances)
rownames(abundances) <- abundances[,1] #Assigning row names from 1st column 
abundances[,1] <- NULL #Removing the first column

## Transform species abundances

abundances <- log10(abundances+1)

#Order species in the community matrix in same order as traits

x<-as.vector(row.names(traits))

abundances <- abundances[,x]

## Check identical species are the same in trait and species matrix

identical(rownames(traits),colnames(abundances))

#######################################
## PLOT SPECIES DENSITY BEFORE/AFTER ##
#######################################

# REPEAT EACH SPECIES ROW AS MANY TIMES AS ITS ABUNDANCE (ON LOG SCALE - MULTIPLIED BY 10 AND ROUNDED)
# TO INCORPORATE ABUNDANCE INTO DENSITY ESTIMATE

n.times <- round(as.numeric(abundances[1,])*10,0)
density_2017 <- traits[rep(seq_len(nrow(traits)), n.times),]

n.times <- round(as.numeric(abundances[2,])*10,0)
density_2019 <- traits[rep(seq_len(nrow(traits)), n.times),]

## Calculate the community centroid (weighthed mean trait values)

## 2017

names(traits)

trait_space <- traits[,c("ThermalMP_5_95","mode.2019")]

centroids.2017 <- functcomp(trait_space, as.matrix(abundances[1,]))

## 2013

names(traits)

trait_space <- traits[,c("ThermalMP_5_95","mode.2019")]

centroids.2019<- functcomp(trait_space, as.matrix(abundances[2,]))

## MARGINAL DENSITY USING COWPLOT ##

## 2017

range(traits$ThermalMP_5_95)
range(traits$MaxLength)

## Relabel factor variables

unique(density_2017$herbivore.guild)

density_2017<-density_2017%>%
  mutate(herbivore.guild=fct_relevel(herbivore.guild,"Browser","Algal farmer","Scraper","Sediment sucker"))%>%
  glimpse()

## Estimate 50%, 75% AND 95% Highest Density Estimates 

getLevel <- function(x,y,prob=c(0.5,0.75,0.95)) {
  kk <- MASS::kde2d(x,y)
  dx <- diff(kk$x[1:2])
  dy <- diff(kk$y[1:2])
  sz <- sort(kk$z)
  c1 <- cumsum(sz) * dx * dy
  approx(c1, sz, xout = 1 - prob)$y
}
L95 <- getLevel(density_2017$ThermalMP_5_95,density_2017$mode.2019)

## Plot 2D density estimates with marginal densities for each trait axis
names(centroids.2017)
d.2017 <- ggplot(density_2017, aes(ThermalMP_5_95,mode.2019)) + 
  ggtitle('(a)')+
  ylab('Body size (cm)')+
  xlab('Thermal affinity (°C)')+
  stat_density_2d(geom = "polygon", aes(alpha = ..level..),show.legend = F,colour="black",breaks=L95) +
  geom_point(data=density_2017, aes(ThermalMP_5_95,mode.2019,
                                    fill=herbivore.guild),colour="black",size=1.5,show.legend = F) +
  geom_point(data=density_2017, aes(ThermalMP_5_95,mode.2019,
                                    color=herbivore.guild),size=1,show.legend = F) +
  geom_point(data=centroids.2017,aes(x=ThermalMP_5_95,y=mode.2019),shape=23,size=3,fill="orange")+
  geom_vline(xintercept = 23,linetype="dashed",size=1)+
  scale_colour_manual(labels = c("Browser","Algal farmer","Scraper","Sediment sucker"),values=c("yellow4","tomato1","steelblue1","orchid2"))+
  xlim(18.0405,30) + ylim(0,75)+
  theme_classic()+
  Theme1

xdens <- axis_canvas(d.2017, axis="x") +
  geom_density(data=density_2017, aes(x=ThermalMP_5_95), fill="black",alpha=0.2)

ydens <- axis_canvas(d.2017, axis="y", coord_flip = TRUE) +
  geom_density(data=density_2017, aes(x=mode.2019), fill="black",alpha=0.2) +
  coord_flip()

p1 <- insert_xaxis_grob(d.2017, xdens,position = "top")
p2 <- insert_yaxis_grob(p1, ydens,position = "right")
ggdraw(p2)

## 2019

range(traits$ThermalMP_5_95)
range(traits$MaxLength)

## Relabel factor variables

unique(density_2019$herbivore.guild)

density_2019<-density_2019%>%
  mutate(herbivore.guild=fct_relevel(herbivore.guild,"Sediment sucker","Browser","Algal farmer","Scraper"))%>%
  glimpse()

## Estimate 50%, 75% AND 95% Highest Density Estimates 

getLevel <- function(x,y,prob=c(0.5,0.75,0.95)) {
  kk <- MASS::kde2d(x,y)
  dx <- diff(kk$x[1:2])
  dy <- diff(kk$y[1:2])
  sz <- sort(kk$z)
  c1 <- cumsum(sz) * dx * dy
  approx(c1, sz, xout = 1 - prob)$y
}
L95 <- getLevel(density_2019$ThermalMP_5_95,density_2019$mode.2019)

## Plot 2D density estimates with marginal densities for each trait axis

d.2019 <- ggplot(density_2019, aes(ThermalMP_5_95,mode.2019)) + 
  ggtitle('(b)')+
  ylab('')+
  xlab('')+
  stat_density_2d(geom = "polygon", aes(alpha = ..level..),show.legend = F,colour="black",breaks=L95) +
  geom_point(data=density_2019, aes(ThermalMP_5_95,mode.2019,
                                    fill=herbivore.guild),colour="black",size=1.5,show.legend = F) +
  geom_point(data=density_2019, aes(ThermalMP_5_95,mode.2019,
                                    color=herbivore.guild),size=1,show.legend = F) +
  geom_point(data=centroids.2019,aes(x=ThermalMP_5_95,y=mode.2019),shape=23,size=3,fill="orange")+
  geom_vline(xintercept = 23,linetype="dashed",size=1)+
  scale_colour_manual(labels = c("Sediment sucker","Browser","Algal farmer","Scraper"),
                      values=c("orchid2","yellow4","tomato1","steelblue1"))+
  xlim(18.0405,30) + ylim(0,75)+
  theme_classic()+
  Theme1

xdens <- axis_canvas(d.2019, axis="x") +
  geom_density(data=density_2019, aes(x=ThermalMP_5_95), fill="black",alpha=0.2)

ydens <- axis_canvas(d.2019, axis="y", coord_flip = TRUE) +
  geom_density(data=density_2019, aes(x=mode.2019), fill="black",alpha=0.2) +
  coord_flip()

p3 <- insert_xaxis_grob(d.2019, xdens,position = "top")
p4 <- insert_yaxis_grob(p3, ydens,position = "right")
ggdraw(p4)

## Arrange plots in a grob

library(ggpubr)

ggarrange(p2,p4,ncol=2,nrow = 1)

## Export plot

setwd(plots.dir)
dir()

name<-"Inter_annual_Trait_Space"

ggsave(paste(name,".png",sep="."),width = 19, height = 8,units = "cm",dpi=600)
ggsave(paste(name,".pdf",sep="."),width = 19, height = 8,units = "cm",dpi=600,useDingbats=FALSE)

## (A.2) Intra-year variability ----

# Only data for 2019 - March and April

dat$month<-as.factor(dat$month)

dat.month<-dat%>%
  filter(!Site%in%c("PGN1","PGN2"))%>%
  filter(year==2019)%>%
  glimpse()
unique(dat.month$month)
unique(dat.month$year)

## Resample data to select 8 transects per month

dat.month$month<-as.factor(dat.month$month)
levels(dat.month$month)

## Check sampling effort per month

test<-dat.month%>%
  group_by(Site,month)%>%
  summarise(N=n_distinct(id))%>%
  glimpse()

## Firt we convert to long format for each year separately 

library(tidyr)

dataframe=data.frame()

number<-1:99
month<-as.vector(unique(dat.month$month))
site<-as.vector(unique(dat.month$Site))

for(j in number) {
  
  for (i in month) {
    
    for(k in site) {
      
      Abun<-dat.month%>%
        filter(month%in%i)%>%
        filter(Site%in%k)%>%
        dplyr::group_by(id,scientific)%>%
        summarise(number=sum(number))%>%
        ungroup()%>%
        pivot_wider(names_from = "scientific",values_from = number) %>%
        mutate_all(~replace_na(., 0))%>%
        sample_n(8,replace = T)%>%
        glimpse()
      
      ## Reconvert to long format
      
      Abun<-Abun%>%
        pivot_longer(cols = c(-1), names_to = "scientific", values_to = "number")%>%
        glimpse()
      
      ## Create metadata 
      
      metadata<-dat.month%>%
        dplyr::select(Location,Site,month,id)%>%
        glimpse()
      
      new.dat<-left_join(Abun,metadata,by="id")%>%
        distinct(id,scientific,.keep_all = TRUE)%>%
        glimpse()
      
      dataframe<-rbind(dataframe,new.dat)
      
    }
  }
}

## Group data and obtain mean 

test<-dataframe%>%
  group_by(month,Site,scientific)%>%
  summarise(number=mean(number))%>%
  glimpse()

## Recheck sampling effort

dat.month<-test

## Clean working directory

rm(Abun,metadata,new.dat,dataframe)

## Summarize data by month

test<-dat.month%>%
  dplyr::group_by(month,scientific)%>%
  dplyr::summarise(number=sum(number))%>%
  glimpse()

## Check the shape of abundance distributions after transformation of species abundances

plot.march<-ggplot(test%>%filter(month%in%c('3'))%>%filter(number>0),aes(x=reorder(scientific,-number),y=number,fill=scientific))+
  geom_bar(stat = "identity",colour="black",show.legend = F)+
  scale_y_continuous(expand = c(0,0),limits = c(0,45))+
  theme_classic()+
  theme(axis.text.x = element_text(angle=90))

plot.october<-ggplot(test%>%filter(month%in%c('10'))%>%filter(number>0),aes(x=reorder(scientific,-number),y=number,fill=scientific))+
  geom_bar(stat = "identity",colour="black",show.legend = F)+
  scale_y_continuous(expand = c(0,0),limits = c(0,45))+
  theme_classic()+
  theme(axis.text.x = element_text(angle=90))

## Arrange plots in a grob

library(ggpubr)

ggarrange(plot.march,plot.october,nrow=1,ncol=3,align = 'hv')

## Check normal histograms

plot.new()
par(mfrow=c(1,2))

hist(test$number)
hist(log10(test$number+1))

## Create a species x abundances matrix for each year

library(tidyr)
library(dplyr)

abundances<-test%>%
  dplyr::select(month,scientific,number)%>%
  pivot_wider(names_from = "scientific",values_from = number) %>%
  mutate_all(~replace_na(., 0))%>%
  glimpse()

## Convert species to column names

abundances<-data.frame(abundances)
rownames(abundances) <- abundances[,1] #Assigning row names from 1st column 
abundances[,1] <- NULL #Removing the first column

## Transform species abundances

abundances <- log10(abundances+1)

#Order species in the community matrix in same order as traits

z<-as.vector(colnames(abundances))

traits<-traits[z,]

## Check identical species are the same in trait and species matrix

identical(rownames(traits),colnames(abundances))

#######################################
## PLOT SPECIES DENSITY BEFORE/AFTER ##
#######################################

# REPEAT EACH SPECIES ROW AS MANY TIMES AS ITS ABUNDANCE (ON LOG SCALE - MULTIPLIED BY 10 AND ROUNDED)
# TO INCORPORATE ABUNDANCE INTO DENSITY ESTIMATE

n.times <- round(as.numeric(abundances[1,])*10,0)
density_3 <- traits[rep(seq_len(nrow(traits)), n.times),]

n.times <- round(as.numeric(abundances[2,])*10,0)
density_10 <- traits[rep(seq_len(nrow(traits)), n.times),]

## Calculate the community centroid (weighthed mean trait values)

## 3

names(traits)

trait_space <- traits[,c("ThermalMP_5_95","mode.2019")]

centroids.3 <- functcomp(trait_space, as.matrix(abundances[1,]))

## 10

names(traits)

trait_space <- traits[,c("ThermalMP_5_95","mode.2019")]

centroids.10<- functcomp(trait_space, as.matrix(abundances[2,]))

## MARGINAL DENSITY USING COWPLOT ##

## 2017

range(traits$ThermalMP_5_95)
range(traits$MaxLength)

## Relabel factor variables

unique(density_3$herbivore.guild)

density_3<-density_3%>%
  mutate(herbivore.guild=fct_relevel(herbivore.guild,"Browser","Algal farmer","Scraper","Sediment sucker"))%>%
  glimpse()

## Estimate 50%, 75% AND 95% Highest Density Estimates 

getLevel <- function(x,y,prob=c(0.5,0.75,0.95)) {
  kk <- MASS::kde2d(x,y)
  dx <- diff(kk$x[1:2])
  dy <- diff(kk$y[1:2])
  sz <- sort(kk$z)
  c1 <- cumsum(sz) * dx * dy
  approx(c1, sz, xout = 1 - prob)$y
}
L95 <- getLevel(density_3$ThermalMP_5_95,density_3$mode.2019)

## Plot 2D density estimates with marginal densities for each trait axis

names(centroids.3)

d.3 <- ggplot(density_3, aes(ThermalMP_5_95,mode.2019)) + 
  ggtitle('(c)')+
  ylab('Body size (cm)')+
  xlab('Thermal affinity (°C)')+
  stat_density_2d(geom = "polygon", aes(alpha = ..level..),show.legend = F,colour="black",breaks=L95) +
  geom_point(data=density_3, aes(ThermalMP_5_95,mode.2019,
                                    fill=herbivore.guild),colour="black",size=1.5,show.legend = F) +
  geom_point(data=density_3, aes(ThermalMP_5_95,mode.2019,
                                    color=herbivore.guild),size=1,show.legend = F) +
  geom_point(data=centroids.3,aes(x=ThermalMP_5_95,y=mode.2019),shape=23,size=3,fill="orange")+
  geom_vline(xintercept = 23,linetype="dashed",size=1)+
  scale_colour_manual(labels = c("Browser","Algal farmer","Scraper","Sediment sucker"),values=c("yellow4","tomato1","steelblue1","orchid2"))+
  xlim(18.0405,30) + ylim(0,75)+
  theme_classic()+
  Theme1

xdens <- axis_canvas(d.3, axis="x") +
  geom_density(data=density_3, aes(x=ThermalMP_5_95), fill="black",alpha=0.2)

ydens <- axis_canvas(d.3, axis="y", coord_flip = TRUE) +
  geom_density(data=density_3, aes(x=mode.2019), fill="black",alpha=0.2) +
  coord_flip()

p1 <- insert_xaxis_grob(d.3, xdens,position = "top")
p2 <- insert_yaxis_grob(p1, ydens,position = "right")
ggdraw(p2)

## 2019

range(traits$ThermalMP_5_95)
range(traits$MaxLength)

## Relabel factor variables

unique(density_10$herbivore.guild)

density_10<-density_10%>%
  mutate(herbivore.guild=fct_relevel(herbivore.guild,"Sediment sucker","Browser","Algal farmer","Scraper"))%>%
  glimpse()

## Estimate 50%, 75% AND 95% Highest Density Estimates 

getLevel <- function(x,y,prob=c(0.5,0.75,0.95)) {
  kk <- MASS::kde2d(x,y)
  dx <- diff(kk$x[1:2])
  dy <- diff(kk$y[1:2])
  sz <- sort(kk$z)
  c1 <- cumsum(sz) * dx * dy
  approx(c1, sz, xout = 1 - prob)$y
}
L95 <- getLevel(density_10$ThermalMP_5_95,density_10$mode.2019)

## Plot 2D density estimates with marginal densities for each trait axis

d.10 <- ggplot(density_10, aes(ThermalMP_5_95,mode.2019)) + 
  ggtitle('(d)')+
  ylab('')+
  xlab('')+
  stat_density_2d(geom = "polygon", aes(alpha = ..level..),show.legend = F,colour="black",breaks=L95) +
  geom_point(data=density_10, aes(ThermalMP_5_95,mode.2019,
                                    fill=herbivore.guild),colour="black",size=1.5,show.legend = F) +
  geom_point(data=density_10, aes(ThermalMP_5_95,mode.2019,
                                    color=herbivore.guild),size=1,show.legend = F) +
  geom_point(data=centroids.10,aes(x=ThermalMP_5_95,y=mode.2019),shape=23,size=3,fill="orange")+
  geom_vline(xintercept = 23,linetype="dashed",size=1)+
  scale_colour_manual(labels = c("Sediment sucker","Browser","Algal farmer","Scraper"),
                      values=c("orchid2","yellow4","tomato1","steelblue1"))+
  xlim(18.0405,30) + ylim(0,75)+
  theme_classic()+
  Theme1

xdens <- axis_canvas(d.10, axis="x") +
  geom_density(data=density_10, aes(x=ThermalMP_5_95), fill="black",alpha=0.2)

ydens <- axis_canvas(d.10, axis="y", coord_flip = TRUE) +
  geom_density(data=density_10, aes(x=mode.2019), fill="black",alpha=0.2) +
  coord_flip()

p3 <- insert_xaxis_grob(d.10, xdens,position = "top")
p4 <- insert_yaxis_grob(p3, ydens,position = "right")
ggdraw(p4)


## Arrange plots in a grob

library(ggpubr)

ggarrange(p2,p4,ncol=2,nrow = 1)

## Export plot

setwd(plots.dir)
dir()

name<-"Intra_annual_Trait_Space"

ggsave(paste(name,".png",sep="."),width = 19, height = 8,units = "cm",dpi=600)
ggsave(paste(name,".pdf",sep="."),width = 19, height = 8,units = "cm",dpi=600,useDingbats=FALSE)

