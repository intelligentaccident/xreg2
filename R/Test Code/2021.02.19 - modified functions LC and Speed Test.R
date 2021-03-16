

###############################################################
# Kim Rand 2021 
# STANDARD SETUP OF WD AND OPTIONS
###############################################################
# clear? 
rm(list = ls())

options(digits=6)
options(scipen=15)

# Change working directory
oldwd <- getwd() 
# Fetch current script path, either from command line or from Rstudio
get_cur_path <- function() {
  caf <- commandArgs()[substr(commandArgs(),1,7) == '--file=']
  caf <- substr(caf, 8, nchar(caf))
  if (length(caf) == 0) caf <- rstudioapi::getSourceEditorContext()$path
  substr(caf, 1,  nchar(caf) - nchar(strsplit(caf, '.*[/|\\]')[[1]][2]))
}



# Set correct working directory
wd <- get_cur_path()

setwd(wd)

# Set to true to (over-)write results
wrt <- T



###############################################################
# Packages and sources
###############################################################

# Install pacman package if necessary
if(class(try(find.package("pacman"), silent = TRUE)) == "try-error") install.packages("pacman")
library(pacman)

# Load packages, install at need
p_load(
  readstata13,  # To read STATA data files
  readxl,       # To read xls
  writexl,      # To write xls
  
  # ggplot2,      # For graphs 
  # cowplot,      # For graphs
  # grid,         # For graphs
  # gridExtra,    # For graphs
  
  parallel,       # For parallelization
  Rfast,          # C++ optimized functions
  microbenchmark, # For benchmarking
  
  GenSA,          # Simulated Annealing
  
  optimx,          # Better optim options
  tidyverse,       # handling data
  functools,   # Vapply (speed optimizied versions of apply functions)
  Rcpp,            # include C++ code
  profvis          # speed tests 
  
)  


# Load xreg, install from github at need
p_load_gh('intelligentaccident/xreg')

# Import functions for easy handling of EQ-5D states
source('2021.02.10 - step_list_functions.R')
source('2021.02.12 - ll_functions.R')
source('2021.02.19 - ll_functions_Rcpp.R')
source('2021.02.19 - LC_tests.R') # in this scrip you find the data simulation and the 'original functions of LC'



#############################################################
#### SPEED UP, MODIFIED FUNCTIONS #### (names are the same as above but with an additional s at the end)
#############################################################

#### FIND GOOD START VALUES ###
# https://www.sciencedirect.com/science/article/abs/pii/S0167947302001639?via%3Dihub


# TBD


### LC_fits ###
LC_fits <- function(par, dtalist, maxit = 10, i = 1) {
  
  stv <- par
  ll_vector <- NULL
  
  while(i <= maxit){ # use while instead of for loop (doesn´t seem to make a real speed difference)
    print(i)
    print(paste0("Run ", i))
    print("STV: ")
    print(stv)    
    
    E_step <- LC_list2s(par = stv, dtalist = dtalist) # use modified LC_list2
    
    print(paste0("Probs: "))
    print(E_step)
    # at the moment no control argument used in ucminf 
    M_step <- ucminf::ucminf(par = stv, fn = LC_list2s, hessian = 0, dtalist = dtalist, probs =E_step)
    stv <- M_step$par
    
    print(paste0("VALUE = ", M_step$value))
    print(M_step$counts)
    
    # stop if convergence is reached
    ll_vector <- c(ll_vector,M_step$value)
    if(i > 1){
      if(abs(ll_vector[i]- ll_vector[i-1]) < 1e-6){
        break
      }
    } 
    
    i <- i +1
  }
  return(M_step)
}



### LC_list2s ###
LC_list2s <- function (par, dtalist, probs = NULL) {
  

  
  tout <- sapply(dtalist, function(dtl) {  # parallelize this with parallel::parSapply ?
    grfn <- paste0(dtl$type,"3")
    do.call(grfn, args = list(dtalist = dtl, par = par, log.p = F, aggregate.p = F))
  })
  if(is.null(probs)) {
    touts <- tout %*% c(1,1)
    tout <- colmeans(tout/as.vector(touts))
  } else {
    tout <- -sum(log(tout %*% probs))
  }
  
  
  return(tout)
}


#### Ll FUNCTION ###


# Negloglik for continuous normal, with censoring, intervals ...
cont_normal3 <- function (par, dtalist, log.p = T, aggregate.p = T) {
  
  
  p_type = dtalist$type
  # Set names to parameters
  names(par) <- dtalist$coeff_names
  
  # Fetch formulas
  formulas <- dtalist$formulas
  
  # Calculations on unique observations
  d_df <- dtalist$unique#[1:10,] ############ ONLY FOR TEST WITH 10 ROWS
  f_df <- dtalist$full
  
  nobs <- NROW(d_df)
  fobs <- NROW(f_df)
  
  # Fetch parameters and fixed values
  args <- c(as.list(par), as.list(dtalist$fixed))
  # print(args)
  
  '
   ### OLD CODE ###
  # Sequence over formulas
  for (formulanum in 1:length(formulas)) {
    formula <- formulas[[formulanum]]
    new_var <- as.character(formula[[2]])
    if (new_var %in% colnames(d_df))
      new_var = "Xb"
    d_df[, new_var] <- with(args, with(d_df, eval(parse(text = formula))))
  }
'
  #### USE ONLY ONE OPTION TO MODIFY THE FOR LOOP ABOVE ######
  ### MODIFICATION OPTION 1
  ### if using this block (line 170 to 178), then comment out line 182 to 183
  
  #  names <- c()
  #  out <- sapply(formulas, function(x){
  #    if(is.element(as.character(x[[2]]),colnames(d_df))){
  #      names <<- append(names,"Xb")
  #    }else{names <<- append(names, as.character(x[[2]]))}
  #   # names <<- if_else(is.element(as.character(x[[2]]),colnames(d_df)),append(names, "Xb"),append(names, as.character(x[[2]])))
  #    with(args, with(d_df, eval(parse(text = x))))})
  #
  #  d_df[,names] <- as.data.frame(out)
  
  
  ### MODIFICATION OPTION 2 
  out <- lapply(formulas, function(x){
    with(args, with(d_df, eval(parse(text = x))))})

  ############################  
  '  
  ### OLD CODE ##
  pvals <- matrix(NA, nrow = nobs, ncol = 1)
  sel <- d_df$internal_type == 1 # Discrete
  pvals[sel,] <-  dnorm(x = d_df$Xb[sel], 
                        mean = d_df$internal_lb[sel], 
                        sd = d_df$sigma_est[sel], 
                        log = log.p)
  sel <- !sel # Intervals
  
  # Assignment appears to be slow, so separating out saves ~50 microseconds per run
  if(log.p) {
    pvals[sel,] <- log(pnorm(q = d_df$Xb[sel], mean = d_df$internal_lb[sel], sd = d_df$sigma_est[sel], lower.tail = T, log.p = F) -
                         pnorm(q = d_df$Xb[sel], mean = d_df$internal_ub[sel], sd = d_df$sigma_est[sel], lower.tail = T, log.p = F))
  } else {
    pvals[sel,] <- pnorm(q = d_df$Xb[sel], mean = d_df$internal_lb[sel], sd = d_df$sigma_est[sel], lower.tail = T, log.p = F) -
      pnorm(q = d_df$Xb[sel], mean = d_df$internal_ub[sel], sd = d_df$sigma_est[sel], lower.tail = T, log.p = F)
  }
'
  # MODIFICATION to create pvals #
  # use lines 209 tp 215 with modificaiton Option 1 (see above) 
  # and lines 217 to 223 with modification Option 2
  
  #  pvals <- create_pval( x = d_df$Xb,
  #             means = d_df$internal_lb,
  #             means2 = d_df$internal_ub, 
  #             sds = d_df$sigma_est,
  #             log = rep(log.p,length(d_df$Xb)),
  #             sel = d_df$internal_type,
  #             lowerT = rep(TRUE,length(d_df$Xb)))
  
  pvals <- create_pval( x = out$main, # problem here is that this is not really dynamically !!!!
                        means = d_df$internal_lb,
                        means2 = d_df$internal_ub, 
                        sds = rep(out$SIGMA,length(out$main)),
                        log = rep(log.p,length(out$main)),
                        sel = d_df$internal_type,
                        lowerT = rep(TRUE,length(out$main)))
  
  
  
  
  
  
  # Separate calculations for right- and left- censored values. 
  # Speed gain is approximately 0 in testing
  # sel <- d_df$internal_type == 0 # Right censored
  # d_df$p[sel] <- pnorm(q = d_df$Xb[sel], mean = d_df$internal_lb[sel], sd = d_df$sigma_est[sel], lower.tail = T, log.p = T)
  # sel <- d_df$internal_type == 2 # Left censored
  # d_df$p[sel] <- pnorm(q = d_df$Xb[sel], mean = d_df$internal_ub[sel], sd = d_df$sigma_est[sel], lower.tail = F, log.p = T)
  
  # Return all p-values if aggregate.p == F
  if(!aggregate.p) return(pvals)
  
  if(log.p) {
    # Set NA to minimum if boundNA
    if("boundNA" %in% names(dtalist)) if(dtalist$boundNA) pvals[is.na(pvals)] <- log(.Machine$double.xmin)
    
    # Move -Inf to log of minimum
    pvals[pvals == -Inf] <- log(.Machine$double.xmin)
    
    # Attach unique-observation p-values to all observations
    
    pvals <- pvals[f_df$umatch,]
    
    # Negloglik
    pvals <- -sum(pvals)
  } else {
    # Set NA to minimum if boundNA
    if("boundNA" %in% names(dtalist)) if(dtalist$boundNA) pvals[is.na(pvals)] <- .Machine$double.xmin
    
    # Move -Inf to log of minimum
    pvals[pvals == 0] <- .Machine$double.xmin
    
    # Attach unique-observation p-values to all observations
    
    pvals <- pvals[f_df$umatch,]
    
    # Negloglik
    pvals <- -prod(pvals)
  }
  return(pvals)
}

########################################################################################

# use these elements to go through cont_normal function manually

#par <- stv
#dtalist <- dummydtalist$var1
#log.p <-  FALSE
#aggregate.p <- FALSE


#  check if results from cont_normal2 and cont_normal3 are the same
cont_normal2(stv,dummydtalist$var1, log.p = FALSE, aggregate.p = FALSE)
cont_normal3(stv,dummydtalist$var1, log.p = FALSE, aggregate.p = FALSE)

# if this is TRUE, then cont_normal2 and cont_normal3 gave same result
all(cont_normal2(stv,dummydtalist$var1, log.p = FALSE, aggregate.p = FALSE) == cont_normal3(stv,dummydtalist$var1, log.p = FALSE, aggregate.p = FALSE)) 


###############################################################################
### Speed Tests ###
###############################################################################


### ll function
microbenchmark(
  cont2 = cont_normal2(stv,dummydtalist$var1, log.p = FALSE, aggregate.p = FALSE),
  cont3 = cont_normal3(stv,dummydtalist$var1, log.p = FALSE, aggregate.p = FALSE)
)


### LC_list2
microbenchmark(
  LC_list2 = LC_list2(stv, dummydtalist, probs = NULL), # with cont_normal2
  LC_list2s = LC_list2s(stv, dummydtalist, probs = NULL) # with cont_normal3 (optimizied) and Vapply 
)



### LC_fit
microbenchmark(
  LC_fit  = LC_fit(stv, dummydtalist, maxit = 5), # with cont_normal2
  LC_fits = LC_fits(stv, dummydtalist, maxit = 5), # with cont_normal3 (optimizied)
  times =1
)




#### PROFVIS ###
# find out where calculation is slow
# If you just type in LC_fit(stv, dummydtalist, maxit = 5) or 
# LC_fits(stv, dummydtalist, maxit = 5)#
# it does not show me a flame graph (I don´t know why)
# therefor I just copied the sourceCode here and used it

profvis({
  par <- stv
  dtalist <- dummydtalist
  maxit = 1
  i = 1
  
  stv <- par
  while (i <= maxit) {
    print(i)
    print(paste0("Run ", i))
    print("STV: ")
    print(stv)
    
    E_step <- LC_list2s(par = stv, dtalist = dtalist) # use LC_list2 or LC_list2s
    
    print(paste0("Probs: "))
    print(E_step)
    M_step <-
      optimr(
        par = stv,
        fn = LC_list2s, # use LC_list2 or LC_list2s
        hess = F,
        method = "BFGS",
        control = ctrl,
        dtalist = dtalist,
        probs = E_step
      )
    stv <- M_step$par
    
    print(paste0("VALUE = ", M_step$value))
    print(M_step$counts)
    i <- i + 1
  }
})


#################### TEST ######################

