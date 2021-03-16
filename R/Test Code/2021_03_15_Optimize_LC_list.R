
library(optimx)
library(microbenchmark)
require(Rcpp)
library(Rfast)


load("~/R/xreg/3_Code/Codebeispiele von Kim/Xreg new/Optimizing_LC.RData")


data <- dummydtalist$var1$unique
dtalist <- dummydtalist$var1
dtalist$formulas <- list(form1 = list(main1 = depvar ~ INTERCEPT + d1 *D1 + d2 * D2, SIGMA1 = sigma_est ~ exp(LN_SIGMA)),
                         form2 = list(main2 = depvar ~ INTERCEPT.1 + d1 * D1.1 + d2 * D2.1, SIGMA2 = sigma_est ~ exp(LN_SIGMA)))

log.p <- F
aggregate.p <- F


LC <- function(par, data,dtalist, probs = NULL){
  
  p_type = dtalist$type
  # Set names to parameters
  names(par) <- dtalist$coeff_names
  
  # Fetch formulas
  formulas <- dtalist$formulas
  
  # Calculations on unique observations
  d_df <- data
  
  nobs <- NROW(d_df)

  
  # Fetch parameters and fixed values
  args <- c(as.list(par), as.list(dtalist$fixed))
  # print(args)
  
  
  ### MODIFICATION OPTION 2 
  out <- lapply(formulas, function(y){
            lapply(y, function(x){
              with(args, with(d_df, eval(parse(text = x))))
              })
           
        })

  
  ### here another create_pval will be needed, if distribution is not cont_normal
  # maybe include an indicator in formulas list which function shpuld be used here?
  tout <- sapply(out,function(x){
    create_pval( x = x[[1]], # problem here is that this is not really dynamically !!!!
                 means = d_df$internal_lb,
                 means2 = d_df$internal_ub, 
                 sds = rep(x[[2]],length(x[[1]])),
                 log = rep(log.p,length(x[[1]])),
                 sel = d_df$internal_type,
                 lowerT = rep(TRUE,length(x[[1]])))
  })

  if(is.null(probs)) {
    touts <- tout %*% c(1,1)
    tout <- colmeans(tout/as.vector(touts))
  } else {
    tout <- -sum(log(tout %*% probs))
  }
  
  return(tout)
}


microbenchmark(new = LC(stv,data,dtalist = dtalist),
               old = LC_list2s(stv,dummydtalist))


E_step <- LC(stv,data,dtalist = dtalist)
microbenchmark(M_step_new =  ucminf::ucminf(par = stv, fn = LC, hessian = 0, dtalist = dtalist, data = data, probs =E_step), 
               M_step_old= ucminf::ucminf(par = stv, fn = LC_list2s, hessian = 0, dtalist = dummydtalist, probs =E_step), 
               times = 5
)





profvis({
 # par <- stv
#  dtalist <- dummydtalist
  maxit = 1
  i = 1
  
  stv <- par
  while (i <= maxit) {
    print(i)
    print(paste0("Run ", i))
    print("STV: ")
    print(stv)
    
    E_step <- LC(stv,data,dtalist = dtalist)
    
    print(paste0("Probs: "))
    print(E_step)
    M_step <- ucminf::ucminf(par = stv, fn = LC, hessian = 0, dtalist = dtalist, data = data, probs =E_step)
    stv <- M_step$par
    
    print(paste0("VALUE = ", M_step$value))
    print(M_step$counts)
    i <- i + 1
  }
})


