library(optimx)
library(microbenchmark)
require(Rcpp)
library(Rfast)


load("~/R/xreg/3_Code/Codebeispiele von Kim/Xreg new/Optimizing_cont_normal3.RData")

par <- stv
dtalist <- dummydtalist$var1
log.p <-  FALSE
aggregate.p <- FALSE


### FUNCTION TO BE OPTIMIZED ###
# cont_normal2 is the original Version of cont_normal3

# INPUT
# par = start Values, named vector
# dtalist = list contaning specific elemts as formulas to be evaluted


# Negloglik for continuous normal, with censoring, intervals ...
cont_normal3 <- function (par, dtalist, log.p = T, aggregate.p = T) {
  
  p_type <- dtalist$type
  
  # Set names to parameters
  names(par) <- dtalist$coeff_names
  
  # Fetch formulas
  formulas <- dtalist$formulas
  
  # Calculations on unique observations
  d_df <- dtalist$unique
  f_df <- dtalist$full
  
  nobs <- NROW(d_df)
  fobs <- NROW(f_df)
  
  # Fetch parameters and fixed values
  args <- c(as.list(par), as.list(dtalist$fixed))
  # print(args)
  

   ### the following for loop seems to be slow ###
  # using lapply is possible, but the dynamic naming (Xb) is a problem then
  
  # Sequence over formulas
  for (formulanum in 1:length(formulas)) {
    formula <- formulas[[formulanum]]
    new_var <- as.character(formula[[2]])
    if (new_var %in% colnames(d_df))
      new_var = "Xb"
    d_df[, new_var] <- with(args, with(d_df, eval(parse(text = formula))))
  }

  
  ### This part (lines 57 to 74) is availabe in Rcpp now ##
  #pvals <- matrix(NA, nrow = nobs, ncol = 1)
  #sel <- d_df$internal_type == 1 # Discrete
  #pvals[sel,] <-  dnorm(x = d_df$Xb[sel], 
  #                      mean = d_df$internal_lb[sel], 
  #                      sd = d_df$sigma_est[sel], 
  #                      log = log.p)
  #sel <- !sel # Intervals
  
  # Assignment appears to be slow, so separating out saves ~50 microseconds per run
  #if(log.p) {
  #  pvals[sel,] <- log(pnorm(q = d_df$Xb[sel], mean = d_df$internal_lb[sel], sd = d_df$sigma_est[sel], lower.tail = T, log.p = F) -
  #                       pnorm(q = d_df$Xb[sel], mean = d_df$internal_ub[sel], sd = d_df$sigma_est[sel], lower.tail = T, log.p = F))
  #} else {
  #  pvals[sel,] <- pnorm(q = d_df$Xb[sel], mean = d_df$internal_lb[sel], sd = d_df$sigma_est[sel], lower.tail = T, log.p = F) -
  #    pnorm(q = d_df$Xb[sel], mean = d_df$internal_ub[sel], sd = d_df$sigma_est[sel], lower.tail = T, log.p = F)
  #}

  # use this in combination with the new Rcpp function create_pval #
  # instead of lines 57 to 74

    pvals <- create_pval( x = d_df$Xb,
               means = d_df$internal_lb,
               means2 = d_df$internal_ub, 
               sds = d_df$sigma_est,
               log = rep(log.p,length(d_df$Xb)),
               sel = d_df$internal_type,
               lowerT = rep(TRUE,length(d_df$Xb)))
  

  
  # Separate calculations for right- and left- censored values. 
  # Speed gain is approximately 0 in testing
  # sel <- d_df$internal_type == 0 # Right censored
  # d_df$p[sel] <- pnorm(q = d_df$Xb[sel], mean = d_df$internal_lb[sel], sd = d_df$sigma_est[sel], lower.tail = T, log.p = T)
  # sel <- d_df$internal_type == 2 # Left censored
  # d_df$p[sel] <- pnorm(q = d_df$Xb[sel], mean = d_df$internal_ub[sel], sd = d_df$sigma_est[sel], lower.tail = F, log.p = T)
  
    
# For the moment it is enough to use aggregate.p = F

    
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






### PARTS in Rcpp ####

### function to use in cont_normal3 (creating pvals matrix)
sourceCpp( code = '#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
arma::vec create_pval(arma::vec x, arma::vec means, arma::vec means2, arma::vec sds, arma::uvec log, arma::vec sel, arma::uvec lowerT) {
  int n = x.size();
  arma::vec res(n);
  for(int i = 0; i < n; i++) {
    if(sel[i] = 1){
      res[i] = R::dnorm(x[i], means[i], sds[i], log[i]);
    } else{ 
            if(log[i]){
                res[i] = log(R::qnorm(x[i], means[i], sds[i], lowerT[i], false) -  R::qnorm(x[i], means2[i], sds[i], lowerT[i], false)) ;
            } else{
                res[i] = R::qnorm(x[i], means[i], sds[i], lowerT[i], false) -  R::qnorm(x[i], means2[i], sds[i], lowerT[i], false) ;
            }
    }
  }
  return res;
}')


### FOR LOOP in Rcpp ### 
# not included in cont_normal3 yet
cppFunction( '
DataFrame getXb(List formulas, DataFrame df){
  int size = df.nrow();
  int lenL = formulas.size();
  List out(lenL);

  for(int i = 0; i < lenL; i++) {
    ExpressionVector exp = formulas[i];
    out[i] = exp.eval(df);
  }
    out.attr("class") = "data.frame";
  out.attr("row.names") = seq(1, size);
  out.attr("names") = CharacterVector::create("main","SIGMA");
return(out);
}
')

# USE THIS TO TEST getXb
# to be able to use the Rcpp Function getXb we need an additional lapply (for parsing)
# and we need to cbind d_df with the values of par

formulas <- dtalist$formulas
d_df <- cbind(dtalist$unique, c(as.list(par), as.list(dtalist$fixed)))
names(par) <- dtalist$coeff_names
head(getXb(formulas = lapply(formulas,function(x)parse(text = x)), df = d_df))




### PROBLEMS ####
# how to parse in Rcpp?
# at the moment formulas have to be parsed in R and can then be used in Rcpp
# actual solution is not faster than using lapply in cont_normal3 instead of the loop
# how to cahnge names of output dynamical (line 176 should not always be names main and SIGMA)





######## TEST ########

#  check if results from cont_normal2 and cont_normal3 are the same
cont_normal2(stv,dummydtalist$var1, log.p = FALSE, aggregate.p = FALSE) # original
cont_normal3(stv,dummydtalist$var1, log.p = FALSE, aggregate.p = FALSE) # modified

# if this is TRUE, then cont_normal2 and cont_normal3 gave same result
all(cont_normal2(stv,dummydtalist$var1, log.p = FALSE, aggregate.p = FALSE) == cont_normal3(stv,dummydtalist$var1, log.p = FALSE, aggregate.p = FALSE)) 



###### USE of cont_normal function(s) ######

LC_list2s <- function (par, dtalist, probs = NULL) {
  
  tout <- sapply(dtalist, function(dtl) {  # parallize this with parallel::parSapply ? seems to be slower?
    grfn <- paste0(dtl$type,"3") # change this to 2, if cont_normal2 should be used
    do.call(grfn, args = list(dtalist = dtl, par = par, log.p = F, aggregate.p = F)) # use of cont_normal2 or cont_normal2
  })
  if(is.null(probs)) {
    touts <- tout %*% c(1,1)
    tout <- colmeans(tout/as.vector(touts))
  } else {
    tout <- -sum(log(tout %*% probs))
  }
  
  
  return(tout)
}



### STEPS IN EM ALGO ###

E_step <- LC_list2s(stv,dummydtalist)
#  M_step <- optimr(par = stv, fn = LC_list2s, hess = F, method = "BFGS",  dtalist = dummydtalist, probs =E_step)

# speed up M step with ucminf (install package ucminf)
M_step <- ucminf::ucminf(par = stv, fn = LC_list2s, hessian = 0, dtalist = dummydtalist, probs =E_step)






#### SPEED TESTS ###

microbenchmark(cont2 = cont_normal2(stv,dummydtalist$var1, log.p = FALSE, aggregate.p = FALSE), # original
               cont3 = cont_normal3(stv,dummydtalist$var1, log.p = FALSE, aggregate.p = FALSE) # modified
)

microbenchmark(LC_list2 = LC_list2(stv, dummydtalist, probs = NULL), # original
               LC_list2s = LC_list2s(stv, dummydtalist, probs = NULL) # modified
)

# in LC_list and LC_list2 sapply loops over dummydtalist (which is a list of lists and contains the data and the formulas)
# it is possible to use different formulas and distribuitions on the different elementes of dummdytalist
# in the M step of EM Algorithm LC_list(s) is the function to be optimized

microbenchmark(M_step_optimr =  optimr(par = stv, fn = LC_list2, hess = F, method = "BFGS",  dtalist = dummydtalist, probs =E_step), # original
               M_step_ucminf = ucminf::ucminf(par = stv, fn = LC_list2s, hessian = 0, dtalist = dummydtalist, probs =E_step), # modified
                times = 5
)




### single core vs. multicore sapply (as used in LC_list to calculate ll) ###

library(parallel)

microbenchmark(
  parallel = {
    n.cores <- detectCores()
    n.cores
    clust <- makeCluster(n.cores)
    clusterExport(clust, varlist = c(ls()))
    
    parSapply(clust, dummydtalist, function(dtl){
      grfn <- paste0(dtl$type,"2") # doesn´t work with cont_normal3 at the moment
      do.call(grfn, args = list(dtalist = dtl, par = par, log.p = F, aggregate.p = F))})
    stopCluster(clust)
  },
  standard = {
    sapply(dummydtalist, function(dtl) {
      grfn <- paste0(dtl$type,"2")
      do.call(grfn, args = list(dtalist = dtl, par = par, log.p = F, aggregate.p = F))})
  },
  times = 5
)





########################## TEST ##############

test_df <- cbind(d_df,args)
with(test_df,eval(parse(text = formula)))
test_df[,"Xb"] <- with(test_df,eval(parse(text = formula)))
0.4 + 4 *3 +2 *0.94014924