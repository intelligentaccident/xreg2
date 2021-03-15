

# Handling of ll-functions
negll_list2 <- function (par, dtalist, ...) {
  sum(sapply(dtalist, function(dtl) do.call(dtl$ll_fun, args = c(list(dtalist = dtl, par = par), list(...)))))
}



# Handling of gradient functions
gradll_list2 <- function (par, dtalist, ...) {
  rowsums(sapply(dtalist, function(dtl) do.call(dtl$grad_fun, args = c(list(dtalist = dtl, par = par), list(...)))))
}


# Forward gradient function for continuous normal
grad_cont_normal2 <- function (par, dtalist, ret_mat = F, eps = 1e-7, centered = T) {
  
  if(centered) par <- par-eps/2
  
  # Calculations on unique observations
  d_df <- dtalist$unique
  f_df <- dtalist$full
  nobs <- NROW(d_df)
  fobs <- NROW(f_df)
  
  
  
  
  # Set names to parameters
  names(par) <- dtalist$coeff_names
  plen <- length(par)
  plen1 <- plen +1
  # plen2 <- 2*plen
  
  parm <- as.data.frame(as.list(par))
  parm <- parm[rep(1, plen+1),]
  
  for(pn in 1:plen) {
    theserows <- pn+1
    parm[theserows, pn] <- parm[theserows, pn]+eps
  }
  par <- parm
  # Fetch formulas
  formulas <- dtalist$formulas
  
  
  # Fetch parameters and fixed values
  args <- c(as.list(par), as.list(dtalist$fixed))
  
  # Sequence over formulas
  for (formulanum in 1:length(formulas)) {
    formula <- formulas[[formulanum]]
    new_var <- as.character(formula[[2]])
    if (new_var %in% c(colnames(d_df), dtalist$valuevar)) new_var = "Xb"
    tmp <- sapply(1:plen1, function(x) {
      theseargs <- c(as.list(par), as.list(dtalist$fixed), as.list(d_df[,colnames(d_df)[colnames(d_df) %in% all.vars(formula)], drop = F]))
      theseargs <- lapply(theseargs, function(thisarg) {
        
        if(!is.null(dim(thisarg))) if(dim(thisarg)[2] == plen1) return(thisarg[,x])
        if(length(thisarg) == plen1) return(thisarg[x])
        thisarg
      })
      with(theseargs, eval(parse(text = formula)))
    })
    if(is.null(dim(tmp)) & length(tmp) == plen1) tmp <- matrix(tmp, nrow = nobs, ncol = plen1, byrow = T)
    d_df[, new_var] <- tmp
  }
  
  # matrix of p-calculations
  pvals <- matrix(NA,nrow = nobs, ncol = plen1)
  
  sel <- d_df$internal_type == 1 # Discrete
  pvals[sel,] <- dnorm(x = d_df$Xb[sel,], 
                       mean = d_df$internal_lb[sel], 
                       sd = d_df$sigma_est[sel,], 
                       log = TRUE)
  sel <- !sel # Intervals and censored
  pvals[sel,] <- log(pnorm(q = d_df$Xb[sel,], mean = d_df$internal_lb[sel], sd = d_df$sigma_est[sel,], lower.tail = T, log.p = F) -
                       pnorm(q = d_df$Xb[sel,], mean = d_df$internal_ub[sel], sd = d_df$sigma_est[sel,], lower.tail = T, log.p = F))
  
  
  
  # Set NA to minimum if boundNA
  if("boundNA" %in% names(dtalist)) if(dtalist$boundNA) pvals[is.na(pvals)] <- log(.Machine$double.xmin)
  
  # Move -Inf to log of minimum
  pvals[pvals == -Inf] <- log(.Machine$double.xmin)
  
  
  
  # Attach unique-observation p-values to all observations
  pvals <- pvals[dtalist$full$umatch,]
  if(ret_mat) return(pvals)
  pvals <- colSums(pvals)
  
  # Subtract modifications from current, set to initial scale
  pvals <- (pvals[1]-pvals[2:plen1])/eps
  
  
  return(pvals)
}


# Forward gradient function for dichotomous logistic
grad_dich_logistic2 <- function (par, dtalist, ret_mat = F, eps = 1e-7, centered = T) {
  
  if(centered) par <- par-eps/2
  # magnitude of change
  # eps <- c(- eps, eps)
  
  # Calculations on unique observations
  d_df <- dtalist$unique
  f_df <- dtalist$full
  nobs <- NROW(d_df)
  fobs <- NROW(f_df)
  
  
  
  
  # Set names to parameters
  names(par) <- dtalist$coeff_names
  plen <- length(par)
  plen1 <- plen +1
  # plen2 <- 2*plen
  
  parm <- as.data.frame(as.list(par))
  parm <- parm[rep(1, plen+1),]
  
  for(pn in 1:plen) {
    theserows <- pn+1
    parm[theserows, pn] <- parm[theserows, pn]+eps
  }
  par <- parm
  # Fetch formulas
  formulas <- dtalist$formulas
  
  
  # Fetch parameters and fixed values
  args <- c(as.list(par), as.list(dtalist$fixed))
  
  # Sequence over formulas
  for (formulanum in 1:length(formulas)) {
    formula <- formulas[[formulanum]]
    new_var <- as.character(formula[[2]])
    if (new_var %in% c(colnames(d_df), dtalist$valuevar)) new_var = "Xb"
    # d_df[, new_var] <- with(args, with(d_df, eval(parse(text = formula))))
    tmp <- sapply(1:plen1, function(x) {
      theseargs <- c(as.list(par[x,]), as.list(dtalist$fixed))
      with(theseargs, with(d_df, eval(parse(text = formula))))
    })
    if(is.null(dim(tmp)) & length(tmp) == plen1) tmp <- matrix(tmp, nrow = nobs, ncol = plen1, byrow = T)
    d_df[, new_var] <- tmp
  }
  
  # matrix of p-calculations
  
  Xb <- d_df$Xb
  if("theta_est" %in% colnames(d_df)) Xb <- Xb*d_df$theta_est
  
  # Calculate logistic
  # logistic_tmp <- (1/(1+exp(-Xb))) # slower than equivalent using hyperbolic tangent
  logistic_tmp <- .5+.5*tanh(Xb/2)
  pvals <- log(d_df$internal_lb *    logistic_tmp + 
                 (1-d_df$internal_lb)* (1-logistic_tmp))
  pvals <- matrix(pvals, nrow = nobs, ncol = plen1)
  
  # Set NA to minimum if boundNA
  if("boundNA" %in% names(dtalist)) if(dtalist$boundNA) pvals[is.na(pvals)] <- log(.Machine$double.xmin)
  
  # Move -Inf to log of minimum
  pvals[pvals == -Inf] <- log(.Machine$double.xmin)
  
  
  
  # Attach unique-observation p-values to all observations
  pvals <- pvals[dtalist$full$umatch,]
  if(ret_mat) return(pvals)
  pvals <- colSums(pvals)
  
  # Subtract modifications from current, set to initial scale
  pvals <- (pvals[1]-pvals[2:plen1])/eps
  
  
  return(pvals)
}


# Negloglik for continuous normal, with censoring, intervals ...
cont_normal2 <- function (par, dtalist, log.p = T, aggregate.p = T, print_pars = F) {
  if(print_pars) message(paste(par, collapse = "\t"))
  
  p_type = dtalist$type
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
  
  # Sequence over formulas
  for (formulanum in 1:length(formulas)) {
    formula <- formulas[[formulanum]]
    new_var <- as.character(formula[[2]])
    if (new_var %in% c(colnames(d_df), dtalist$valuevar)) new_var = "Xb"
    d_df[, new_var] <- with(args, with(d_df, eval(parse(text = formula))))
  }
  
  
  
  
  pvals <- matrix(NA, nrow = nobs, ncol = 1)
  sel <- d_df$internal_type == 1 # Discrete
  pvals[sel,] <- dnorm(x = d_df$Xb[sel], 
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

attr(x = cont_normal2, which = "required") <- list(sigma_est = list(formula = formula(sigma_est ~ exp(LN_SIGMA), env = globalenv()), 
                                                                    start_values = c(LN_SIGMA = 0)))





# Negloglik for logistic
dich_logistic2 <- function (par, dtalist) {
  
  
  p_type = dtalist$type
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
  
  # Sequence over formulas
  for (formulanum in 1:length(formulas)) {
    formula <- formulas[[formulanum]]
    new_var <- as.character(formula[[2]])
    if (new_var %in% c(colnames(d_df), dtalist$valuevar)) new_var = "Xb"
    d_df[, new_var] <- with(args, with(d_df, eval(parse(text = formula))))
  }
  
  
  
  Xb <- d_df$Xb
  if("theta_est" %in% colnames(d_df)) Xb <- Xb*d_df$theta_est
  
  # Calculate logistic
  # logistic_tmp <- (1/(1+exp(-test2$Xb))) # slower than equivalent using hyperbolic tangent
  logistic_tmp <- .5+.5*tanh(Xb/2)
  pvals <- log(d_df$internal_lb *    logistic_tmp + 
                 (1-d_df$internal_lb)* (1-logistic_tmp))
  pvals <- matrix(pvals, nrow = nobs, ncol = 1)
  
  
  
  # Set NA to minimum if boundNA
  if("boundNA" %in% names(dtalist)) if(dtalist$boundNA) pvals[is.na(pvals)] <- log(.Machine$double.xmin)
  
  # Move -Inf to log of minimum
  pvals[pvals == -Inf] <- log(.Machine$double.xmin)
  
  # Attach unique-observation p-values to all observations
  
  pvals <- pvals[f_df$umatch,]
  
  # Negloglik
  tout <- -sum(pvals)
  return(tout)
}



# Negloglik for continuous normal, with random intercept in internal_id
cont_r_normal2 <- function (par, dtalist) {
  
  # Set names to parameters
  names(par) <- dtalist$coeff_names
  
  # Fetch formulas
  formulas <- dtalist$formulas
  
  # Calculations on unique observations
  d_df <- dtalist$unique
  f_df <- dtalist$full
  
  # Fetch parameters and fixed values
  args <- c(as.list(par), as.list(dtalist$fixed))
  
  # Sequence over formulas
  for (formulanum in 1:length(formulas)) {
    formula <- formulas[[formulanum]]
    new_var <- as.character(formula[[2]])
    if (new_var %in% c(colnames(d_df), dtalist$valuevar)) new_var = "Xb"
    d_df[, new_var] <- with(args, with(d_df, eval(parse(text = formula))))
  }
  
  nobs <- NROW(d_df)
  fobs <- NROW(f_df)
  
  
  if(!sum(c("LN_OMEGA", "OMEGA") %in% names(par))== 1) stop("Ambiguous or non-unique OMEGA/LN_OMEGA paramter defined.")
  if("LN_OMEGA" %in% names(par)) betwSD <- exp(par['LN_OMEGA'])
  if("OMEGA" %in% names(par)) betwSD <- par['OMEGA']
  
  uid <- sort(unique(f_df$internal_id))
  nuid <- NROW(uid)
  
  pvals <- rep(0, nuid)
  
  
  GH_order <- 8L
  if(!is.environment((xregenv <- getOption('xreg')))) {
    options(xreg = (xregenv = new.env(parent = globalenv())))
  }
  if (exists("Gauss_Hermite_Order", envir = xregenv)) {
    GH_order <- get("Gauss_Hermite_Order", envir = xregenv)
  }  else {
    assign(x = "Gauss_Hermite_Order", value = GH_order, envir = xregenv)
  }
  ws <- as.matrix(xreg:::xreg_GHS[[GH_order]]$w)
  xs <- xreg:::xreg_GHS[[GH_order]]$x
  mus <- sqrt(2) * betwSD * xs
  tmpm <- matrix(NA, nrow = nobs, ncol = GH_order)
  Xbm <- Outer(y= d_df$Xb, x = mus, "+")
  
  sel <- d_df$internal_type == 1 # Discrete
  sel2 <- !sel
  tmpm[sel, ] <- dnorm(x = d_df$internal_lb[sel], 
                       mean = Xbm[sel,], 
                       sd = d_df$sigma_est[sel], 
                       log = T)
  tmpm[sel2, ] <- Log(pnorm(q = (d_df$internal_lb[sel2]), 
                            mean = Xbm[sel2,], 
                            sd = d_df$sigma_est[sel2], lower.tail = F) - 
                        pnorm(q = (d_df$internal_ub[sel2]), 
                              mean = Xbm[sel2,], 
                              sd = d_df$sigma_est[sel2], 
                              lower.tail = F))
  
  # Expand to full dataset
  tmpm <- tmpm[f_df$umatch,]
  
  # Calculate individual likelihood contributions
  pvals <- Log(1/sqrt(pi) * 
                 exp(apply(X = tmpm, 
                           MARGIN = 2, 
                           FUN = group.sum, 
                           ina = f_df$internal_id) + 
                       rep_row(Log(ws), nuid)) %*% rep(1, GH_order)) 
  
  
  
  # Set NA to minimum if boundNA
  if("boundNA" %in% names(dtalist)) if(dtalist$boundNA) pvals[is.na(pvals)] <- log(.Machine$double.xmin)
  
  # Move -Inf to log of minimum
  pvals[pvals == -Inf] <- log(.Machine$double.xmin) 
  
  # Negloglik
  tout <- -sum(pvals)
  return(tout)
}

attr(x = cont_r_normal2, which = "required") <- list(sigma_est = list(formula = formula(sigma_est ~ exp(LN_SIGMA), env = globalenv()), 
                                                                      start_values = c(LN_SIGMA = 0)),
                                                     omega_est = list(formula = formula(omega_est ~ exp(LN_OMEGA), env = globalenv()), 
                                                                      start_values = c(OMEGA = 0)))




# Forward gradient function for continuous normal with random intercept in internal_id
grad_cont_r_normal1 <- function (par, dtalist, ret_mat = F, eps = 1e-7, centered = T) {
  
  if(centered) par <- par-eps/2
  # magnitude of change
  # eps <- c(- eps, eps)
  
  # Calculations on unique observations
  d_df <- dtalist$unique
  f_df <- dtalist$full
  nobs <- NROW(d_df)
  fobs <- NROW(f_df)
  
  
  # Set names to parameters
  names(par) <- dtalist$coeff_names
  
  # Fetch formulas
  formulas <- dtalist$formulas
  
  # Fetch parameters and fixed values
  args <- c(as.list(par), as.list(dtalist$fixed))  
  
  
  # Set names to parameters
  plen <- length(par)
  plen1 <- plen +1
  
  parm <- as.data.frame(as.list(par))
  parm <- parm[rep(1, plen+1),]
  
  for(pn in 1:plen) {
    theserows <- pn+1
    parm[theserows, pn] <- parm[theserows, pn]+eps
  }
  
  
  
  
  # Fetch parameters and fixed values
  # CHECK! length of fixed...
  
  
  # Sequence over formulas
  for (formulanum in 1:length(formulas)) {
    formula <- formulas[[formulanum]]
    new_var <- as.character(formula[[2]])
    if (new_var %in% c(colnames(d_df), dtalist$valuevar)) new_var = "Xb"
    # d_df[, new_var] <- with(args, with(d_df, eval(parse(text = formula))))
    tmp <- sapply(1:plen1, function(x) {
      theseargs <- c(as.list(parm[x,]), as.list(dtalist$fixed))
      with(theseargs, with(d_df, eval(parse(text = formula))))
    })
    if(is.null(dim(tmp)) & length(tmp) == plen1) tmp <- matrix(tmp, nrow = nobs, ncol = plen1, byrow = T)
    d_df[, new_var] <- tmp
  }
  
  # matrix of p-calculations
  
  Xb <- d_df$Xb
  
  
  
  betwSD <- exp(parm[['LN_OMEGA']])
  uid <- sort(unique(f_df$internal_id))
  nuid <- NROW(uid)
  
  pvals <- rep(0, nuid)
  
  
  GH_order <- 8L
  if(!is.environment((xregenv <- getOption('xreg')))) {
    options(xreg = (xregenv = new.env(parent = globalenv())))
  }
  if (exists("Gauss_Hermite_Order", envir = xregenv)) {
    GH_order <- get("Gauss_Hermite_Order", envir = xregenv)
  }  else {
    assign(x = "Gauss_Hermite_Order", value = GH_order, envir = xregenv)
  }
  ws <- as.matrix(xreg:::xreg_GHS[[GH_order]]$w)
  xs <- xreg:::xreg_GHS[[GH_order]]$x
  mus <- sqrt(2) * Outer(betwSD, xs)
  # 3D-array, Xb, layers shifted by mus
  # pa <-Xba <- array(Outer(mus, Xb, "+"), c(dim(Xb), GH_order))
  
  Xba <- array(Outer(rep(0, GH_order), Xb, "+"), c(dim(Xb), GH_order))
  musa <- aperm(array(Outer(rep(0, nobs), mus, "+"), c(dim(mus), nobs)))
  pa <- Xba <- Xba+musa
  
  
  sel <- d_df$internal_type == 1 # Discrete
  sel2 <- !sel
  pa[sel,, ] <- dnorm(x = d_df$internal_lb[sel], 
                      mean = Xba[sel,,], 
                      sd = d_df$sigma_est[sel,], 
                      log = T)
  pa[sel2,, ] <- log((pnorm(q = (d_df$internal_lb[sel2]), 
                            mean = Xba[sel2,,], 
                            sd = d_df$sigma_est[sel2,], 
                            lower.tail = F) - 
                        pnorm(q = (d_df$internal_ub[sel2]), 
                              mean = Xba[sel2,,], 
                              sd = d_df$sigma_est[sel2,], 
                              lower.tail = F)))
  
  # Expand to full dataset
  pa <- pa[f_df$umatch,,]
  
  
  # Gauss-Hermite quadrature, calculate loglik contributions per id per variant
  pvals <- apply(X = pa, MARGIN = 2, FUN = function(slice) {
    Log(1/sqrt(pi) * 
          exp(apply(X = slice, 
                    MARGIN = 2, 
                    FUN = group.sum, 
                    ina = f_df$internal_id) + 
                rep_row(Log(ws), nuid)) %*% rep(1, GH_order)) 
  })
  
  
  # Set NA to minimum if boundNA
  if("boundNA" %in% names(dtalist)) if(dtalist$boundNA) pvals[is.na(pvals)] <- log(.Machine$double.xmin)
  
  # Move -Inf to log of minimum
  pvals[pvals == -Inf] <- log(.Machine$double.xmin)
  
  
  
  if(ret_mat) return(pvals)
  pvals <- colSums(pvals)
  
  # Subtract modifications from current, set to initial scale
  pvals <- (pvals[1]-pvals[2:plen1])/eps
  
  
  return(pvals)
}

# Forward gradient function for continuous normal with random intercept in internal_id. 
# Using sapply over each parameter appears to be faster than generating and handling 3D-arrays.
grad_cont_r_normal2 <- function (par, dtalist, ret_mat = F, eps = 1e-7, centered = T) {
  
  if(centered) par <- par-eps/2
  # Calculations on unique observations
  d_df <- dtalist$unique
  f_df <- dtalist$full
  nobs <- NROW(d_df)
  fobs <- NROW(f_df)
  
  
  # Set names to parameters
  names(par) <- dtalist$coeff_names
  
  # Fetch formulas
  formulas <- dtalist$formulas
  
  # Fetch parameters and fixed values
  args <- c(as.list(par), as.list(dtalist$fixed))  
  
  
  # Set names to parameters
  plen <- length(par)
  plen1 <- plen +1
  
  parm <- as.data.frame(as.list(par))
  parm <- parm[rep(1, plen+1),]
  
  for(pn in 1:plen) {
    theserows <- pn+1
    parm[theserows, pn] <- parm[theserows, pn]+eps
  }
  
  
  
  
  # Fetch parameters and fixed values
  # CHECK! length of fixed...
  
  
  # Sequence over formulas
  for (formulanum in 1:length(formulas)) {
    formula <- formulas[[formulanum]]
    new_var <- as.character(formula[[2]])
    if (new_var %in% c(colnames(d_df), dtalist$valuevar)) new_var = "Xb"
    # d_df[, new_var] <- with(args, with(d_df, eval(parse(text = formula))))
    tmp <- sapply(1:plen1, function(x) {
      theseargs <- c(as.list(parm[x,]), as.list(dtalist$fixed))
      with(theseargs, with(d_df, eval(parse(text = formula))))
    })
    if(is.null(dim(tmp)) & length(tmp) == plen1) tmp <- matrix(tmp, nrow = nobs, ncol = plen1, byrow = T)
    d_df[, new_var] <- tmp
  }
  
  # matrix of p-calculations
  uid <- sort(unique(f_df$internal_id))
  nuid <- NROW(uid)
  
  pvals <- rep(0, nuid)
  
  
  GH_order <- 8L
  if(!is.environment((xregenv <- getOption('xreg')))) {
    options(xreg = (xregenv = new.env(parent = globalenv())))
  }
  if (exists("Gauss_Hermite_Order", envir = xregenv)) {
    GH_order <- get("Gauss_Hermite_Order", envir = xregenv)
  }  else {
    assign(x = "Gauss_Hermite_Order", value = GH_order, envir = xregenv)
  }
  ws <- as.matrix(xreg:::xreg_GHS[[GH_order]]$w)
  xs <- xreg:::xreg_GHS[[GH_order]]$x
  
  
  
  
  pvals <- sapply(1:plen1, function(x) {
    Xb <- d_df$Xb[,x]
    
    betwSD <- exp(parm['LN_OMEGA'][x,])
    
    mus <- sqrt(2) * betwSD * xs
    tmpm <- matrix(NA, nrow = nobs, ncol = GH_order)
    Xbm <- Outer(y= Xb, x = mus, "+")
    
    sel <- d_df$internal_type == 1 # Discrete
    sel2 <- !sel
    tmpm[sel, ] <- dnorm(x = d_df$internal_lb[sel], 
                         mean = Xbm[sel,], 
                         sd = d_df$sigma_est[sel,x], 
                         log = T)
    tmpm[sel2, ] <- Log(pnorm(q = (d_df$internal_lb[sel2]), 
                              mean = Xbm[sel2,], 
                              sd = d_df$sigma_est[sel2,x], lower.tail = F) - 
                          pnorm(q = (d_df$internal_ub[sel2]), 
                                mean = Xbm[sel2,], 
                                sd = d_df$sigma_est[sel2,x], 
                                lower.tail = F))
    
    # Expand to full dataset
    tmpm <- tmpm[f_df$umatch,]
    
    pvals <- Log(1/sqrt(pi) * 
                   exp(apply(X = tmpm, 
                             MARGIN = 2, 
                             FUN = group.sum, 
                             ina = f_df$internal_id) + 
                         rep_row(Log(ws), nuid)) %*% rep(1, GH_order)) 
    
    
  })
  
  
  # Set NA to minimum if boundNA
  if("boundNA" %in% names(dtalist)) if(dtalist$boundNA) pvals[is.na(pvals)] <- log(.Machine$double.xmin)
  
  # Move -Inf to log of minimum
  pvals[pvals == -Inf] <- log(.Machine$double.xmin)
  
  
  
  if(ret_mat) return(pvals)
  pvals <- colSums(pvals)
  
  # Subtract modifications from current, set to initial scale
  pvals <- (pvals[1]-pvals[2:plen1])/eps
  
  
  return(pvals)
}

grad_default <- function(func, x, method = "simple",...) numDeriv::grad(func = func, x = x, method = method, ...)



