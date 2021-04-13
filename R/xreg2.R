#returns match vector for two dataframes with the same number of columns
multimatch <- function(df1, df2) {
  if(ncol(df1) != ncol(df2)) stop("Data.frames must have same number of columns.")
  a <- apply(df1, 1, function(x) paste(x, collapse = "_"))
  b <- apply(df2, 1, function(x) paste(x, collapse = "_"))
  
  return(match(a, b))
}

xregControl2 <- function (formulas, 
                          start_values = numeric(), 
                          fixed_values = numeric(), 
                          p_fun = cont_normal2,
                          grad_fun = function(par, ...) numDeriv::grad(func = p_fun, x = par, method = "simple", ...),
                          p_aggregation_fun = function(d_df) return(-d_df$p * d_df$internal_count), 
                          weights_var = NA, 
                          name = NA, 
                          censor_bounds = c(-Inf, Inf), 
                          lower = NA, 
                          upper = NA) 
{
  if (class(formulas)[1] %in% c("xregControl", "xregControlList")) 
    return(c(formulas))
  
  if (class(formulas)[1] == "formula") 
    formulas <- list(formulas)
  required <- attr(x = p_fun,which = "required")
  
  defined_vars <- sapply(formulas, function(x){as.character(x[[2]])})
  if (length(required)) {
    
    req <- as.list(rep(0, length(required)))
    names(req) <- names(required)
    # print(req)
    # print(required)
    for (i in 1:length(required)) {
      if (!names(req)[i] %in% defined_vars) {
        
        defined_vars <- c(defined_vars, names(req)[i])
        
        req_return <- required[[i]]
        # print(req_return)
        formulas <- c(formulas, list(req_return$formula))
        if (exists("start_values", req_return))
          start_values <- c(start_values, req_return$start_values[!req_return$start_values %in% start_values])
        if (exists("fixed_values", req_return))
          fixed_values <- c(fixed_values, req_return$fixed_values[!req_return$fixed_values %in% fixed_values])
      }
    }
  }
  
  
  start_values <- start_values[unique(names(start_values))]
  retObj <- mget(c(names(formals()), "defined_vars"))
  retObj$columnNames <- character()
  for (formula in formulas) {
    retObj$columnNames <- c(retObj$columnNames, all.vars(formula)[!all.vars(formula) %in% 
                                                                    retObj$columnNames])
  }
  retObj$columnNames <- retObj$columnNames[!retObj$columnNames %in% 
                                             names(start_values)]
  class(retObj) <- c("xregControl", "list")
  retObj <- c(retObj)
  return(retObj)
}



c.xregControl <- function(...) {
  args <- list(...)
  
  if(is.null(names(args) )) names(args) <- NA
  for(i in 1:length(args)) {
    if(is.na(names(args)[i])) {
      if(is.na(args[[i]]$name)) {
        names(args)[i] <- i
      } else {
        names(args)[i] <- args[[i]]$name
      }
    }
    if(!"xregControl" %in% class(args[[i]])) stop("At least one non-xregControl object among arguments.")
  }
  class(args) <- c("xregControlList", "list")
  
  
  columnNames <- character()
  startValues <- numeric()
  fixedValues <- numeric()
  lower <- numeric()
  upper <- numeric()
  for(xregControl in args) {
    columnNames <- c(columnNames, xregControl$columnNames[!xregControl$columnNames %in% columnNames])
    startValues <- c(startValues, xregControl$start_values[!names(xregControl$start_values) %in% names(startValues)])
    fixedValues <- c(fixedValues, xregControl$fixed_values[!names(xregControl$fixed_values) %in% names(fixedValues)])
    lower <- c(lower, xregControl$lower[!xregControl$lower %in% lower])
    upper <- c(upper, xregControl$upper[!xregControl$upper %in% upper])
  }
  
  
  attr(args, which = "columnNames") <- columnNames
  attr(args, which = "startValues") <- startValues
  attr(args, which = "fixedValues") <- fixedValues
  attr(args, which = "lower") <- lower
  attr(args, which = "upper") <- upper
  
  
  return(args)
}

#' @export
c.xregControlList <- function(...) {  #return(do.call(c.xregControl, do.call(c.list, list(...))))
  args <- list(...)
  for(narg in 1:length(args)) {
    class(args[[narg]]) <- "list"
  }
  #return(args)
  return(do.call(c.xregControl, do.call(c, args)))
}




xreg2 <- function (controlList, 
                     dataList = NULL, 
                     start_values = numeric(), 
                     fixed_values = numeric(), 
                     return_type = "fit", 
                     print_sum = F,
                     method = "ucminf",
                     hessian = T,
                     bounds = F,
                     ...) {
  
  xreg2_time <- Sys.time()
  if("xreg2_obj" %in% class(controlList)) return(xreg2Optim(xreg2obj = controlList, method = method, hessian = hessian, ...))
  if(is.null(dataList)) stop('Argument "dataList" is missing, with no default.')
  
  
  dot_args <- list(...)
  dot_args[c("p_fun", "p_aggrgation_fun")] <- NULL
  this_call <- match.call()

  
  #####################
  #### PREPARATION ####
  #####################
  
  startValues <- attr(controlList, "startValues") # if controlList is an xregControl Object it will find start values
  fixedValues <- numeric()
  
  
  # find start values and formula in xreg object entered as controlList(formula)
  if ("xreg" %in% class(controlList)) {
    if (is.na(start_values[1])) 
      start_values <- controlList$coef
    controlList <- controlList$controlList 
    # is the formula always saved as controlList in an xreg object?
    # we need the formula here saved as controlList because next lines work with controlList as a formula
  }
  

  formulas <- controlList
  controlList <- do.call(xregControl2, match_formals(fun = xregControl2, 
                                                     formulas = controlList)) # ll functions need to be loaded
  
  
  ### START VALUES ####
  
  if ("xreg" %in% class(start_values)) 
    startValues <- start_values$coef
  if ("numeric" %in% class(start_values) & length(start_values) > 0) 
    startValues <- start_values
  # use only first given value for each parameter
  if(length(unique(names(startValues))) < length(startValues)) {
    warning(paste0("More than one start value assigned for parameters (", 
                   paste(names(startValues)[-match(unique(names(startValues)), 
                                                   names(startValues))], collapse = ", "), "). First value used."))
    startValues <- startValues[unique(names(startValues))]
  }
  
  
  #### FIXED VALUES ####
  
  if (length(fixed_values)) {
    if ("xreg" %in% class(fixed_values)) {
      fixed_values <- fixed_values$pars
    }
    if ("data.frame" %in% class(fixed_values)) {
      fixed_df <- fixed_values[, 1:2] # this is not dynamically
      fixed_values <- fixed_values[, 1] # # this is not dynamically
      names(fixed_values) <- rownames(fixed_df)
    }else {
      fixed_df <- data.frame(Estimate = fixed_values, 
                             `Std. Error` = rep(NA, NROW(fixed_values)))
    }
  }else {
    fixed_df <- data.frame(Estimate = numeric(), `Std. Error` = numeric())
  }
  
  
  colnames(fixed_df) <- c("Estimate", "Std. Error")  # can be left out, if `Std Error` works above
  rownames(fixed_df) <- names(fixed_values)
  
  
  
  #### TRANSFORM DATALIST TO A LIST ####
  
  if ("data.frame" %in% class(dataList)) {
    dataList <- list(dataList)
    names(dataList)[1] <- names(controlList)[1]
    warning(paste0("Dataframe provided where list expected in dataList argument. The dataframe has been enclosed in a list, and named \"", 
                   names(controlList)[1], "\" to match the first xregControl-object."))
  }
  
  if (!all(names(dataList) %in% names(controlList))) 
    stop(paste("No control list provided for ", paste(names(dataList)[!names(dataList) %in% names(controlList)], collapse = ", ")))
  if (is.null(names(dataList))) {
    if (length(dataList) <= length(controlList)) {
      warning("Provided datalist was unnamed. Names coerced from controlList")
      names(dataList) <- names(controlList)[1:length(dataList)]
    }
    else {
      stop("Provided datalist was unnamed, and longer than provided controllist. Exiting.")
    }
  }
  
  
  
  ### WHAT DOES THIS DO ?? ####
  #  aggr <- TRUE
  return_lik_df <- (return_type == "df")
  return_first <- (return_type == "first")
  if (return_type == "predict") {
    #    aggr <- FALSE
    return_lik_df = TRUE
    return_type = "first_df"
  }
  if (return_type == "first_df") {
    return_first <- TRUE
    return_first_df <- TRUE
  }else {
    return_first_df <- FALSE
  }
  if (return_type == "precalc_df") {
    return_first <- TRUE
    return_lik_df <- FALSE
    return_first_df <- TRUE
  }
  
  
  
  # Preparing Loop of dataList elements to find relevant columns
  # and start values
  newDataList <- list()
  all_defined_vars <- character()
  undefined_vars <- character()
  valuevars <- list()
  
  
  ### LOOP OVER DATA LIST ###
  
  for (dataName in names(dataList)) {
    rel_col <- vector()
    thisControl <- controlList[[dataName]]
    data_df <- dataList[[dataName]]
    # orig_cols <- colnames(data_df)
    valueVar <- NA
    defined_vars <- colnames(data_df)
    local_vars <- character()
    
    
    # loop over different formulas to select target variable
    # and relevant and local variables: not everything calculated here is really needed later
    # maybe do this with lapply??
    for (fnum in 1:length(thisControl$formulas)) {
      formula <- thisControl$formulas[[fnum]]
      targetName <- as.character(formula[[2]])
      rel_col <- c(rel_col, all.vars(formula)[!all.vars(formula) %in% 
                                                rel_col])
      if (is.na(valueVar)) 
        if (targetName %in% colnames(data_df) | 
            (length(grep(paste0(targetName, "\\."), colnames(data_df), value = TRUE)) == 2)) {
          valuevars[[dataName]] <- targetName
          valueVar <- targetName
          
        }
      rhs_vars <- all.vars(formula[[3]])
      undefined_vars <- c(undefined_vars, rhs_vars[!rhs_vars %in% 
                                                     defined_vars])
      all_defined_vars <- c(all_defined_vars, undefined_vars[!undefined_vars %in% 
                                                               all_defined_vars])
      if (!targetName %in% c(colnames(data_df), valueVar)) 
        defined_vars <- c(defined_vars, targetName)
      local_vars <- c(local_vars, all.vars(formula), targetName)
    }
    
    local_vars <- local_vars[!local_vars %in% colnames(data_df)]
    controlList[[dataName]][["defined_vars"]] <- local_vars
    
    
    ### BOUNDS FOR CENSORING ####
    
    if(any(thisControl$censor_bounds != c(-Inf,Inf))){
      censor_bounds <- thisControl$censor_bounds
      grep_vals <- grep(paste0(valueVar, "\\."), colnames(data_df), value = TRUE)
      if (length(grep_vals) > 2) 
        stop("More than two columns matching target variable.")
      if (length(grep_vals) == 2) {
        if (all((data_df[, grep_vals[1]] - data_df[, grep_vals[2]]) >= 0)) {
          upper_bound_var <- grep_vals[1]
          lower_bound_var <- grep_vals[2]
        }
        else if (all((data_df[, grep_vals[1]] - data_df[, grep_vals[2]]) <= 0)) {
          upper_bound_var <- grep_vals[2]
          lower_bound_var <- grep_vals[1]
        }
        else stop(paste0("Two columns matched the left-hand side: ", 
                         grep_vals[1], " and ", grep_vals[2], ". However, neither was allways >= the other, so they cannot be used as upper/lower bounds. Quitting."))
      }
      else {
        if (!valueVar %in% colnames(data_df)) {
          stop("No value columns specified. ")
        }
        else {
          upper_bound_var <- lower_bound_var <- valueVar
        }
      }
      
      controlList[[dataName]]$upper_bound_var <- upper_bound_var
      controlList[[dataName]]$lower_bound_var <- lower_bound_var
      controlList[[dataName]]$valueVar <- valueVar
      data_df[, c("internal_ub", "internal_lb")] <- data_df[, c(upper_bound_var, lower_bound_var)]
      
      formula[[2]] <- substitute(Xb) # WHY???
      e_params <- list(...)
      data_df$internal_classgroup <- paste0(dataName, ".", 1:NROW(data_df))
      
      if("id" %in% colnames(data_df)) data_df$internal_id <- as.integer(data_df$id)
      data_df <- within(data_df, {
        both_na <- is.na(internal_ub) * is.na(internal_lb)
        internal_ub[internal_ub >= max(thisControl$censor_bounds)] <- Inf
        internal_lb[internal_lb <= min(thisControl$censor_bounds)] <- -Inf
        internal_lb[internal_lb >= max(thisControl$censor_bounds)] <- max(thisControl$censor_bounds)
        internal_ub[internal_ub <= min(thisControl$censor_bounds)] <- min(thisControl$censor_bounds)
        internal_ub[is.na(internal_ub)] <- Inf
        internal_lb[is.na(internal_lb)] <- -Inf
        internal_type <- rep(3, NROW(internal_ub))
        internal_type[internal_ub == internal_lb] <- 1
        internal_type[internal_ub == Inf] <- 0 # should this be more dynamically?
        internal_type[internal_lb == -Inf] <- 2 # should this be more dynamically?
        internal_count <- 1
      })
      if (sum(data_df$both_na)) {
        warning(paste0(sum(data_df$both_na), " rows in which ", valueVar, " was NA were removed."))
        data_df <- data_df[!data_df$both_na, -(NCOL(data_df))]
      }
      
    }else{ # if we don?t have bounds, we don´t need to calcute everything above
      formula[[2]] <- substitute(Xb) # WHY???
      e_params <- list(...)
      
      controlList[[dataName]]$upper_bound_var <- valueVar
      controlList[[dataName]]$lower_bound_var <- valueVar
      controlList[[dataName]]$valueVar <- valueVar
      data_df[, c("internal_ub", "internal_lb")] <- data_df[, c(valueVar, valueVar)]
      data_df[,"internal_type"] <- 1
      data_df[,"internal_count"] <- 1
      data_df$both_na <- is.na(data_df$internal_ub) * is.na(data_df$internal_lb)
      if(sum(data_df$both_na)){
        warning(paste0(sum(data_df$both_na), " rows in which ", valueVar, " was NA were removed."))
        data_df <- data_df[!data_df$both_na, -(NCOL(data_df))]
      }
    }
    
    #### WEIGHTS ####
    if (!is.na((weights_var <- thisControl$weights_var)) & weights_var %in% colnames(data_df)) 
      data_df$internal_count <- data_df[, weights_var]
    
    relevant_columns <- c(rel_col, colnames(data_df)[grep("internal_",colnames(data_df))])
    
    # I Think we don?t need this part any longer
    # while creating xreg2lobject we save the unique dataset
    
    #   if (aggr) {
    #      new_df <- data_df[, is.element(colnames(data_df),relevant_columns)]
    #      new_df$internal_unique <- as.numeric(as.factor(apply(new_df,MARGIN = 1, FUN = "paste", collapse = "")))
    #      unique_df <- unique(new_df)
    #      unique_df$internal_count  <- new_df %>% group_by(internal_unique) %>% summarise(n = n()) %>% .$`n`
    #      newDataList[[dataName]] <- unique_df
    #    }
    #    else {
    newDataList[[dataName]] <- data_df[, is.element(colnames(data_df),relevant_columns)]
    #    }
  
    
    
    # Maybe this can be changed? 
    # does it rellay have to be in the loop or can this be found in the xregControl Object directly?
    startValues <- c(startValues, thisControl$start_values[!is.element(names(thisControl$start_values), names(startValues))])
    
    fixedValues <- c(fixedValues, thisControl$fixed_values[!is.element(names(thisControl$fixed_values),names(fixedValues))])
    
    all_defined_vars <- c(all_defined_vars, defined_vars[!is.element(defined_vars,all_defined_vars)])
    
    controlList[[dataName]]$obs_types <- c(Uncensored = sum(data_df$internal_count *(data_df$internal_ub == data_df$internal_lb)), 
                                           `Left censored` = sum(data_df$internal_count * (data_df$internal_ub > data_df$internal_lb) * (data_df$internal_ub < Inf) * (data_df$internal_lb == -Inf)), 
                                           `Right censored` = sum(data_df$internal_count * (data_df$internal_ub > data_df$internal_lb) * (data_df$internal_lb > -Inf) * (data_df$internal_ub == Inf)),
                                           Intervals = sum(data_df$internal_count * (data_df$internal_ub > data_df$internal_lb) * (data_df$internal_lb >-Inf) * (data_df$internal_ub < Inf)))
  } # end of loop over dataList
  ### END LOOP DATA LIST ###
  
  
  ### PREPARE START AND FIXED VALUES ###
  
  if (length(fixedValues)) {
    fixed_values[names(fixedValues)] <- fixedValues
    tmp <- rep(NA, NROW(fixed_values))
    names(tmp) <- names(fixed_values)
    tmp[rownames(fixed_df)] <- fixed_df[, 2]
    fixed_df <- data.frame(Estimate = fixed_values, `Std. Error` = tmp)
    colnames(fixed_df) <- c("Estimate", "Std. Error")
  }
  
  
  startValues[names(fixed_values)[names(fixed_values) %in% 
                                    names(startValues)]] <- fixed_values[names(fixed_values) %in% 
                                                                           names(startValues)]
  undefined_vars <- unique(undefined_vars[!undefined_vars %in% 
                                            names(startValues)])
  undefined_vars <- undefined_vars[!undefined_vars %in% c("Xb",names(fixedValues))]
  
  
  
  ### SETTING START VALUES WHICH ARE NOT GIVEN ####
  
  if (length(undefined_vars)) {
    tmp <- paste(undefined_vars, collapse = ", ")
    used_from <- "xreg"
    if ("run_from" %in% names(dot_args)) 
      used_from <- dot_args[["run_from"]]
    if (!used_from == "hyreg") 
      warning(paste0("Used formulas contained variables not found in start_values or column names (", 
                     tmp, "). These will be initiated with a starting value of 0.1. Defining start values either in xregControl or in the call to xreg is highly recommended."))
    tmp <- rep(0.1, length(undefined_vars))
    names(tmp) <- undefined_vars
    startValues <- c(startValues, tmp)
  }
  
  if (length((tmp <- names(startValues[!names(startValues) %in% all_defined_vars])))) {
    startValues <- startValues[!names(startValues) %in% tmp]
    tmp <- paste(tmp, collapse = ", ")
    warning(paste0("Start values defined for variables (", 
                   tmp, ") not found in column names or formulas. These have been removed."))
  }
  
  startValues[names(fixed_values)] <- fixed_values
  relcols <- relevant_columns[!is.element(relevant_columns, names(startValues))]
  
  
  ### E_PARAMS ###
  # What is this for?
  if ("optim_control" %in% names(e_params)) {
    ctrl <- e_params[["optim_control"]]
    dot_args["optim_control"] <- NULL
  }
  
  ctrl <- switch(method,
                 "BFGS" = list(maxit = 10000, 
                               abstol = 0.000000001, 
                               reltol = 0.00000000001),
                 "L-BFGS-B" = list(maxit = 10000, 
                                   factr = 0.00000000001, 
                                   lmm = 10))
  
  
  ###########################################
  ### CREATE XREG2LOBJ AND ESTIMATE MODEL ###
  ###########################################
  
  dtns <- names(newDataList)
  names(dtns) <- dtns
  
  greek <- c("ALPHA", "BETA", "GAMMA", "DELTA", "EPSILON", "ZETA", "ETA", "THETA", "IOTA", "KAPPA", "LAMBDA", "MU", "NU", "XI", "OMIKRON", "PI", "RHO", "SIGMA", "TAU", "UPSILON", "PHI", "CHI", "PSI", "OMEGA")
  
  # change it to lapply over the list itself?
  xreg2lobj <- lapply(X = dtns, FUN = function(dtn) {
    thisfull <- newDataList[[dtn]]
    thisfull <- thisfull[, colnames(thisfull)[!colnames(thisfull) %in% c("internal_classgroup")]]
    tmpfull <- thisfull[, colnames(thisfull)[!colnames(thisfull) %in% c("internal_id", "internal_count", "internal_unique")]]
    thisunique <- unique(tmpfull)
    # Add rows in unique data.frame to full df
    thisfull$umatch <- multimatch(tmpfull, thisunique)
    startValues <- startValues[!names(startValues) %in% names(fixed_values)]
    
    
    
    # Special value names
    spvns <- unique(unlist(lapply(c(greek, "INTERCEPT", "SLOPE"), function(x) grep(pattern = x, x = names(startValues), value = T))))
    
    # Special values
    spvs <- startValues[spvns]
    
    # Start values minus special values
    startValues <- startValues[!names(startValues) %in% spvns]
    
    # # Start values sorted
    # startValues <- startValues[sort(names(startValues))]
    
    # Start values with special values at the end
    startValues <- c(startValues, spvs)
    
    
    return(list(full = thisfull,
                unique = thisunique,
                coeff_names = names(startValues),
                formulas = controlList[[dtn]]$formulas,
                random_formula = NULL,
                fixed = fixed_values,
                ll_fun = controlList[[dtn]]$p_fun,
                grad_fun = controlList[[dtn]]$grad_fun,
                valuevar = valuevars[[dtn]],
                startValues = startValues,
                ctrl = ctrl))
    
  })
  class(xreg2lobj) <- c("xreg2_obj", "list")
  # message(xreg2lobj[[1]]$ctrl)
  
  ### ESTIMATE MODEL ###
  
  if(return_type == "control") {
    return(xreg2lobj)
  } else {
    
    oout <- xreg2Optim(xreg2lobj, method = method, hessian = hessian, ...)
  }
  
  oout$xreg2_time <- Sys.time()-xreg2_time
  
  return(oout)
}




xreg2Optim <- function(xreg2obj, method = NULL, hessian = T, ...) {
  # ctrl <- list(maxit = 10000, 
  #              abstol = 0.000000001, 
  #              reltol = 0.00000000001)
  
  optim_time_start <- Sys.time()
  
  if(tolower(method) == "ucminf") {
    oout <- ucminf::ucminf(par = xreg2obj[[1]]$startValues, fn = negll_list2, gr = gradll_list2, dtalist = xreg2obj, hessian = hessian, control = xreg2obj[[1]]$ctrl,...)
  } else {
    oout <- optimr(par = xreg2obj[[1]]$startValues, fn = negll_list2, gr = gradll_list2, method = method, dtalist = xreg2obj, hessian = hessian, control = xreg2obj[[1]]$ctrl,...)  
  }
  
  
  oout$SE <- rep(NA, length(oout$par))
  names(oout$SE) <- names(oout$par) <- xreg2obj[[1]]$coeff_names
  oout$fixed <- xreg2obj[[1]]$fixed
  
  if(hessian) {
    # oout$hessian <- numDeriv::hessian(func = negll_list2, x = oout$par, dtalist = xreg2obj)
    oout$vcov <-     solve(oout$hessian)
    oout$SE[] <- sqrt(diag(oout$vcov))
  }
  
  coef <- oout$par
  
  ncoef <- length(coef)
  nfixed <- length(oout$fixed)
  min <- oout$value
  
  
  oout$coef <- coef
  
  oout$coef_wfixed <- fullcoef <- data.frame(Estimate = c(coef, oout$fixed))
  if (NROW(fullcoef)) {
    fullcoef[, "Std. Error"] <- NA
    fullcoef$type <- rep(c("Fitted", "Fixed"), c(ncoef, nfixed))
    if (hessian) 
      fullcoef[1:ncoef, "Std. Error"] <- oout$SE
  }
  
  oout$full_coef <- fullcoef
  
  oout$method = method
  oout$optim_time <- Sys.time()-optim_time_start
  oout
}




