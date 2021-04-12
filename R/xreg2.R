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
                   latent_classes = 0, 
                   latent_class_parameters = character(), 
                   latent_id_colname = character(), 
                   return_type = "fit", 
                   print_sum = F,
                   method = "ucminf",
                   hessian = T,
                   ...) {
  xreg2_time <- Sys.time()
  if("xreg2_obj" %in% class(controlList)) return(xreg2Optim(xreg2obj = controlList, method = method, hessian = hessian, ...))
  if(is.null(dataList)) stop('Argument "dataList" is missing, with no default.')
  
  dot_args <- list(...)
  dot_args[c("p_fun", "p_aggrgation_fun")] <- NULL
  if (latent_classes < 2) 
    latent_class_parameters <- character()
  if (!is.character(latent_class_parameters)) {
    warning("Latent_class_parameters not character. Latent classes will be disregarded.")
    latent_class_parameters <- character()
    latent_classes <- 0
    latent_id_colname <- character()
  }
  if (!is.character(latent_id_colname)) {
    warning("latent_id_colname not character. Latent classes will be disregarded.")
    latent_class_parameters <- character()
    latent_classes <- 0
    latent_id_colname <- character()
  }
  if (length(latent_id_colname)) 
    latent_id_colname <- latent_id_colname[1]
  lik_fun = xlik
  startValues <- attr(controlList, "startValues")
  if ("xreg" %in% class(controlList)) {
    if (is.na(start_values[1])) 
      start_values <- controlList$coef
    controlList <- controlList$controlList
  }
  formulas <- controlList
  controlList <- do.call(xregControl, match_formals(fun = xregControl, 
                                                    formulas = controlList))
  if ("xreg" %in% class(start_values)) 
    startValues <- start_values$coef
  if ("numeric" %in% class(start_values) & length(start_values) > 0) 
    startValues <- start_values
  if ("data.frame" %in% class(dataList)) {
    dataList <- list(dataList)
    names(dataList)[1] <- names(controlList)[1]
    warning(paste0("Dataframe provided where list expected in dataList argument. The dataframe has been enclosed in a list, and named \"", 
                   names(controlList)[1], "\" to match the first xregControl-object."))
  }
  if (!all(names(dataList) %in% names(controlList))) 
    stop(paste("No control list provided for ", paste(names(dataList)[!names(dataList) %in% 
                                                                        names(controlList)], collapse = ", ")))
  if (is.null(names(dataList))) {
    if (length(dataList) <= length(controlList)) {
      warning("Provided datalist was unnamed. Names coerced from controlList")
      names(dataList) <- names(controlList)[1:length(dataList)]
    }
    else {
      stop("Provided datalist was unnamed, and longer than provided controllist. Exiting.")
    }
  }
  if (length(unique(names(startValues))) < length(startValues)) {
    warning(paste0("More than one start value assigned for parameters (", 
                   paste(names(startValues)[-match(unique(names(startValues)), 
                                                   names(startValues))], collapse = ", "), "). First value used."))
    startValues <- startValues[unique(names(startValues))]
  }
  if (length(fixed_values)) {
    if ("xreg" %in% class(fixed_values)) {
      fixed_values <- fixed_values$pars
    }
    if ("data.frame" %in% class(fixed_values)) {
      fixed_df <- fixed_values[, 1:2]
      fixed_values <- fixed_values[, 1]
      names(fixed_values) <- rownames(fixed_df)
    }
    else {
      fixed_df <- data.frame(Estimate = fixed_values, 
                             `Std. Error` = rep(NA, NROW(fixed_values)))
      colnames(fixed_df) <- c("Estimate", "Std. Error")
      rownames(fixed_df) <- names(fixed_values)
    }
  }
  else {
    fixed_df <- data.frame(a = numeric(), b = numeric())
  }
  colnames(fixed_df) <- c("Estimate", "Std. Error")
  fixedValues <- numeric()
  this_call <- match.call()
  initControlList <- controlList
  aggr <- TRUE
  return_lik_df <- (return_type == "df")
  return_first <- (return_type == "first")
  if (return_type == "predict") {
    aggr <- FALSE
    return_lik_df = TRUE
    return_type = "first_df"
  }
  if (return_type == "first_df") {
    return_first <- TRUE
    return_first_df <- TRUE
  }
  else {
    return_first_df <- FALSE
  }
  if (return_type == "precalc_df") {
    return_first <- TRUE
    return_lik_df <- FALSE
    return_first_df <- TRUE
  }
  start_vals = vector()
  prep_out <- list()
  newDataList <- list()
  all_defined_vars <- character()
  undefined_vars <- character()
  valuevars <- list()
  for (dataName in names(dataList)) {
    rel_col <- vector()
    thisControl <- controlList[[dataName]]
    data_df <- dataList[[dataName]]
    orig_cols <- colnames(data_df)
    valueVar <- NA
    defined_vars <- colnames(data_df)
    local_vars <- character()
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
    # message(paste(colnames(data_df), collapse = ", "))
    controlList[[dataName]]$upper_bound_var <- upper_bound_var
    controlList[[dataName]]$lower_bound_var <- lower_bound_var
    controlList[[dataName]]$valueVar <- valueVar
    data_df[, c("internal_ub", "internal_lb")] <- data_df[, c(upper_bound_var, lower_bound_var)]
    formula[[2]] <- substitute(Xb)
    e_params <- list(...)
    if (length(latent_id_colname)) {
      if (latent_id_colname %in% colnames(data_df)) 
        data_df$internal_classgroup <- data_df[, latent_id_colname]
    }
    else data_df$internal_classgroup <- paste0(dataName, ".", 1:NROW(data_df))
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
      internal_type[internal_ub == Inf] <- 0
      internal_type[internal_lb == -Inf] <- 2
      internal_count <- 1
    })
    if (sum(data_df$both_na)) {
      warning(paste0(sum(data_df$both_na), " rows in which ", 
                     valueVar, " was NA were removed."))
      data_df <- data_df[!data_df$both_na, -(NCOL(data_df))]
    }
    if (!is.na((weights_var <- thisControl$weights_var)) & 
        weights_var %in% colnames(data_df)) 
      data_df$internal_count <- data_df[, weights_var]
    relevant_columns <- c(rel_col, colnames(data_df)[grep("internal_", 
                                                          colnames(data_df))])
    if (aggr) {
      new_df <- data_df[, colnames(data_df) %in% relevant_columns]
      new_df$internal_unique <- as.numeric(as.factor(apply(new_df, 
                                                           MARGIN = 1, FUN = "paste", collapse = "")))
      unique_df <- unique(new_df)
      unique_df$internal_count <- aggregate(new_df$internal_count, 
                                            by = list(new_df$internal_unique), FUN = "sum")[unique_df$internal_unique, 
                                                                                            2]
      newDataList[[dataName]] <- unique_df
    }
    else {
      newDataList[[dataName]] <- data_df
    }
    controlList[[dataName]]$orig_cols <- orig_cols
    startValues <- c(startValues, thisControl$start_values[!names(thisControl$start_values) %in% 
                                                             names(startValues)])
    fixedValues <- c(fixedValues, thisControl$fixed_values[!names(thisControl$fixed_values) %in% 
                                                             names(fixedValues)])
    all_defined_vars <- c(all_defined_vars, defined_vars[!defined_vars %in% 
                                                           all_defined_vars])
    controlList[[dataName]]$obs_types <- c(Uncensored = sum(data_df$internal_count * 
                                                              (data_df$internal_ub == data_df$internal_lb)), `Left censored` = sum(data_df$internal_count * 
                                                                                                                                     (data_df$internal_ub > data_df$internal_lb) * (data_df$internal_ub < 
                                                                                                                                                                                      Inf) * (data_df$internal_lb == -Inf)), `Right censored` = sum(data_df$internal_count * 
                                                                                                                                                                                                                                                      (data_df$internal_ub > data_df$internal_lb) * (data_df$internal_lb > 
                                                                                                                                                                                                                                                                                                       -Inf) * (data_df$internal_ub == Inf)), Intervals = sum(data_df$internal_count * 
                                                                                                                                                                                                                                                                                                                                                                (data_df$internal_ub > data_df$internal_lb) * (data_df$internal_lb > 
                                                                                                                                                                                                                                                                                                                                                                                                                 -Inf) * (data_df$internal_ub < Inf)))
  }
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
                                            c(names(startValues), latent_class_parameters)])
  undefined_vars <- undefined_vars[!undefined_vars %in% c("Xb", 
                                                          names(fixedValues))]
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
  if (latent_classes > 1) {
    for (lc in latent_class_parameters) {
      tmp <- 0.1
      names(tmp) <- lc
      if (lc %in% names(startValues)) {
        tmp <- startValues[lc]
        startValues <- startValues[-which(names(startValues) == 
                                            lc)]
        all_defined_vars <- all_defined_vars[-which(all_defined_vars == 
                                                      lc)]
      }
      for (i in 1:latent_classes) {
        lci <- paste(lc, i, sep = ".")
        if (!lci %in% names(startValues)) {
          startValues <- c(startValues, tmp)
          names(startValues)[length(startValues)] <- lci
        }
        if (!lci %in% all_defined_vars) {
          all_defined_vars <- c(all_defined_vars, lci)
        }
      }
    }
  }
  if (length((tmp <- names(startValues[!names(startValues) %in% all_defined_vars])))) {
    startValues <- startValues[!names(startValues) %in% tmp]
    tmp <- paste(tmp, collapse = ", ")
    warning(paste0("Start values defined for variables (", 
                   tmp, ") not found in column names or formulas. These have been removed."))
  }
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
  
  
  
  startValues[names(fixed_values)] <- fixed_values
  
  relcols <- relevant_columns[!relevant_columns %in% names(startValues)]
  
  dtns <- names(newDataList)
  names(dtns) <- dtns
  xreg2lobj <- lapply(X = dtns, FUN = function(dtn) {
    thisfull <- newDataList[[dtn]]
    thisfull <- thisfull[, colnames(thisfull)[!colnames(thisfull) %in% c("internal_classgroup")]]
    tmpfull <- thisfull[, colnames(thisfull)[!colnames(thisfull) %in% c("internal_id", "internal_count", "internal_unique")]]
    thisunique <- unique(tmpfull)
    # Add rows in unique data.frame to full df
    thisfull$umatch <- multimatch(tmpfull, thisunique)
    startValues <- startValues[!names(startValues) %in% names(fixed_values)]
    
    greek <- c("ALPHA", "BETA", "GAMMA", "DELTA", "EPSILON", "ZETA", "ETA", "THETA", "IOTA", "KAPPA", "LAMBDA", "MU", "NU", "XI", "OMIKRON", "PI", "RHO", "SIGMA", "TAU", "UPSILON", "PHI", "CHI", "PSI", "OMEGA")
    
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
  if(return_type == "control") {
    return(xreg2lobj)
  } else {
    
    oout <- xreg2Optim(xreg2lobj, method = method, hessian = hessian, ...)
  }
  
  oout$xreg2_time <- Sys.time()-xreg2_time
  
  return(oout)
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  if (return_first) {
    ret <- do.call(lik_fun, c(startValues, list(controlList = controlList, 
                                                dataList = newDataList, print_sum = F, return_lik_df = return_lik_df, 
                                                return_first_df = return_first_df, latent_classes = latent_classes, 
                                                latent_class_parameters = latent_class_parameters, 
                                                dot_args)))
    for (nm in names(ret)) {
    }
    return(ret)
  }
  testmle <- do.call(x_mle, c(list(minuslogl = lik_fun, start = as.list(startValues), 
                                   controlList = controlList, dataList = newDataList, control = control, 
                                   print_sum = print_sum, fixed = as.list(fixed_values), 
                                   latent_classes = latent_classes, latent_class_parameters = latent_class_parameters), 
                              dot_args))
  dfs <- do.call(lik_fun, c(c(testmle$coef, fixed_values), 
                            list(controlList = controlList, dataList = newDataList, 
                                 print_sum = F, return_lik_df = T, return_first_df = T, 
                                 latent_classes = latent_classes, latent_class_parameters = latent_class_parameters, 
                                 dot_args)))
  pars <- cbind(testmle$fullcoef, type = rep("Fitted", NROW(testmle$fullcoef)))
  testmle$fixed_values <- fixed_df
  testmle$fixed_values <- fixed_df
  pars <- rbind(pars, cbind(fixed_df, type = rep("Fixed", 
                                                 NROW(fixed_df))))
  res_types <- list()
  pars[, paste0("n_", names(dfs))] <- 0
  obs_types <- as.data.frame(matrix(rep(0, length(names(dfs)) * 
                                          4), nrow = 4))
  colnames(obs_types) <- names(dfs)
  rownames(obs_types) <- c("Uncensored", "Left-censored", 
                           "Right-censored", "Intervals")
  for (nm in names(dfs)) {
    res_types[["p_sum"]][nm] <- sum(controlList[[dataName]]$p_aggregation_fun(dfs[[nm]]))
    pars[rownames(pars) %in% controlList[[nm]]$defined_vars, 
         paste0("n_", nm)] <- res_types[["counts"]][nm] <- sum(dfs[[nm]]$internal_count)
    obs_types[, nm] <- controlList[[nm]]$obs_types
  }
  obs_types$total <- rowSums(obs_types)
  pars[, "n_sum"] <- rowSums(pars[, paste0("n_", names(dfs)), 
                                  F])
  pars[, "t value"] <- abs(pars$Estimate)/pars$"Std. Error"
  pars[, "Pr(>|t|)"] <- 2 * pt(pars$"t value", df = pars$n_sum - 
                                 1, lower.tail = F)
  testmle$pars <- pars
  testmle$pars_v <- pars[, 1]
  testmle$coef <- testmle$pars[, 1]
  names(testmle$coef) <- rownames(testmle$pars)
  res_types[["counts"]][["total_count"]] <- sum(res_types[["counts"]])
  res_types[["p_sum"]]["total"] <- sum(res_types[["p_sum"]])
  testmle$minima <- res_types$p_sum
  testmle$obs_types <- obs_types
  testmle$info <- rbind(obs_types, total_obs = res_types$counts, 
                        minima = res_types$p_sum, logLik = -res_types$p_sum)
  names(testmle$pars_v) <- rownames(pars)
  res <- list(controlList = controlList, start_values = start_values, 
              mle_obj = testmle, coef = testmle$coef, coefficients = testmle$coef, 
              full_coef = testmle$pars, valueVar = valueVar, fixed_values = fixed_df, 
              minima = res_types[["p_sum"]], counts = res_types[["counts"]], 
              obs_types = obs_types, pars = pars, pars_v = testmle$pars_v)
  class(res) <- c("xreg", "list")
  return(res)
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




