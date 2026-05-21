# Regression mixture model with fixed-intercept york regression components
# for a single number of components k
# 
# Arguments:
# * data: data frame with 4 columns for x, se(x), y, and se(y), in that order
# * re: scalar, correlation of sampling error between x and y, assumed to be same for all observation pairs
# * k: scalar, number of  mixture components 
# * estimator: which slope estimator to use, one of MLE, alt_MLE, or numeric
# * maxIterations: maximum number of iterations for fitting; applies to both EM and M-step regression iterations
# * EMTolerance: convergence criteria for EM algorithm, in terms of relative change in loglik
# * InnerTolerance: convergence criteria for fitting slope at M step, in terms of relative change in loglik
# * verbose: whether to print loglik at each iteration; setting to >1 additional prints ll for inner slope iterations
#
# Returns:
# structured list of class 'magicmap'


magicmap_1k <- function(data, re, k, estimator="MLE", maxIterations=1000, EMTolerance=1e-12, InnerTolerance=1e-12, verbose=FALSE){
  
  # sanity checks
  if(class(data)[1] != "data.frame"){
    stop("data must be an object of class data.frame, try converting data to a data.frame \'data = as.data.frame(data)\'")
  }
  
  if(!(estimator %in% c("MLE","alt_MLE","numeric"))){
    stop("estimator must be one of \'MLE\', \'alt_MLE\', or \'numeric\'.")
  }
  
  
  dat <- data
  colnames(dat) <- c("x","sx","y","sy")
  dat$re <- rep(re,nrow(dat))
  
  ##################  
  
  # single component
  if(k==1){
    york <- suppressWarnings(yorktools::york(dat,intercept=0,method = "ML",tol=InnerTolerance,maxit=maxIterations,verbose=verbose))
    b_est <- york$b[1]
    p_est <- 1
    post_est <- matrix(1, nrow(dat), 1)
    x_est <- yorktools::pred.york(dat, model=york, intercept=0)$xpred
    ll_comp <- ll_comp_func(dat, b_est, p_est, x_est)
    em_converge <- york$converge$converged
    empty_comp <- F
    iter <- NA
    em_dec <- NA
    york_dec <- NA
    
  }else{
  
    # initial groups to get inital parameters
    ratios <- dat$y/dat$x
    ratio_var <- dat$sy^2/dat$x^2 + 
                    dat$y^2*dat$sx^2/dat$x^4 -
                    2*re*dat$y*dat$sy*dat$sx/dat$x^3
    init_clust <- Ckmeans.1d.dp::Ckmeans.1d.dp(ratios,k=k,y=1/ratio_var)$cluster
    
    # initial parameters
    b_init <- rep(NA,k)
    p_init <- rep(NA,k)
    for(i in 1:k){
      b_init[i] <- stats::lm(dat$y[init_clust==i] ~ 0 + dat$x[init_clust==i], 
                             weights=1/((dat$sx[init_clust==i]^2)+(dat$sy[init_clust==i]^2))
                             )$coefficients[1]
      p_init[i] <- mean(init_clust==i)
    }
    xpred_init <- dat$x
    
    
    b_est <- b_init
    p_est <- p_init
    x_est <- xpred_init
    ll_comp <- ll_comp_func(dat, b_est, p_est, x_est)
  
    # convergence and warning flags
    em_converge <- F
    # =========================================================================
    # BEGIN SQUAREM INTEGRATION
    # =========================================================================
    
    # 1. Define the Objective Function (Calculates the Marginal Log-Likelihood)
    em_objective <- function(theta, data, k_comps, ...) {
      b_curr <- theta[1:k_comps]
      p_curr <- theta[(k_comps+1):(2*k_comps)]
      x_curr <- theta[(2*k_comps+1):length(theta)]
      return(ll_comp_func(data, b_curr, p_curr, x_curr))
    }
    
    # 2. Define the Fixed-Point Function
    em_fixed_point <- function(theta, data, k_comps, est_method, inn_tol, max_it, ...) {
      
      # Unpack parameters
      b_curr <- theta[1:k_comps]
      p_curr <- theta[(k_comps+1):(2*k_comps)]
      x_curr <- theta[(2*k_comps+1):length(theta)]
      
      # --- E STEP ---
      post_curr <- estep_posteriors(data, b_curr, p_curr, x_curr)
      
      # --- M STEP ---
      p_new <- colMeans(post_curr)
      
      # Prevent hard zeroes which crash the likelihood
      if(any(p_new == 0)) {
        p_new <- pmax(p_new, 1e-10)
        p_new <- p_new / sum(p_new)
      }
      
      b_new <- b_curr
      x_new <- x_curr
      
      # Execute your existing M-Step estimators (The York Regression)
      if(est_method == "numeric") {
        par <- optim(c(b_curr, x_curr), 
                     fn=mstep_optim_func, gr=mstep_gradient_func,
                     data=data, priorp=pmax(p_new, 1e-16), postp=post_curr, method="BFGS",
                     control=list(fnscale=-1, maxit=max_it, reltol=inn_tol))
        b_new <- par$par[1:k_comps]
        x_new <- par$par[(k_comps+1):(nrow(data)+k_comps)]
        
      } else {
        # Internal York iterations
        ll_Q_y_old <- -Inf
        york_converge_inner <- F
        y_iter <- 0
        
        while(!york_converge_inner && y_iter < max_it) {
          y_iter <- y_iter + 1
          
          if(est_method == "MLE") {
            b_new <- slope_est_func(data, b_new, x_new, post_curr)
          } else if(est_method == "alt_MLE") {
            b_new <- slope_est_func2(data, b_new, post_curr)
          }
          
          x_new <- xpred_est_func(data, b_new, post_curr)
          ll_Q_y_new <- sum(Q_loglik_matrix(data, b_new, p_new, x_new, post_curr))
          
          # Check inner York convergence
          if(abs((ll_Q_y_new - ll_Q_y_old) / ll_Q_y_new) < inn_tol) {
            york_converge_inner <- T
          } else {
            ll_Q_y_old <- ll_Q_y_new
          }
        }
      }
      
      # Repack and return updated parameters to SQUAREM
      return(c(b_new, p_new, x_new))
    }
    
    # 3. Bundle initial parameters and run SQUAREM accelerator
    theta_init <- c(b_est, p_est, x_est)
    
    sq_fit <- SQUAREM::squarem(
      par = theta_init, 
      fixptfn = em_fixed_point, 
      objfn = em_objective,
      control = list(tol = EMTolerance, maxiter = maxIterations, trace = verbose, minimize = FALSE),
      # Arguments passed down to the fixpt and obj functions:
      data = dat, k_comps = k, est_method = estimator, inn_tol = InnerTolerance, max_it = maxIterations
    )
    
    # 4. Unpack the final, accelerated parameters back into your script's variables
    b_est <- sq_fit$par[1:k]
    p_est <- sq_fit$par[(k+1):(2*k)]
    x_est <- sq_fit$par[(2*k+1):length(sq_fit$par)]
    
    # Run the E-step one final time with the optimized parameters so post_est is ready for output
    post_est <- estep_posteriors(dat, b_est, p_est, x_est)
    ll_comp <- sq_fit$value.objfn
    iter <- sq_fit$fpevals # Number of times the fixed-point function (EM cycle) was evaluated
    
    # Check if SQUAREM hit the tolerance or maxed out iterations
    # Set convergence flags expected by the rest of the script
    empty_comp <- any(p_est < 1e-8)
    em_converge <- sq_fit$convergence
    em_dec <- NA     
    york_dec <- NA
    
    # =========================================================================
    # END SQUAREM INTEGRATION
    # =========================================================================
  }
  # model df for k-1 priors and k slopes
  df <- 2*k-1
  
  # fit stats
  BIC <- df*log(nrow(dat)) - 2*ll_comp
  AIC <- 2*df - 2*ll_comp
    
  # build output
  # order components by slope
  if(k>1){
    idx <- order(b_est, decreasing=T)
    b_est <- b_est[idx]
    p_est <- p_est[idx]
    post_est <- post_est[,idx]
  }
  
  colnames(post_est) <- paste0("ProbClass",1:k)
  rownames(post_est) <- rownames(data)
  posteriors <- cbind(data,post_est)
    
  # Genimi code, passed initial inspection
  # =================================================================
  # BEGIN PER-COMPONENT METRICS (D2)
  # =================================================================
  raw_mahalanobis_d2 <- vector("list", k)

  if (!empty_comp) {
    # 2. Pre-calculate the inverse covariance matrix elements for ALL SNPs
    # This vectorizes the math so we don't need slow loops!
    det_Sigma <- dat$sx^2 * dat$sy^2 - (dat$re * dat$sx * dat$sy)^2
    inv_11 <- dat$sy^2 / det_Sigma
    inv_22 <- dat$sx^2 / det_Sigma
    inv_12 <- -(dat$re * dat$sx * dat$sy) / det_Sigma

    for(c in 1:k) {
      # 3. Calculate D^2 for ALL variants relative to this component's specific line
      resid_x <- dat$x - x_est
      resid_y <- dat$y - (b_est[c] * x_est)

      # The Mahalanobis quadratic form
      d2_all <- (resid_x^2 * inv_11) + (resid_y^2 * inv_22) + (2 * resid_x * resid_y * inv_12)
      
      # Save the raw vector into the list
      raw_mahalanobis_d2[[c]] <- d2_all
    }
  }

  # Build the comprehensive parameters dataframe with placeholders for the post-hoc SEs
  parameters <- data.frame(
    b = b_est, 
    se = NA
  )
  
  # Attach the raw D2 list as a list-column (data.frames don't like lists in the init function)
  parameters$raw_mahalanobis_d2 <- raw_mahalanobis_d2
  
  # End Gemini code
  # =================================================================
  
  notes <- ""
  
  if(any(abs(diff(b_est,lag=1)) < 0.01)){
    warning(paste("Slopes appear collapsed to fewer than ",k,"components"))
    notes <- paste0(notes, "Slopes not well differentiated; ")
    em_converge <- NA
  }
  if(!is.na(em_dec)){
    notes <- paste0(notes, paste0("EM LL decreased ",signif(em_dec,4),"; "))
  }
  if(!is.na(york_dec)){
    notes <- paste0(notes, paste0("Inner LL decreased, max step ",signif(york_dec,4),"; "))
  }
  if(ll_comp == -Inf){
    em_converge <- NA
    notes <- paste0(notes, paste0("LL not finite; "))
  }
  
  
  out <- list(
    fit_stats = data.frame(k_components=k, 
                           logLik=ll_comp, 
                           model_df=df, 
                           BIC=BIC, 
                           AIC=AIC,
                           converged=em_converge,
                           notes=notes,
                           iter=iter),
    component_proportions = p_est,
    slopes = parameters,
    predicted_target_betas = x_est,
    posteriors = posteriors
  )
  
  if(!is.na(em_converge) && !em_converge){
    warning(paste("EM algorithm for",k,"components failed to converge in",iter,"iterations"))
  }
  
  if(empty_comp){
    warning(paste("EM algorithm for",k,"components encountered empty component"))
    em_converge <- FALSE
    
    if(notes==""){
      notes <- paste("EM halted due to empty component")
    }else{
      notes <- paste0(paste("EM halted due to empty component"),"; ",notes)
    }
    
    out <- list(
      fit_stats = data.frame(k_components=k, 
                             logLik=NA, 
                             model_df=df, 
                             BIC=NA, 
                             AIC=NA,
                             converged=F,
                             notes=notes,
                             iter=iter),
      component_proportions = NA,
      slopes = NA,
      predicted_target_betas = NA,
      posteriors = NA
    )
  }
  
  class(out) <- 'magicmap'
  
  return(out)
  
}


#' Fit MAGICMAP model
#'
#' @description
#' Fits MAGICMAP (Mixtures Aggregating Genes Into Coordinated Modules And Pathways), a regression mixture model where each component is a york regression with the intercept fixed to zero.
#'
#' @details
#' For a full description of the model, see the vignette pdf.
#' 
#' @param data
#' Data frame with the specified column names containing genetic effect sizes and standard errors
#' 
#' @param betaTargetX 
#' String, giving the column name in \code{data} for the effect sizes (betas) on the target (x-axis) trait
#' 
#' @param sdTargetX, 
#' String, giving the column name in \code{data} for the standard error (se) of the effect sizes on the target (x-axis) trait
#' 
#' @param betaComparatorY
#' String, giving the column name in \code{data} for the effect sizes (betas) on the comparator (y-axis) trait

#' @param sdComparatorY,
#' String, giving the column name in \code{data} for the standard error (se) of the effect sizes on the comparator (y-axis) trait
#' 
#' @param k
#' Number of components to fit. Can be a scalar or a vector.
#' 
#' @param ids 
#' (Optional) string giving the column name in \code{data} of variant IDs to use for labeling results tables
#' 
#' @param CovIntercept 
#' Intercept from LDSC (Bulik-Sullivan et al. 2015) genetic correlation analysis (\code{gcov_int}) comparing the two traits. Default = 0, which assumes the two sets of effect sizes do not have correlated error (i.e. no sample overlap).
#' 
#' @param TargetXIntercept
#' Intercept from LDSC heritability analysis of the target (x axis) trait. Default = 1, which assumes no population stratification or other confounding.
#' 
#' @param ComparatorYIntercept
#' Intercept from LDSC heritability analysis of the comparator (y axis) trait. Default = 1, which assumes no population stratification or other confounding.
#' 
#' @param inflateSEs
#' Logical, whether to apply LSDC h2 intercept adjustment to the input standard errors. By default, input SEs are used as given.
#' 
#' @param estimator
#' Which estimator to use for the regression slopes. One of "MLE", "alt_MLE", or "numeric".
#' 
#' @param maxIterations
#' Maximum number of iterations for fitting each mixture model. Applies separately to EM iterations and slope fitting iterations in each M step. Default = 1000.
#' 
#' @param EMTolerance
#' Tolerance to use for EM convergence, in terms of relative change in the log likelihood. Default = 1e-12.
#' 
#' @param InnerTolerance
#' Tolerance to use for convergence of slope estimator at M step, in terms of relative change in the log likelihood. Default = 1e-12.
#' 
#' @param verbose 
#' Whether to print verbose logging. 0=no logging, 1=logging of EM iterations, 2=logging of EM and M step iterations.
#' 
#'
#' @return
#' A list of class 'magicmap' containing:
#' 
#' \describe{
#'
#'   \item{\code{scoutjoy_test}}{
#'     Global test of heterogeneity from SCOUTJOY (Elliott et al., 2024). P-value is for null hypothesis that the effect sizes have a single homogeneous relationship for all variants.
#'   }
#'   
#'   \item{\code{fit_stats}}{
#'     Table of fit statistics for the mixture model with each requested value of \code{k}. Includes:
#'     \itemize{
#'       \item \code{logLik}: log of the complete data (marginal) likelihood
#'       \item \code{model_df}: the number of free paramters in the model, excluding nuisance parameters
#'       \item \code{BIC}: Bayesian information criterion; lower values indicate a better combination of fit and complexity
#'       \item \code{AIC}: Akaike information criterion; lower values indicate a better combination of fit and complexity
#'       \item \code{converged}: boolean indicator of whether the model fitting converged; set to NA if fitting encountered issues other than convergence (see notes field)
#'       \item \code{notes}: text description of additional considerations
#'       \item \code{iter}: number of EM iterations required for model convergence
#'     }
#'   }
#'   
#'   \item{\code{component_proportions}}{
#'     Estimated prior probabilities for membership in each mixture component.
#'   }
#'   
#'   \item{\code{slopes}}{
#'     Estimated slope in each mixture component. SEs are currently always NA.
#'   }
#'   
#'   \item{\code{predicted_target_betas}}{
#'     Estimated values of the genetic effects on the target (x axis) trait.
#'   }
#'   
#'   \item{\code{posteriors}}{
#'     Data frame containing the original input data (i.e. the X effect size, X standard error, Y effect size, and Y standard error) and the estimated posterior probability of the observations being a member of each of the \code{k} mixture components.
#'   }
#'   
#'
#' }
#' 
#' The \code{component_proportions}, \code{slopes}, \code{predicted_target_betas}, and \code{posteriors} are each structured as a list with entries named for the number of fitted components (e.g. \code{mix2components}, \code{mix3components}, etc). Within each model, mixture components are ordered from highest to lowest slope.
#' 
#' 
#' @references
#' \itemize{
#'   \item Elliott, A. \emph{et al}. Distinct and shared genetic architectures of gestational diabetes mellitus and type 2 diabetes. \emph{Nat. Genet.} 56(3), 377-382 (2024). PMCID: PMC10937370.
#'   \item Bulik-Sullivan, B. \emph{et al}. An atlas of genetic correlations across human diseases and traits. \emph{Nat. Genet.} 47(11): 1236-1241 (2015). PMCID: PMC4797329.
#' }
#'
#' @export 
#'

magicmap <- function(data, betaTargetX, sdTargetX, betaComparatorY, sdComparatorY, k, ids = NULL, CovIntercept = 0, TargetXIntercept=1, ComparatorYIntercept=1, inflateSEs=FALSE, estimator="MLE", maxIterations=1000, EMTolerance=1e-12, InnerTolerance=1e-12, verbose=FALSE){

  # sanity checks
  if((length(betaTargetX)!=1) || (length(sdTargetX)!=1) || (length(betaComparatorY)!=1) || (length(sdComparatorY)!=1)){
    stop("betaTargetX, sdTargetX, betaComparatorY, and sdComparatorY must each be a single column name")
  }
  
  if(class(data)[1] != "data.frame"){
    stop("data must be an object of class data.frame, try converting data to a data.frame \'data = as.data.frame(data)\'")
  }
  
  if(!(estimator %in% c("MLE","alt_MLE","numeric"))){
    stop("estimator must be one of \'MLE\', \'alt_MLE\', or \'numeric\'.")
  }
  
  
  # format data
  if(!is.null(ids)){
    row.names(data) <- as.character(data[,ids])
  }  
  data <- data[, c(betaTargetX, sdTargetX, betaComparatorY, sdComparatorY)]
  if(any(is.na(data))){
    na_rows <- apply(data,1,function(a) any(is.na(a)))
    warning(paste("Removing",na_rows,"rows with missing values"))
    data <- data[!na_rows, ]
  }
  
  data[, c(betaComparatorY, betaTargetX)] <- data[, c(betaComparatorY, betaTargetX)] * sign(data[, betaTargetX[1]])
  
  # handle covariance
  re <- CovIntercept / sqrt(TargetXIntercept*ComparatorYIntercept)
  
  if(inflateSEs){
    data[,sdTargetX] <- data[,sdTargetX]*sqrt(TargetXIntercept)
    data[,sdComparatorY] <- data[,sdComparatorY]*sqrt(ComparatorYIntercept)
  }
  
  scout <- SCOUTJOY::scoutjoy(betaComparatorY, betaTargetX, sdComparatorY, sdTargetX, data, CovIntercept = CovIntercept, ExposureIntercept = TargetXIntercept, OutcomeIntercept=ComparatorYIntercept, printProgress = FALSE)
  

  fits <- NULL
  slopes <- list()
  priors <- list()
  xpreds <- list()
  posts <- list()
  stop_k_too_high <- F
  
  for(i in seq_along(k)){
    if(verbose){
      print(paste0("Testing k = ",k[i],"...")) 
    }

    mod <- magicmap_1k(data=data, 
                       re=re, 
                       k=k[i], 
                       estimator=estimator, 
                       maxIterations=maxIterations, 
                       EMTolerance=EMTolerance, 
                       InnerTolerance=InnerTolerance, 
                       verbose=verbose)
    
    fits <- rbind(fits, mod$fit_stats)
    priors[[i]] <- mod$component_proportions
    slopes[[i]] <- mod$slopes
    xpreds[[i]] <- mod$predicted_target_betas
    posts[[i]] <- mod$posteriors
    
    if(any(grepl("Slopes not well differentiated",mod$fit_stats$notes))){
      stop_k_too_high <- T
      if(i < length(k)){
        k <- k[1:i]
        warning(paste0("Components collapsed at k=",k[i],"; skipping remaining k values"))
      }
      break
    }
    
        
  }

  names(priors) <- paste0("mix",k,"components")
  names(slopes) <- paste0("mix",k,"components")
  names(xpreds) <- paste0("mix",k,"components")
  names(posts) <- paste0("mix",k,"components")

  # Gemini code, past initial inspection. I haven't verified the math yet
  # =================================================================
  # POST-HOC STANDARD ERROR CALCULATION (Lowest BIC Model Only)
  # =================================================================
    # 1. Identify the "winning" model based on lowest BIC
    best_i <- which.min(fits$BIC)
    best_k <- fits$k_components[best_i]
    print(paste0("Calculating SE for the model with ", best_k, " components"))
    
    # Only calculate if the winning model actually converged
    if (!is.na(fits$converged[best_i]) && fits$converged[best_i] && best_k > 0) {
      
      # 2. Extract the converged parameters for this specific model
      b_best <- slopes[[best_i]]$b
      x_best <- xpreds[[best_i]]
      p_best <- priors[[best_i]]
      # Extract just the probability columns and convert to matrix
      post_best <- as.matrix(posts[[best_i]][, -c(1:4)]) 
      
      # Recreate the formatted data for the gradient function
      dat_hess <- data
      colnames(dat_hess) <- c("x","sx","y","sy")
      dat_hess$re <- rep(re, nrow(dat_hess))
      
      # 3. Calculate the Hessian
      theta_hat <- c(b_best, x_best)
      n_obs <- nrow(dat_hess)
      
      H_start_time <- Sys.time()
      H_full <- numDeriv::jacobian(func = mstep_gradient_func, 
                                   x = theta_hat, 
                                   data = dat_hess, 
                                   priorp = pmax(p_best, 1e-16), 
                                   postp = post_best)
      
      # 4. Schur Complement
      H_bb <- H_full[1:best_k, 1:best_k, drop=FALSE]
      H_xx_diag <- diag(H_full[(best_k+1):(n_obs+best_k), (best_k+1):(n_obs+best_k)])
      H_bx <- H_full[1:best_k, (best_k+1):(n_obs+best_k), drop=FALSE]
      H_xb <- H_full[(best_k+1):(n_obs+best_k), 1:best_k, drop=FALSE]
      
      penalty_term <- H_bx %*% (H_xb / H_xx_diag)
      info_matrix_b <- -H_bb + penalty_term
      
      # 5. Invert and assign
      tryCatch({
        cov_b <- solve(info_matrix_b)
        se_best <- sqrt(diag(cov_b))
        
        # Overwrite the NAs in the slopes list for the winning model!
        slopes[[best_i]]$se <- se_best          
      }, error = function(e) {
        # Printing the actual error message 'e' helps debug typos vs. actual singular matrices!
        warning(paste("Error calculating SE for best model (k=", best_k, "): ", e$message, sep=""))
      })
      H_end_time <- Sys.time()
      print(paste0('Calculation of the variance-covariance matrix took ',as.numeric(H_end_time - H_start_time, units = "secs"), ' secs.'))
    }
  # End Gemini code
  # =================================================================
  
  out <- list(
    scoutjoy_test=scout$`Global Test`,
    fit_stats = fits,
    component_proportions = priors,
    slopes = slopes,
    predicted_target_betas = xpreds,
    posteriors = posts
  )
  
  class(out) <- 'magicmap'
  
  if(scout$`Global Test`$Pvalue > .05){
    warning("!! Scoutjoy global test is not significant, cannot reject 1 component solution !!")
    out$fit_stats$notes <- paste("SCOUTJOY does not reject k=1;",out$fit_stats$notes)
  }
  
  return(out)
  
}
  



#' Plot MAGICMAP results
#'
#' 
#' @param model
#' output of \code{\link{magicmap}}()
#' 
#' @param which_plot
#' numeric, which model to plot (when multiple models were fit). Default is the model with lowest BIC.
#' 
#' @param class_thresh
#' threshold of posterior probability to use for assigning variants to components
#' 
#' @param label_comp
#' Numeric, which mixture component to label points from. Unclassified points are treated as the k+1 component.
#' 
#' @param colors
#' set of k+1 colors to use for plotting. Last color is used for unclassified variants.
#' 
#' @param se_bars
#' logical, whether to include standard error bars
#' 
#' @param legend
#' logical, whether to include a legend
#' 
#' @param comp_names
#' list of names to use for mixture components (only affects legend)
#' 
#' @param se_length
#' length for caps on SE intervals; passed to \code{\link[graphics]{arrows}}()
#'
#' @param ...
#' additional parameters to be passed directly to \code{\link[base]{plot}}()
#'
#'
#' @export plot.magicmap
#' @export
#'

plot.magicmap <- function(model, which_plot=NULL, class_thresh=0.95, label_comp=NULL, colors=NULL, se_bars=TRUE, legend=TRUE, comp_names=NULL, se_length=0.025, conf_region=TRUE, ...){  
  args <- list(...)
  
  #conf_region=FALSE # plotting code available, but not currently part of magicmap estimation
  
  if(nrow(model$fit_stats)==1){
    
    if(is.na(model$fit_stats$converged) | !model$fit_stats$converged){
      warning("plotted model failed to converge")
    }
    
    x <- model$posteriors[[1]][,1]
    y <- model$posteriors[[1]][,3]
    sx <- model$posteriors[[1]][,2]
    sy <- model$posteriors[[1]][,4]
    
    if(!("xlab" %in% names(args))){
      args$xlab <- colnames(model$posteriors[[1]])[1]
    }
    
    if(!("ylab" %in% names(args))){
      args$ylab <- colnames(model$posteriors[[1]])[3]
    }
    
    k <- model$fit_stats$k[1]
    
    if(model$scoutjoy_test$Pvalue > .05 && k!=1){
      warning(paste0("plotting k=",k," but SCOUTJOY did not reject k=1"))
    }
    
    slopes <- model$slopes[[1]][,"b"]
    slopese <- model$slopes[[1]][,"se"]
    
    if(k==1){
      classes=rep(1, nrow(model$posteriors[[1]]))
    }else{
      classes <- apply(model$posteriors[[1]][,-c(1:4)], MARGIN = 1, which.max)
      maxprob <- apply(model$posteriors[[1]][,-c(1:4)], MARGIN = 1, max)
      classes[maxprob < 0.95] <- k+1
    }
    
  }else{
    
    if(is.null(which_plot)){
      which_plot <- which.min(model$fit_stats$BIC)
    }
    
    if(is.na(model$fit_stats$converged[which_plot])){
      warning("model selected for plotting had issues with fitting; check notes in fit_stats table")
    }else if(!model$fit_stats$converged[which_plot]){
      warning("model selected for plotting failed to converge")
    }else if(any(is.na(model$fit_stats$converged)) | any(!model$fit_stats$converged) ){
      warning("model selected for plotting was fit successfully, but some comparison models did not")
    }
    
    k <- model$fit_stats$k[which_plot]
    
    if(model$scoutjoy_test$Pvalue > .05 && k!=1){
      warning(paste0("plotting k=",k," but SCOUTJOY did not reject k=1"))
    }
    
    x <- model$posteriors[[paste0("mix",k,"components")]][,1]
    y <- model$posteriors[[paste0("mix",k,"components")]][,3]
    sx <- model$posteriors[[paste0("mix",k,"components")]][,2]
    sy <- model$posteriors[[paste0("mix",k,"components")]][,4]
    
    if(!("xlab" %in% names(args))){
      args$xlab <- colnames(model$posteriors[[paste0("mix",k,"components")]])[1]
    }
    
    if(!("ylab" %in% names(args))){
      args$ylab <- colnames(model$posteriors[[paste0("mix",k,"components")]])[3]
    }
    
    slopes <- model$slopes[[paste0("mix",k,"components")]][,"b"]
    slopese <- model$slopes[[paste0("mix",k,"components")]][,"se"]
    
    if(k==1){
      classes=rep(1, nrow(model$posteriors[[paste0("mix",k,"components")]]))
    }else{
      classes <- apply(model$posteriors[[paste0("mix",k,"components")]][,-c(1:4)], MARGIN = 1, which.max)
      maxprob <- apply(model$posteriors[[paste0("mix",k,"components")]][,-c(1:4)], MARGIN = 1, max)
      classes[maxprob < class_thresh] <- k+1
    }
  }
  
  if(is.null(colors)){
    if(k > 7){
      stop(paste("Selected model has",k,"components, but there are only default colors for 7. Please specify colors."))
    }
    # colorbrewer dark2
    colors <- c("#1b9e77","#d95f02","#7570b3","#e7298a","#66a61e","#e6ab02","#a6761d")
    colors <- c(colors[1:k],"#666666")
  }
  
  
  if(!("xlim" %in% names(args))){
    args$xlim <- c(min(0,min(x-sx)), 1.1*max(x+sx))
  }
  
  if(!("ylim" %in% names(args))){
    args$ylim <- c(min(0,min(y-sy)), 1.1*max(y+sy))
  }
  
  
  if(!("main" %in% names(args))){
    args$main <- ""
  }
  
  if(!("pch" %in% names(args))){
    args$pch <- 20
  }
  
  if(!("cex" %in% names(args))){
    args$cex <- 0.8
  }
  
  if(!("col" %in% names(args))){
    args$col <- colors[classes]
  }
  
  if(!("lwd" %in% names(args))){
    args$lwd <- 2
  }
  
  if(!("las" %in% names(args))){
    args$las <- 1
  }
  
  if(!("mar" %in% names(args))){
    mar <- c(3.5, 3.5, 1.5, 1.5)
  }else{
    mar <- args$mar
  }
  
  if(!("mgp" %in% names(args))){
    mgp <- c(2.25,1,0)
  }else{
    mgp <- args$mgp
  }
  
  mar_bak <- par()$mar
  mgp_bak <- par()$mgp
  
  par(mar=mar, mgp=mgp)
  
  args$x <- 0
  args$y <- 0
  plotcol <- args$col
  args$col <- rgb(0,0,0,0)
  args$bty='l'
  do.call(plot, args)
  args$col <- plotcol
  
  if(conf_region){
    for(i in 1:k){
      
      # ---> SAFETY CHECK <---
      if(!is.na(slopese[i])) {
        
        if(min(x-sx)<0){
          polygon(x=c(
            0,
            min(x-sx), 
            min(x-sx),
            0,
            max(x+sx), 
            max(x+sx)
          ),
          y=c(
            0,
            min(x-sx)*(slopes[i]-qnorm(.025,lower=F)*slopese[i]), 
            min(x-sx)*(slopes[i]+qnorm(.025,lower=F)*slopese[i]), 
            0,
            max(x+sx)*(slopes[i]+qnorm(.025,lower=F)*slopese[i]), 
            max(x+sx)*(slopes[i]-qnorm(.025,lower=F)*slopese[i])
          ),
          border=FALSE,
          col=adjustcolor(colors[i], alpha.f=0.3))
        }else{
          polygon(x=c(
            min(x-sx), 
            min(x-sx),
            max(x+sx), 
            max(x+sx)
          ),
          y=c(
            min(x-sx)*(slopes[i]-qnorm(.025,lower=F)*slopese[i]), 
            min(x-sx)*(slopes[i]+qnorm(.025,lower=F)*slopese[i]), 
            max(x+sx)*(slopes[i]+qnorm(.025,lower=F)*slopese[i]), 
            max(x+sx)*(slopes[i]-qnorm(.025,lower=F)*slopese[i])
          ),
          border=FALSE,
          col=adjustcolor(colors[i], alpha.f=0.3))
        }
        
      } # ---> CLOSE THE SAFETY CHECK <---
    }
  }
  
  abline(h=0,col="grey20",lwd=0.5)
  abline(v=0,col="grey20",lwd=0.5)
  abline(0,1,col="grey80",lty=2)
  abline(0,-1,col="grey80",lty=2)
  
  
  if(se_bars){    
    arrows(x0=x, x1=x, y0=y-sy, y1=y+sy, length=se_length, angle=90, code=3, col="grey50")
    arrows(x0=x-sx, x1=x+sx, y0=y, y1=y, length=se_length, angle=90, code=3, col="grey50")
  }
  
  points(x, y,
         pch=args$pch,
         col=args$col,
         cex=args$cex)
  
  for(i in 1:k){
    abline(0, slopes[i], col=colors[i], lwd=args$lwd)
  }
  
  
  if(legend){
    if(is.null(comp_names)){
      comp_names <- c(paste("Component ",1:k), "Unclassified")
    }else if(length(comp_names)==k){
      comp_names <- c(comp_names, "Unclassified")
    }
    
    legend("topright", 
           fill=colors,
           border=colors,
           # lwd=rep(lwd,k+1), 
           cex = 0.9*args$cex, 
           bg="white",
           box.col="white", 
           legend = comp_names)
  }
  
  if(!is.null(label_comp)){
      
    if(nrow(model$fit_stats)==1){
      df_post <- model$posteriors
    }else{
      df_post <- model$posteriors[[paste0("mix",k,"components")]]
    }
    
    for(c in label_comp){
      
      text(x=df_post[classes==c, 1], 
           y=df_post[classes==c, 3], 
           labels=rownames(df_post)[classes==c], 
           cex=0.8, adj=c(-.1,-.25))
    }
  }

  par(mar=mar_bak, mgp=mgp_bak)
  
}

# detach("package:MAGICMAP", unload = T)
# remove.packages("MAGICMAP")
# setwd("~/Documents/Code/github/MAGICMAP")
# devtools::document()
# devtools::build(manual=T)
# install.packages("~/Documents/Code/github/MAGICMAP_0.1.1.1.tar.gz")
# require(MAGICMAP)
# mm <- MAGICMAP::magicmap(dat, "x","sx","y","sy",k=1:3,estimator="MLE")
# mm <- MAGICMAP::magicmap(data.frame(x=x_obs,y=y,sx=x_se,sy=x_se), "x","sx","y","sy",k=1:3, estimator="MLE")
