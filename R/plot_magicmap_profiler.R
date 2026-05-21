#' Generate MAGICMAP Diagnostic Profiler Plot
#'
#' Evaluates model convergence, calculates comprehensive diagnostic metrics for each
#' component, generates a dynamic diagnostic PNG plot, and returns the updated model.
#' 
#' @param model The fitted MAGICMAP model object
#' @param prefix String, the base filename used to name the output PNG
#' @return The MAGICMAP model object, updated with profiler_metrics_raw and profiler_logic_dictionary
#'
#' @export

generate_profiler_plot <- function(model, prefix, gene_label = TRUE) {
  
  # =========================================================================
  # 1. MODEL SELECTION & VALIDATION
  # =========================================================================
  best_row_idx <- which.min(model$fit_stats$BIC)
  n_fitted_comp <- model$fit_stats$k_components[best_row_idx]
  fit_converged <- model$fit_stats$converged[best_row_idx]
  comp_key <- paste0("mix", n_fitted_comp, "components")
  
  cat(sprintf("  -> Profiler evaluating best model: k=%d (Converged: %s)\n", n_fitted_comp, fit_converged))
  
  if (!isTRUE(fit_converged)) {
    cat(sprintf("[WARNING] Model collapsed or failed to fit. Skipping diagnostics and plotting.\n"))
    return(model) # Return original model unmodified
  }
  
  best_model_params <- model$slopes[[comp_key]]
  estimated_slopes <- best_model_params[, "b"]
  
  if (length(estimated_slopes) == 0) {
    cat("No estimated slopes found. Skipping plotting.\n")
    return(model)
  }
  
  # =========================================================================
  # 2. EXTRACT DATA & POSTERIORS
  # =========================================================================
  # MAGICMAP posteriors inherently store the original data in columns 1-4 (x, sx, y, sy)
  full_post_df <- model$posteriors[[comp_key]]
  x_vals <- full_post_df[, 1]
  sx_vals <- full_post_df[, 2]
  
  post_data <- as.matrix(full_post_df[, -c(1:4), drop=FALSE])
  classes <- apply(post_data, MARGIN = 1, which.max)
  raw_d2_list <- best_model_params$raw_mahalanobis_d2
  
  # =========================================================================
  # 3. DIAGNOSTIC LOGIC & HELPER FUNCTION
  # =========================================================================
  diagnostic_logic_dict <- list(
    "Insufficient N"     = "n_confident < 4",
    "Low N"              = "n_confident < 7",
    "No Spread"          = "!is.na(x_spread) && x_spread < 0.05",
    "Leveraged"          = "!is.na(lev_entropy) && lev_entropy < 0.40",
    "Somewhat Leveraged" = "!is.na(lev_entropy) && lev_entropy > 0.40 && lev_entropy < 0.60",
    "Weak Inst."         = "!is.na(mean_f) && mean_f < 50.0",
    "Fuzzy"              = "cred_ratio < 0.40",
    "Somewhat Fuzzy"     = "cred_ratio >= 0.40 && cred_ratio <= 0.60",
    "High mean D2"       = "!is.na(m1) && m1 > 2.0",
    "Overfit"            = "!is.na(ad_p) && ad_p < 0.05 && !is.na(m1) && m1 < 1.0",
    "Loose Fit"          = "!is.na(ad_p) && ad_p < 0.05 && (is.na(m1) || m1 >= 1.0)",
    "WtD2-Reject"        = "!is.na(wt_pval) && wt_pval < 0.05"
  )
  
  get_cluster_diagnostics <- function(w_i, assigned_idx, raw_d2, x_vals, sx_vals) {
    n_assigned <- length(assigned_idx)
    cred_ratio <- 0; cont_cred <- 0; x_spread <- NA; lev_entropy <- NA
    m1 <- NA; m2 <- NA; ad_p <- NA; wt_d2 <- NA; wt_pval <- NA
    df_eff <- NA; mean_f <- NA
    
    n_confident <- sum(w_i > 0.95)
    eff_n <- sum(w_i)
    
    if (eff_n > 0) {
      cred_ratio <- sum(w_i[w_i > 0.95]) / eff_n
      valid_idx <- is.finite(raw_d2) & is.finite(w_i) & is.finite(x_vals) & is.finite(sx_vals)
      if (sum(valid_idx) > 0) {
        w_valid <- w_i[valid_idx]; d2_valid <- raw_d2[valid_idx]
        x_valid <- x_vals[valid_idx]; sx_valid <- sx_vals[valid_idx]
        eff_n_valid <- sum(w_valid)
        
        wt_d2 <- sum(w_valid * d2_valid)
        df_eff <- max(1e-5, eff_n_valid - 1)
        wt_pval <- pchisq(wt_d2, df = df_eff, lower.tail = FALSE)
        mean_f <- sum(w_valid * (x_valid / sx_valid)^2) / eff_n_valid
        cont_cred <- sum(w_valid^2) / eff_n_valid
      }
    }
    
    if (n_assigned > 0) {
      x_assigned <- x_vals[assigned_idx]
      sx_assigned <- sx_vals[assigned_idx]
      if (n_assigned > 1) {
        obs_variance <- var(x_assigned)
        noise_variance <- mean(sx_assigned^2)
        x_spread <- max(0, 1 - (noise_variance / obs_variance)) 
        if (n_assigned >= 4) {
          sq_dist_x <- (x_assigned)^2 
          sum_sq_x <- sum(sq_dist_x)
          if (sum_sq_x > 0) {
            safe_shares <- pmax(sq_dist_x / sum_sq_x, 1e-15) 
            lev_entropy <- (-sum(safe_shares * log(safe_shares))) / log(n_assigned) 
          } else { lev_entropy <- 0 }
        }
      } else { x_spread <- 0 }
      
      cluster_d2 <- raw_d2[assigned_idx]
      cluster_d2 <- cluster_d2[is.finite(cluster_d2)]
      if (length(cluster_d2) > 2) {
        m1 <- mean(cluster_d2); m2 <- var(cluster_d2)
        ad_p <- suppressWarnings(goftest::ad.test(cluster_d2, "pchisq", df = 1))$p.value
      }
    }
    
    if (eval(parse(text = diagnostic_logic_dict[["Insufficient N"]]))) {
      decision <- "Insufficient N"
    } else {
      failures <- c()
      for (rule_name in names(diagnostic_logic_dict)) {
        if (rule_name == "Insufficient N") next
        if (eval(parse(text = diagnostic_logic_dict[[rule_name]]))) failures <- c(failures, rule_name)
      }
      decision <- if (length(failures) == 0) "Currently Homogeneous" else paste(failures, collapse = " + ")
    }
    
    return(list(n_confident=n_confident, cred_ratio=cred_ratio, cont_cred=cont_cred, 
                x_spread=x_spread, lev_entropy=lev_entropy, ad_pval=ad_p, m1=m1, m2=m2, 
                wt_d2=wt_d2, wt_pval=wt_pval, df_eff=df_eff, mean_f=mean_f, diagnostic_decision=decision))
  }
  
  # =========================================================================
  # 4. COMPUTE METRICS & UPDATE MODEL OBJECT
  # =========================================================================
  diag_results <- lapply(1:n_fitted_comp, function(c) {
    get_cluster_diagnostics(w_i = post_data[, c], assigned_idx = which(classes == c), 
                            raw_d2 = raw_d2_list[[c]], x_vals = x_vals, sx_vals = sx_vals)
  })
  
  best_model_params <- cbind(best_model_params, dplyr::bind_rows(diag_results))
  
  model$profiler_metrics_raw <- best_model_params %>%
    dplyr::mutate(component = dplyr::row_number()) %>%
    dplyr::select(component, slope = b, n_confident, cred_ratio, cont_cred, x_spread, 
                  lev_entropy, ad_pval, m1, m2, wt_d2, wt_pval, mean_f, diagnostic_decision)
  
  model$profiler_logic_dictionary <- diagnostic_logic_dict
  
  slopes_df <- best_model_params %>%
    dplyr::rowwise() %>%
    dplyr::mutate(
      annotation_font_size = 0.4,
      ad_safe = ifelse(is.na(ad_pval), 1.0, ad_pval),
      lev_safe = ifelse(is.na(lev_entropy), 0.0, lev_entropy),
      wt_p_safe = ifelse(is.na(wt_pval), 1.0, wt_pval),
      label_text = sprintf("N_conf: %d | IGX2: %.2f | LevEnt: %.2f\nCred: %.2f | wtF-stat: %.1f |\nAD p: %.3f | m1: %.2f\nwtD2: %.2f (df=%.1f,p=%.3f)",
                           n_confident, x_spread, lev_safe, cred_ratio, mean_f, ad_safe, 
                           m1, wt_d2, df_eff, wt_p_safe)
    ) %>% dplyr::ungroup()
    
  model$profiler_metrics_formatted <- slopes_df

  # =========================================================================
  # 5. EXECUTE PLOTTING
  # =========================================================================
  base_colors <- c("#1b9e77", "#d95f02", "#7570b3", "#e7298a",
                   "#66a61e", "#e6ab02", "#a6761d", "#e41a1c",
                   "#377eb8", "#984ea3")
                   
  master_plot_name <- paste0(prefix, "_line_profiles.png")
  png(master_plot_name, width = 12, height = 8, units = "in", res = 200)

  par(oma = c(0, 0, 5, 0), mar = c(4.5, 5.5, 3, 16))
  
  wrapped_decisions <- sapply(slopes_df$diagnostic_decision, function(text) {
    paste(strwrap(text, width = 25), collapse = "\n  ")
  })
  custom_legend <- sprintf("Comp %d:\n  %s", seq_len(n_fitted_comp), wrapped_decisions)

  has_noise <- "ProbNoise" %in% colnames(post_data)

  if (has_noise) {
    magic_colors <- c(base_colors[seq_len(n_fitted_comp)], "#cccccc", "#666666")
    unassigned_idx <- n_fitted_comp + 2; noise_idx <- n_fitted_comp + 1
  } else {
    magic_colors <- c(base_colors[seq_len(n_fitted_comp)], "#666666") 
    unassigned_idx <- n_fitted_comp + 1
  }

  wrapper_classes <- apply(post_data, MARGIN = 1, which.max)
  wrapper_maxprob <- apply(post_data, MARGIN = 1, max)
  wrapper_classes[wrapper_classes <= n_fitted_comp & wrapper_maxprob < 0.95] <- unassigned_idx
  
  unassigned_label <- sprintf("Unassigned (N=%d)", sum(wrapper_classes == unassigned_idx))

  if (has_noise) {
    noise_label <- sprintf("Unstructured (N=%d)", sum(wrapper_classes == noise_idx))
    custom_legend <- c(custom_legend, noise_label, unassigned_label)
    leg_lty <- c(rep(1, n_fitted_comp), NA, NA)  
    leg_lwd <- c(rep(3, n_fitted_comp), NA, NA)  
    leg_pch <- c(rep(NA, n_fitted_comp), 15, 15) 
  } else {
    custom_legend <- c(custom_legend, unassigned_label)
    leg_lty <- c(rep(1, n_fitted_comp), NA)  
    leg_lwd <- c(rep(3, n_fitted_comp), NA)  
    leg_pch <- c(rep(NA, n_fitted_comp), 15) 
  }

  total_legend_items <- length(custom_legend)
  leg_y_gap <- ifelse(total_legend_items > 6, 1.4, 2.5)
  leg_txt_size <- ifelse(total_legend_items > 6, 0.65, 0.8)
  leg_font <- ifelse(total_legend_items > 6, 1, 2) 

  if (gene_label == TRUE){
      gene_label_name = seq_len(n_fitted_comp)
  } else {
      gene_label_name = NULL
  }
  plot(model, colors = magic_colors, mgp = c(3.5, 1, 0), mar = c(4.5, 5.5, 3, 16), label_comp = gene_label_name, legend = FALSE)
  
  par(xpd = NA) 
  usr <- par("usr"); xmin=usr[1]; xmax=usr[2]; ymin=usr[3]; ymax=usr[4]
  
  legend("topright", inset = c(-0.35, 0), legend = custom_legend, col = magic_colors, lty = leg_lty, lwd = leg_lwd, pch = leg_pch, pt.cex = leg_txt_size * 2.2, text.col = magic_colors, bty = "n", cex = leg_txt_size, text.font = leg_font, y.intersp = leg_y_gap)
          
  mx <- (xmax - xmin) * 0.02; my <- (ymax - ymin) * 0.02

  for (j in seq_len(nrow(slopes_df))) {
    cand_x <- seq(max(0, xmin), xmax, length.out=1000); cand_y <- slopes_df$b[j] * cand_x
    valid <- which(cand_y >= ymin & cand_y <= ymax)
    if(length(valid)>0){ idx <- valid[max(1, round(length(valid)*0.85))]; ax=cand_x[idx]; ay=cand_y[idx]
    } else { ax=(xmax+max(0,xmin))/2; ay=(ymax+ymin)/2 }
    
    adj_x <- 1; adj_y <- ifelse(ay > (ymin + ymax)/2, 1, 0)
    label_cex <- slopes_df$annotation_font_size[j]
    tw <- strwidth(slopes_df$label_text[j], cex=label_cex, font=2)
    th <- strheight(slopes_df$label_text[j], cex=label_cex, font=2)
    
    bx1 <- ax - tw; bx2 <- ax; by1 <- ifelse(adj_y == 1, ay - th, ay); by2 <- ifelse(adj_y == 1, ay, ay + th)
    if (bx2 > xmax - mx) ax <- ax - (bx2 - (xmax - mx))
    if (bx1 < xmin + mx) ax <- ax + ((xmin + mx) - bx1)
    if (by2 > ymax - my) ay <- ay - (by2 - (ymax - my))
    if (by1 < ymin + my) ay <- ay + ((ymin + my) - by1)
    
    text(ax, ay, labels=slopes_df$label_text[j], adj=c(adj_x, adj_y), cex=label_cex, font=2, col=magic_colors[j])
  }
  par(xpd = FALSE)
  
  sj_pval <- model$scoutjoy_test$Pvalue
  sj_reject <- ifelse(sj_pval < 0.05, "TRUE", "FALSE")
  title_line1 <- sprintf("%s | SJ Reject: %s", prefix, sj_reject)

  base_cex <- 1.2
  max_width_allowed <- 0.7 
  current_width <- strwidth(title_line1, units="figure", cex=base_cex)
  fit_cex <- if(current_width > max_width_allowed) base_cex * (max_width_allowed / current_width) else base_cex
  fit_cex <- max(fit_cex, 0.6)

  par(oma = c(0, 0, 5, 0)) 
  mtext(title_line1, outer = TRUE, side = 3, line = 1, font = 2, cex = fit_cex, adj = 0.5) 
  saved_plot <- recordPlot()
  dev.off()
  cat(sprintf("  -> Saved Scatter Plot: %s\n", master_plot_name))
  
  # Return the fully updated model!
  return(list(
      model = model,
      plot = saved_plot
  ))
}
