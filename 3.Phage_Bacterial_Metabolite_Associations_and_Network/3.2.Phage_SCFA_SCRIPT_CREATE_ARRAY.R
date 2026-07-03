## Author: 
    # Panayiotis Louca 

## Clear environment 
    rm(list = ls()) 

## Set seed 
    set.seed(1234)
      
## load up packages: 
    
    ### core 
    library(tidyverse)
    
    # mixed effect modelling 
    library(glmmTMB)
    library(lmerTest)
    library(broom.mixed)
    
    library(parallel)
    library(future.apply)
    
# -------------------------------------------------------------------------- # 
    
# ************************* # 
#   IMPORT & PREP DATA   ---- 
# ************************* # 
    
    ##   Dataset ---- 
    
    # HPC 
    df <- read_rds(file.path("/scratch/users/k2480753/resistome/data/Revised_phages_microbes_metabs_DATASET.rds")) %>% 
      as.data.frame() %>%
      tibble::remove_rownames() 
    
    
    # vectorise columns 
    phage_columns <- grep("^phage_[0-9]+", names(df), value = TRUE)
    scfa_columns = grep("serum_SCFA_Acetate|serum_SCFA_Butyrate|serum_SCFA_Propionate|stool_SCFA_Acetate|stool_SCFA_Butyrate|stool_SCFA_Propionate", names(df), value = TRUE)
    
    covar_columns = c('gm_age','gm_sex','gm_BMI','gm_batch')
    
# -------------------------------------------------------------------------- #  
    
    df_analyse <- df %>%
      select(
        iid,
        all_of(covar_columns),
        fid,
        # phages 
        all_of(phage_columns),
        # scfa  
        all_of(scfa_columns)
      )
    
# -------------------------------------------------------------------------- #  
    
    # Scale/normalise 
    #- qunatile normalisation function 
    rankit = function(vect) {
      # Nan position 
      vect.position = is.na(vect)
      # Rank non Nans 
      qnormvect = qnorm((rank(vect[!vect.position]) - 1/2) / length(vect[!vect.position]))
      # Replace non-Nans with rank transformed values 
      vect[!vect.position] = qnormvect
      return(vect)
    }
    
    # z score age 
    df_analyse$gm_age = (df_analyse$gm_age - mean(df_analyse$gm_age))/sd(df_analyse$gm_age)
    df_analyse$gm_BMI = (df_analyse$gm_BMI - mean(df_analyse$gm_BMI))/sd(df_analyse$gm_BMI)
    
    # CLR normalise phages and taxa 
    
    # Create a function to CLR normalise 
    CLRnorm <- function(features) {
      if (!is.data.frame(features) & !is.matrix(features)) {
        stop("Input must be a data frame or matrix.")
      }
      features_norm <- as.matrix(features)
      features_norm[features_norm == 0] <- 1e-6
      features_CLR <- chemometrics::clr(features_norm)
      as.data.frame(features_CLR, col.names = colnames(features))
    }
    
# -------------------------------------------------------------------------- #  
    
    ###   setup data ---- 
    
      # gut phage data only 
    dat_phage <- df_analyse %>%
      select(iid,
             all_of(phage_columns)
      )
    
    dat_phage <- # move iid to rownames 
      dat_phage %>%
      tibble::column_to_rownames('iid')
    
    # Apply CLR transformation to phage columns  
    dat_phage <- CLRnorm(dat_phage)
    
    # Apply inverse normalisation to phage columns 
    dat_phage <- dat_phage %>%
      mutate(across(everything(), rankit))
    
        # metabolite data only 
    dat_metabs <- df_analyse %>%
      select(iid,
             all_of(scfa_columns))
    
    dat_metabs <- dat_metabs %>%
      tibble::column_to_rownames('iid')

    dat_metabs <- lapply(dat_metabs, function(x) (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE))
    
# -------------------------------------------------------------------------- #  
    
    # Function for error and warning reporting 
    # REF: https://stackoverflow.com/a/4952908/10360530 
    catch <- function(expr) {
      warn <- err <- NULL
      value <- withCallingHandlers(
        tryCatch(expr, error=function(e) {
          err <<- e
          NULL
        }), warning=function(w) {
          warn <<- w
          invokeRestart("muffleWarning")
        })
      list(v=value, w=warn, e=err)
    }
    
# -------------------------------------------------------------------------- #  
    
    set.seed(1234)
    
    ######   CREATE LMER FUNCTION  ---- 
    
    # Create formulas for unadjusted and adjusted models 
    formula_unadjusted <- "y ~ x"
    formula_adjusted <- "y ~ x + gm_age + gm_sex + gm_BMI + (1 | fid) + (1 | gm_batch)"
    
    # REGRESSION FUNCTION 
    run_regression_parallel <- function(dependent_vars,
                                        independ_vars,
                                        covars = NULL,
                                        cores = parallel::detectCores() - 1) {
      start_time <- Sys.time()
      
      dependent_vars <- as.data.frame(dependent_vars)
      independ_vars <- as.data.frame(independ_vars)
      
      # Setup parallel backend 
      future::plan("multisession", workers = cores)
      progressr::handlers(global = TRUE)
      
      # All (i,j) combinations 
      combos <- expand.grid(
        depend_idx = seq_along(dependent_vars),
        independ_idx = seq_along(independ_vars)
      )
      
      results_list <- progressr::with_progress({
        p <- progressr::progressor(along = 1:nrow(combos))
        future_lapply(1:nrow(combos), function(k) {
          i <- combos$depend_idx[k]
          j <- combos$independ_idx[k]
          p()  # update progress 
          
          depend_var_name <- colnames(dependent_vars)[i]
          independ_var_name <- colnames(independ_vars)[j]
          y <- dependent_vars[[i]]
          x <- independ_vars[[j]]
          
          df_model <- covars
          df_model$y <- y
          df_model$x <- x
          
          n_total <- nrow(df_model)
          n_not_0_depend_var <- sum(y != 0, na.rm = TRUE)
          n_not_0_independ_var <- sum(x != 0, na.rm = TRUE)
          
          # Skip if too few observations 
          if (n_total < 10) {
            return(data.frame(
              depend_var = depend_var_name,
              independ_var = independ_var_name,
              estimate_unadjusted = NA,
              std_err_unadjusted = NA,
              ci_low_unadjusted = NA,
              ci_upper_unadjusted = NA,
              pval_unadjusted = NA,
              estimate_adjusted = NA,
              std_err_adjusted = NA,
              ci_low_adjusted = NA,
              ci_upper_adjusted = NA,
              pval_adjusted = NA,
              formula_unadjusted = formula_unadjusted,
              formula_adjusted = formula_adjusted,
              n_total = n_total,
              n_not_0_depend_var = n_not_0_depend_var,
              n_not_0_independ_var = n_not_0_independ_var,
              error_message_unadjusted = "Too few observations",
              error_message_adjusted = "Too few observations",
              warning_message_unadjusted = NA,
              warning_message_adjusted = NA,
              stringsAsFactors = FALSE
            ))
          }
          
          # UNADJUSTED model 
          unadj_result <- tryCatch({
            model_unadj <- glmmTMB(as.formula(formula_unadjusted),
                                   data = df_model, REML = FALSE)
            model_output_unadj <- broom.mixed::tidy(model_unadj, conf.int = TRUE) %>%
              dplyr::filter(term == "x") %>%
              dplyr::select(estimate, std.error, conf.low, conf.high, p.value) %>%
              dplyr::rename(
                estimate_unadjusted = estimate,
                std_err_unadjusted = std.error,
                ci_low_unadjusted = conf.low,
                ci_upper_unadjusted = conf.high,
                pval_unadjusted = p.value
              )
            if (nrow(model_output_unadj) == 0)
              stop("No unadjusted model output for variable")
            list(output = model_output_unadj, error = NA, warning = NA)
          }, warning = function(w) {
            list(output = data.frame(
              estimate_unadjusted = NA,
              std_err_unadjusted = NA,
              ci_low_unadjusted = NA,
              ci_upper_unadjusted = NA,
              pval_unadjusted = NA
            ), error = NA, warning = w$message)
          }, error = function(e) {
            list(output = data.frame(
              estimate_unadjusted = NA,
              std_err_unadjusted = NA,
              ci_low_unadjusted = NA,
              ci_upper_unadjusted = NA,
              pval_unadjusted = NA
            ), error = e$message, warning = NA)
          })
          
          # ADJUSTED model 
          adj_result <- tryCatch({
            model_adj <- glmmTMB(as.formula(formula_adjusted),
                                 data = df_model, REML = FALSE)
            model_output_adj <- broom.mixed::tidy(model_adj, conf.int = TRUE) %>%
              dplyr::filter(term == "x") %>%
              dplyr::select(estimate, std.error, conf.low, conf.high, p.value) %>%
              dplyr::rename(
                estimate_adjusted = estimate,
                std_err_adjusted = std.error,
                ci_low_adjusted = conf.low,
                ci_upper_adjusted = conf.high,
                pval_adjusted = p.value
              )
            if (nrow(model_output_adj) == 0)
              stop("No adjusted model output for variable")
            list(output = model_output_adj, error = NA, warning = NA)
          }, warning = function(w) {
            list(output = data.frame(
              estimate_adjusted = NA,
              std_err_adjusted = NA,
              ci_low_adjusted = NA,
              ci_upper_adjusted = NA,
              pval_adjusted = NA
            ), error = NA, warning = w$message)
          }, error = function(e) {
            list(output = data.frame(
              estimate_adjusted = NA,
              std_err_adjusted = NA,
              ci_low_adjusted = NA,
              ci_upper_adjusted = NA,
              pval_adjusted = NA
            ), error = e$message, warning = NA)
          })
          
          # Combine results
          combined_results <- cbind(
            depend_var = depend_var_name,
            independ_var = independ_var_name,
            unadj_result$output,
            adj_result$output,
            formula_unadjusted = formula_unadjusted,
            formula_adjusted = formula_adjusted,
            n_total = n_total,
            n_not_0_depend_var = n_not_0_depend_var,
            n_not_0_independ_var = n_not_0_independ_var,
            error_message_unadjusted = unadj_result$error,
            error_message_adjusted = adj_result$error,
            warning_message_unadjusted = unadj_result$warning,
            warning_message_adjusted = adj_result$warning
          )
          return(combined_results)
        })
      })
      
      # Combine all into single dataframe 
      return(bind_rows(results_list))
    }
    
# -------------------------------------------------------------------------- #  
    
    ######   SPLIT INTO TASKS  ---- 
    
    options(future.globals.maxSize = 2000 * 1024^2)  # Increase to 2 GiB 
    
    cores <- as.integer(Sys.getenv("SLURM_NTASKS")) # get n cores 
    
    args <- commandArgs(trailingOnly = TRUE)
    index <- as.integer(args[1])
    
    # define ranges of columns to test for each job 
    
    # split n phages to test into 7 batches 
    n_jobs <- 7 # Number of array jobs 
    
    n_cols <- length(phage_columns) # Number of phage columns 
    
    cols_per_job_base <- floor(n_cols / n_jobs) 
    remainder <- n_cols %% n_jobs  
    
    # Create a list to store the ranges 
    ranges <- list()
    
    # Distribute columns across 7 jobs 
    current_start <- 1
    for (i in 1:7) {
      
      cols_this_job <- cols_per_job_base + ifelse(i <= remainder, 1, 0)
      
      start_idx <- current_start
      end_idx <- current_start + cols_this_job - 1
      ranges[[as.character(i)]] <- start_idx:end_idx
      
      current_start <- end_idx + 1
    }
    
    # Verify the split 
    sapply(ranges, length)  
    sum(sapply(ranges, length))  
   
# -------------------------------------------------------------------------- #  
          
          selected_columns <- ranges[[index]]
          if (is.null(selected_columns)) stop("Invalid index: no column range defined.")
          
          # Run the regression 
          regression_output <- run_regression_parallel(
            independ_vars = dat_phage[ ,selected_columns],
            dependent_vars = dat_metabs,
            covars = df_analyse %>% select(gm_age, gm_sex, gm_BMI, fid, gm_batch),
            cores = cores
          )
          
# -------------------------------------------------------------------------- #  
          
          # Save results for this job 
          path <- file.path("/scratch/users/k2480753/resistome/Phages_SCFA_and_Bile_acids/Phages_SCFA/results", paste0("phages_SCFA_RESULTS_", index, ".csv"))
          write.csv(regression_output,
                    path,
                    row.names = FALSE)
          
          
############################################################################################### 
##################################### END ##################################################### 
############################################################################################### 
