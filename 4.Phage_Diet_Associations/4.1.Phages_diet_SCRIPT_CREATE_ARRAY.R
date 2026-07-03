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
    
    # HPC 
    df <- read_rds(file.path("/scratch/users/k2480753/resistome/data/Revised_phages_microbes_metabs_DATASET.rds")) %>%
      as.data.frame()
    
    # vectorise column names 
    phage_columns <- grep("^phage_[0-9]+", names(df), value = TRUE)
    diet_columns = grep("diet_", names(df), value = TRUE)
    
    covar_columns = c('gm_age','gm_sex','gm_BMI','gm_batch')
    
# -------------------------------------------------------------------------- #  
    
    # subset dataset to necessary variables 
    df_analyse <- df %>%
      select(
        iid,
        all_of(covar_columns),
        fid,
        # diet 
        all_of(diet_columns),
        # phages 
        all_of(phage_columns)
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
    
    # z score age & BMI 
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
    
    # dietary data only 
    dat_diet <- df_analyse %>%
      select(iid,
             all_of(diet_columns))
    
    dat_diet <- # move iid to rownames 
      dat_diet %>%
      tibble::column_to_rownames('iid')
    
    dat_diet <- dat_diet %>%
      mutate(across(everything(), rankit))
    
    
    # gut phage data only 
    dat_phage <- df_analyse %>%
      select(iid,
             all_of(phage_columns))
    
    dat_phage <- # move iid to rownames 
      dat_phage %>%
      tibble::column_to_rownames('iid')
    
    dat_phage <- CLRnorm(dat_phage)
    
    dat_phage <- dat_phage %>%
      mutate(across(everything(), rankit))
    
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
    
    # create formula 
    formula_template <- "y ~ x + gm_age + gm_BMI + gm_sex + (1 | fid) + (1 | gm_batch)"
    
    # Mixed-effect regression analysis function with parallelisation and error/warning reporting 
    run_regression_parallel <- function(dependent_vars, independ_vars, covars = NULL, cores = parallel::detectCores() - 1) {
      start_time <- Sys.time()
      
      # Ensure inputs are data frames 
      dependent_vars <- as.data.frame(dependent_vars)
      independ_vars <- as.data.frame(independ_vars)
      
      
      # Setup parallel backend 
      future::plan("multisession", workers = cores)
      progressr::handlers(global = TRUE)
      
      
      # Create all (i,j) combinations 
      combos <- expand.grid(
        depend_idx = seq_along(dependent_vars),
        independ_idx = seq_along(independ_vars)
      )
      
      progressr::with_progress({
        p <- progressr::progressor(along = 1:nrow(combos))
        
        results_list <- future_lapply(1:nrow(combos), function(k) {
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
          
          n_total = nrow(df_model)
          n_not_0_depend_var <- sum(y != 0, na.rm = TRUE)
          n_not_0_independ_var <- sum(x != 0, na.rm = TRUE)
          
    
          # Skip if too few observations 
          if (n_total < 10) {
            return(data.frame(
              depend_var = depend_var_name,
              independ_var = independ_var_name,
              estimate = NA, std_err = NA, ci_low = NA, ci_upper = NA, pval = NA,
              formula = formula_template,
              n_total = n_total,
              n_not_0_depend_var = n_not_0_depend_var,
              n_not_0_independ_var = n_not_0_independ_var,
              error_message = "Too few observations",
              warning_message = NA,
              stringsAsFactors = FALSE
            ))
          }
          
          # Build model formula 
          formula <- as.formula(formula_template)
          
          
          result <- tryCatch({
            model <- glmmTMB(formula, data = df_model, REML = FALSE)
            
            model_output <- broom.mixed::tidy(model, conf.int = TRUE) %>%
              dplyr::filter(term == "x") %>%
              dplyr::select(term, estimate, std.error, conf.low, conf.high, p.value) %>%
              dplyr::rename(independ_var = term, std_err = std.error, ci_low = conf.low,
                            ci_upper = conf.high, pval = p.value) %>%
              dplyr::mutate(independ_var = independ_var_name)
            
            if (nrow(model_output) == 0) stop("No model output for variable")
            
            cbind(
              depend_var = depend_var_name,
              model_output,
              formula = formula_template,
              n_total = n_total,
              n_not_0_depend_var = n_not_0_depend_var,
              n_not_0_independ_var = n_not_0_independ_var,
              error_message = NA,
              warning_message = NA
            )
          }, warning = function(w) {
            data.frame(
              depend_var = depend_var_name,
              independ_var = independ_var_name,
              estimate = NA, std_err = NA, ci_low = NA, ci_upper = NA, pval = NA,
              formula = formula_template,
              n_total = n_total,
              n_not_0_depend_var = n_not_0_depend_var,
              n_not_0_independ_var = n_not_0_independ_var,
              error_message = NA,
              warning_message = w$message,
              stringsAsFactors = FALSE
            )
          }, error = function(e) {
            data.frame(
              depend_var = depend_var_name,
              independ_var = independ_var_name,
              estimate = NA, std_err = NA, ci_low = NA, ci_upper = NA, pval = NA,
              formula = formula_template,
              n_total = n_total,
              n_not_0_depend_var = n_not_0_depend_var,
              n_not_0_independ_var = n_not_0_independ_var,
              error_message = e$message,
              warning_message = NA,
              stringsAsFactors = FALSE
            )
          })
          
          return(result)
        })
      })
      
      final_results <- data.table::rbindlist(results_list, fill = TRUE)
      
      end_time <- Sys.time()
      message("Code took: ", round(difftime(end_time, start_time, units = "mins"), 2), " minutes")
      
      return(final_results)
    }

# -------------------------------------------------------------------------- #  
    
    options(future.globals.maxSize = 2000 * 1024^2)  # Increase to 2 GiB 
    
    cores <- as.integer(Sys.getenv("SLURM_NTASKS")) # get n cores 
    
    args <- commandArgs(trailingOnly = TRUE)
    index <- as.integer(args[1])
    
    # define ranges of columns to test for each job (nutrients are in columns 1:45) 
    ranges <- list(
      "1" = 1:5,
      "2" = 6:10,
      "3" = 11:15,
      "4" = 16:20,
      "5" = 21:25,
      "6" = 26:30,
      "7" = 31:35,
      "8" = 36:40,
      "9" = 41:45
    )
    
    selected_columns <- ranges[[index]]
    if (is.null(selected_columns)) stop("Invalid index: no column range defined.")
    
    # Run the regression 
    regression_output <- run_regression_parallel(
      independ_vars = dat_phage,
      dependent_vars = dat_diet[ , selected_columns],
      covars = df_analyse %>% select(gm_age, gm_sex, fid, gm_batch),
      cores = cores
    )
    
# -------------------------------------------------------------------------- #  
    
    # save file 
    path <- file.path("/scratch/users/k2480753/resistome/Phages_diet/results", paste0("phages_diet_RESULTS_", index, ".csv"))
    write.csv(regression_output,
              path,
              row.names = FALSE)
    
############################################################################################### 
##################################### END ##################################################### 
############################################################################################### 
    