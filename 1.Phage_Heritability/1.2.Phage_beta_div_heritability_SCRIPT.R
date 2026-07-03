## Author: 
# Panayiotis Louca 

## Clear environment 
    rm(list = ls()) 

## Set seed 
    set.seed(1234)
      
## load up packages: 

    ### core 
    library(tidyverse)
    
    # for diversity 
    library(vegan)
    
    # for heritability 
    library(mets)
    library(lava) # to compare models 

# -------------------------------------------------------------------------- # 
  
  
# ******************** # 
#   Import dataset  ---- 
# ******************** # 
   
    df <- read_rds(file.path("/Users/panayiotislouca/Documents/KCL_files/Phages_metagenome_characterisation/Data/Revised_all_phages_microbes_metabs_DATASET.rds")) %>%
      as.data.frame() %>%
      tibble::remove_rownames() 
    
    # vectorise column names 
    phage_columns <- grep("^phage_[0-9]+", names(df), value = TRUE)
    
# -------------------------------------------------------------------------- # 
    
# --------------------------------------------- # 
######   Distance-based reduncancy appraoch  ---- 
# --------------------------------------------- # 
    
    # Split phages 
    phageome <- df %>% select(all_of(phage_columns))
    
    # Bray–Curtis distances from phage abundance matrix 
    phage_dist <- vegdist(phageome, method = "bray")
  
    # Adonis with zygosity and covariates 
    adonis_res <- adonis2( 
      phage_dist ~ gm_age + gm_BMI + gm_sex + gm_batch + gm_zygosity,
      data = df,
      permutations = 999,
      strata = df$fid # keep MZ/DZ pairs together 
    )
    
    adonis_res_terms <- adonis2(
      phage_dist ~ gm_zygosity + gm_age + gm_BMI + gm_sex + gm_batch,
      data = df,
      permutations = 999,
      strata = df$fid
    )
    
# -------------------------------------------------------------------------- # 
    
# ------------------------------------------------------------- # 
######   Standard SEM approach using PCOA & SEM on first PC  ---- 
# ------------------------------------------------------------- # 

    # keep top 5 axes 
    phage_ord <- cmdscale(phage_dist, k = 5, eig = TRUE)
    phage_axes <- as.data.frame(phage_ord$points)
    colnames(phage_axes) <- paste0("PCo", 1:5)
    
    # Eigenvalues 
    eig_vals <- phage_ord$eig
    
    # Proportion of variance explained per axis 
    var_explained <- eig_vals / sum(eig_vals[eig_vals > 0])
    
    # Cumulative variance explained 
    cum_var <- cumsum(var_explained)
    
    axis_info <- data.frame(
      Axis = paste0("PCo", 1:length(var_explained)),
      Eigenvalue = eig_vals,
      Proportion = round(var_explained, 4)*100,
      Cumulative = round(cum_var, 4)*100
    )
    
    head(axis_info, 10)  # show first 10 axes 
    
    
    # Build phenotype dataframe 
    df_analyse <- df %>%
      select(iid, fid, gm_zygosity, gm_age, gm_BMI, gm_sex)
    
    # Merge pheno and axes 
    df_axes <- cbind(df_analyse, phage_axes)
    
    # combine PCo1 & PCo2 into new variable
    df_axes$PCo1_2 <- df_axes$PCo1 + df_axes$PCo2
    
# -------------------------------------------------------------------------- #  
    
          # ******************************** # 
          ######   Estimate heritability  ---- 
          # ******************************** # 
    
      phen <- data.frame(
        fid = df_axes$fid,
        iid = df_axes$iid,
        zygosity = df_axes$gm_zygosity,
        age = df_axes$gm_age,
        bmi = df_axes$gm_BMI,
        sex = df_axes$gm_sex
      )
      
      dat <- t(data.frame(PCo1 = df_axes[['PCo1']],
                          PCo2 = df_axes[['PCo2']],
                          PCo1_2 = df_axes[['PCo1_2']]))
      colnames(dat) <- df_axes$fid
      dat <- as.data.frame(dat)
      
      id <- phen$fid
      age <- phen$age
      zyg <- factor(phen$zygosity)
      bmi <- phen$bmi
      sex <- phen$sex
    
################################################ data input ###################################################
    
    res<-{}
    cor_sat <- {}
    model_tab <-{}
    
    for (i in 1:nrow(dat)) {
      
      methy <- as.numeric(dat[i, ])
      
      data <-	data.frame(methy, zyg, id, age, bmi)
      
      
####################################### Assign formula and run ############################################# 
      
    formula <- as.formula("methy ~ age + bmi + sex")

    ace <- summary(twinlm(formula, id = "id", DZ = "DZ", zyg = "zyg", data = data, missing = FALSE, type = "ace"))
    ae <- summary(twinlm(formula, id = "id", DZ = "DZ", zyg = "zyg", data = data, missing = FALSE, type = "ae"))
    ce <- summary(twinlm(formula, id = "id", DZ = "DZ", zyg = "zyg", data = data, missing = FALSE, type = "ce"))
      
    e <- summary(twinlm(formula, id = "id", DZ = "DZ", zyg = "zyg", data = data, missing = FALSE, type = "e"))
      
    sat <- summary(twinlm (formula, id = "id", DZ = "DZ", zyg = "zyg", data = data, missing = FALSE, type = "sat"))
      
############################### P value for models ##############################################
      
    ace2<- twinlm(formula, id = "id", DZ = "DZ", zyg = "zyg", data = data, missing = FALSE, type = "ace")
    ae2<- twinlm(formula, id = "id", DZ = "DZ", zyg = "zyg", data = data, missing = FALSE, type = "ae")
    ce2<- twinlm(formula, id = "id", DZ = "DZ", zyg = "zyg", data = data, missing = FALSE, type = "ce")
      
    e2<- twinlm(formula, id = "id", DZ = "DZ", zyg = "zyg", data = data, missing = FALSE, type = "e")
      
    sat2<- twinlm (formula, id = "id", DZ = "DZ", zyg = "zyg", data = data, missing = FALSE, type = "sat")
      
##################################### Compare models ##################################################### 
    
    comparison_model<- sat2

    pval_ACE<- compare(ace2,sat2)
    Rpval_ACE<- pval_ACE$p.value
      
      
    pval_AE<- compare(ae2,comparison_model)
    Rpval_AE<- pval_AE$p.value
      
    pval_CE<- compare(ce2,comparison_model)
    Rpval_CE<- pval_CE$p.value
      
    pval_E<- compare(e2,comparison_model)
    Rpval_E<- pval_E$p.value
      
#################################   Generate comparison metrics    ###########################################

      ### Saturated model ###
      sat_2LL<- -2*logLik(sat2)
      sat_AIC<- sat$AIC 
      sat_df<- attr(logLik(sat2), "df")
      
      ### ACE ###
      ACE_2LL<- -2*logLik(ace2)
      ACE_AIC<- ace$AIC 
      ACE_df<- attr(logLik(ace2), "df")
      ACE_DEL_2LL<- ACE_2LL - sat_2LL
      ACE_DEL_df<- pval_ACE$parameter[1]
      
      ACE_CHI<- pval_ACE$statistic[1]
      
      ### AE ###
      
      AE_2LL<- -2*logLik(ae2)
      AE_AIC<- ae$AIC 
      AE_df<- attr(logLik(ae2), "df")
      AE_DEL_2LL<- AE_2LL - sat_2LL
      AE_DEL_df<- pval_AE$parameter[1]
      AE_CHI<- pval_AE$statistic[1]
      
      ### CE ### 
      CE_2LL <- -2*logLik(ce2)
      CE_AIC <- ce$AIC
      CE_df <- attr(logLik(ce2),"df")
      CE_DEL_2LL <- CE_2LL - sat_2LL
      CE_DEL_df <- pval_CE$parameter[1]
      CE_CHI <- pval_CE$statistic[1]
      
      ### E ###
      E_2LL<- -2*logLik(e2)
      E_AIC<- e$AIC 
      E_df<- attr(logLik(e2), "df")
      E_DEL_2LL<- E_2LL - sat_2LL
      E_DEL_df<- pval_E$parameter[1]
      E_CHI<- pval_E$statistic[1]
      
################################# Generate and format result table ###########################################
      
      res1<- rbind(
        c(row.names(dat)[i],"ACE",ace$acde[1,],ace$acde[2,],"-","-","-",ace$acde[3,],ace$corMZ,ace$corDZ,ace$zyg[1],ace$zyg[2],"ACE",ACE_2LL,ACE_df,ACE_AIC,ACE_DEL_2LL,ACE_DEL_df,ACE_CHI,Rpval_ACE),
        c(row.names(dat)[i],"AE", ae$acde[1,],"-","-","-","-","-","-",ae$acde[2,],ae$corMZ,ae$corDZ,ae$zyg[1],ae$zyg[2],"AE",AE_2LL,AE_df,AE_AIC,AE_DEL_2LL,AE_DEL_df,AE_CHI,Rpval_AE),
        c(row.names(dat)[i],"E", "-","-","-","-","-","-","-","-","-",e$acde[1,],e$corMZ,e$corDZ,e$zyg[1],e$zyg[2],"E",E_2LL,E_df,E_AIC,E_DEL_2LL,E_DEL_df,E_CHI,Rpval_E),
        c(row.names(dat)[i],"CE", "-","-","-",ce$acde[1,],"-","-","-",ce$acde[2,],ce$corMZ,ce$corDZ,ce$zyg[1],ce$zyg[2],"CE",CE_2LL,CE_df,CE_AIC,CE_DEL_2LL,CE_DEL_df,CE_CHI,Rpval_CE)
      )
      
      res<- rbind(res, res1)
      
      cor_sat1<- rbind(c(sat$corMZ,sat$corDZ))
      
      cor_sat<- rbind(cor_sat,c(row.names(dat)[i],cor_sat1))
      
      mod_tab1 <- rbind(
        c(row.names(dat)[i],"saturated",sat_2LL,sat_df,sat_AIC,"-","-","-","-"),
        c(row.names(dat)[i],"ACE",ACE_2LL,ACE_df,ACE_AIC,ACE_DEL_2LL,ACE_DEL_df,ACE_CHI,Rpval_ACE),
        c(row.names(dat)[i],"AE",AE_2LL,AE_df,AE_AIC,AE_DEL_2LL,AE_DEL_df,AE_CHI,Rpval_AE),
        c(row.names(dat)[i],"E",E_2LL,E_df,E_AIC,E_DEL_2LL,E_DEL_df,E_CHI,Rpval_E)
      )
      
      model_tab<-rbind(model_tab,mod_tab1)
      
    }
    
########################################## Output result table ###############################################
    
    head <- c("Trait",
              "Model","A","A_2.50%","A_97.50%","C","C_2.50","C_97.50%","D","D_2.50","D_97.50%","E","E_2.50%","E_97.50%","Corr_MZ","Corr_MZ_2.50%","Corr_MZ_97.50%","Corr_DZ","Corr_DZ_2.50%","Corr_DZ_97.50%","MZ:pairs/singletons","DZ:pairs/singletons","Model","-2LL","df","AIC","D-2LL","D-df","chi","Pval")
    resF <- rbind(head,res)
    
    
    head_sat <- c("Trait",
                  "Corr_MZ","Corr_MZ_2.50%","Corr_MZ_97.50%",
                  "Corr_DZ","Corr_DZ_2.50%","Corr_DZ_97.50%")
    res_sat <- rbind(head_sat,cor_sat)
    
    
    head_mod <- c("Trait",
                  "Model","-2LL","df","AIC","D-2LL","D-df","chi","Pval")
    res_mod <- rbind(head_mod,model_tab)
    
########################################## Save result table ###############################################
    
    write.table(resF,
                file = "/Users/panayiotislouca/Documents/KCL_files/Phages_metagenome_characterisation/Analysis/Phages_heritability/Beta_diversity/heritability_results/Output_heritibility_beta_diversity.txt",
                col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")
    
    write.table(res_sat, 
                file = "/Users/panayiotislouca/Documents/KCL_files/Phages_metagenome_characterisation/Analysis/Phages_heritability/Beta_diversity/heritability_results/beta_diversity_CORR_SAT.txt",
                col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")
    
    write.table(res_mod, 
                file = "/Users/panayiotislouca/Documents/KCL_files/Phages_metagenome_characterisation/Analysis/Phages_heritability/Beta_diversity/heritability_results/beta_diversity_MOD.txt",
                col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")
    
############################################################################################### 
##################################### END ##################################################### 
############################################################################################### 