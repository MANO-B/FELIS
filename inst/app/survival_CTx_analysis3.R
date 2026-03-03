# =========================================================================
# Tamura & Ikegami Model Ver 2.3.2 (Real-Data Port)
#  - censor == 1 means EVENT (death)
#  - Keep output variable names: time_pre, time_all, censor, Group
#  - Replaced downstream analysis/plotting with Proposed Simulation Pipeline:
#    Age-stratified external OS curve matching via weights on T_obs (time_all),
#    followed by weighted AFT modeling on observed data.
# =========================================================================

library(dplyr)
library(survival)
library(flexsurv)
library(ggplot2)
library(rlang)
library(patchwork)

# =========================================================================
# 0) Helpers
# =========================================================================

# Robust age parsing (e.g., "60代" -> 60)
extract_age_num <- function(x) {
  suppressWarnings(as.numeric(gsub("[^0-9.]", "", as.character(x))))
}

make_age_class_6grp <- function(age_num) {
  dplyr::case_when(
    is.na(age_num) ~ "全年齢",
    age_num < 40 ~ "40未満",
    age_num < 50 ~ "40代",
    age_num < 60 ~ "50代",
    age_num < 70 ~ "60代",
    age_num < 80 ~ "70代",
    age_num >= 80 ~ "80以上",
    TRUE ~ "全年齢"
  )
}

# Build registry/manual reference survival list (1-5y) with safe fallback
build_ref_surv_list_realdata <- function(input) {
  ref_surv_list <- list()
  age_groups_all <- c("40未満", "40代", "50代", "60代", "70代", "80以上", "全年齢")
  if (!is.null(input$survival_data_source)) {
    if (input$survival_data_source == "registry") {
      cancer_data <- Data_age_survival_5_year[[input$registry_cancer_type]]
      fallback_surv <- cancer_data[["全年齢"]]
      for (ag in age_groups_all) {
        if (!is.null(cancer_data[[ag]]) && length(cancer_data[[ag]]) == 5) {
          ref_surv_list[[ag]] <- as.numeric(cancer_data[[ag]])
        } else if (!is.null(fallback_surv) && length(fallback_surv) == 5) {
          ref_surv_list[[ag]] <- as.numeric(fallback_surv)
        }
      }
    } else if (input$survival_data_source == "manual_all") {
      vals <- c(input$surv_all_1y, input$surv_all_2y, input$surv_all_3y, input$surv_all_4y, input$surv_all_5y)
      for (ag in age_groups_all) ref_surv_list[[ag]] <- vals
    } else if (input$survival_data_source == "manual_age") {
      ref_surv_list[["40未満"]] <- c(input$surv_u40_1y, input$surv_u40_2y, input$surv_u40_3y, input$surv_u40_4y, input$surv_u40_5y)
      ref_surv_list[["40代"]]   <- c(input$surv_40s_1y, input$surv_40s_2y, input$surv_40s_3y, input$surv_40s_4y, input$surv_40s_5y)
      ref_surv_list[["50代"]]   <- c(input$surv_50s_1y, input$surv_50s_2y, input$surv_50s_3y, input$surv_50s_4y, input$surv_50s_5y)
      ref_surv_list[["60代"]]   <- c(input$surv_60s_1y, input$surv_60s_2y, input$surv_60s_3y, input$surv_60s_4y, input$surv_60s_5y)
      ref_surv_list[["70代"]]   <- c(input$surv_70s_1y, input$surv_70s_2y, input$surv_70s_3y, input$surv_70s_4y, input$surv_70s_5y)
      ref_surv_list[["80以上"]] <- c(input$surv_80s_1y, input$surv_80s_2y, input$surv_80s_3y, input$surv_80s_4y, input$surv_80s_5y)
      ref_surv_list[["全年齢"]] <- ref_surv_list[["60代"]] # fallback
    }
  }

  # hard fallback if empty or partially empty
  if (length(ref_surv_list) == 0) {
    for (ag in age_groups_all) ref_surv_list[[ag]] <- c(50, 30, 20, 15, 10)
  } else {
    for (ag in age_groups_all) {
      if (is.null(ref_surv_list[[ag]]) || length(ref_surv_list[[ag]]) != 5) {
        ref_surv_list[[ag]] <- ref_surv_list[[1]]
      }
    }
  }

  ref_surv_list
}

# =========================================================================
# Simulation Pipeline Port: External survival llogis smoothing
# =========================================================================
fit_llogis_from_S15 <- function(S_y_percent) {
  t_years <- 1:5
  S <- pmax(pmin(as.numeric(S_y_percent) / 100, 0.999), 0.001)

  obj <- function(par) {
    log_shape <- par[1]
    log_scale <- par[2]
    shape <- exp(log_shape)
    scale <- exp(log_scale)
    S_hat <- 1 / (1 + (t_years / scale)^shape)
    sum((S_hat - S)^2)
  }

  init <- c(log(1.5), log(3))
  opt <- tryCatch(
    optim(init, obj, method = "Nelder-Mead", control = list(maxit = 2000)),
    error = function(e) NULL
  )
  if (is.null(opt) || any(!is.finite(opt$par))) {
    return(list(shape = 1.5, scale = 3.0))
  }
  list(shape = exp(opt$par[1]), scale = exp(opt$par[2]))
}

llogis_surv <- function(t_years, shape, scale) {
  t_years <- pmax(t_years, 1e-8)
  1 / (1 + (t_years / scale)^shape)
}

# =========================================================================
# Simulation Pipeline Port: Direct T_obs Curve Matching (Age-stratified)
# IMPORTANT: censor==1 means death event, time_all is T_obs
# =========================================================================
calculate_obs_matching_weights_real <- function(data, ref_surv_list) {
  data$iptw <- 1.0

  # 0.5-year grid up to 5 years (10 points)
  t_points_years <- seq(0.5, 5.0, by = 0.5)
  t_points_days  <- t_points_years * 365.25

  if (!("age_class" %in% names(data))) stop("age_class missing")
  if (!all(c("time_all", "censor") %in% names(data))) stop("time_all/censor missing")

  for (ag in unique(data$age_class)) {
    if (!(ag %in% names(ref_surv_list))) next

    ag_idx <- which(data$age_class == ag)
    if (length(ag_idx) < 10) next
    ag_data <- data[ag_idx, , drop = FALSE]

    # ---- build target curve from external 1–5y points (smoothed by llogis) ----
    pars <- fit_llogis_from_S15(ref_surv_list[[ag]])
    S_target <- llogis_surv(t_points_years, pars$shape, pars$scale)
    S_target[!is.finite(S_target)] <- 0
    S_target <- pmax(pmin(S_target, 0.999), 0.001)

    # ---- smooth weight function of time_all (polynomial of log time) ----
    log_t <- log(pmax(ag_data$time_all / 365.25, 1e-6))
    log_t2 <- log_t^2
    log_t3 <- log_t^3

    obj_func <- function(theta) {
      w <- exp(theta[1] * log_t + theta[2] * log_t2 + theta[3] * log_t3)
      w <- pmax(w, 1e-6)
      w <- w / mean(w)

      km <- tryCatch(
        survival::survfit(survival::Surv(time_all, censor) ~ 1, data = ag_data, weights = w),
        error = function(e) NULL
      )
      if (is.null(km)) return(1e9)

      S_est <- approx(km$time, km$surv, xout = t_points_days, method = "constant", f = 0, rule = 2)$y
      S_est[!is.finite(S_est)] <- 0
      S_est <- pmax(pmin(S_est, 0.999), 0.001)

      loss <- sum((S_est - S_target)^2)
      loss + 0.001 * sum(theta^2)
    }

    opt <- tryCatch(
      optim(c(0, 0, 0), obj_func, method = "Nelder-Mead", control = list(maxit = 2000)),
      error = function(e) list(par = c(0, 0, 0))
    )
    th <- opt$par

    w_opt <- exp(th[1] * log_t + th[2] * log_t2 + th[3] * log_t3)
    w_opt <- pmax(w_opt, 1e-6)

    # winsorize
    lb <- stats::quantile(w_opt, 0.01, na.rm = TRUE)
    ub <- stats::quantile(w_opt, 0.99, na.rm = TRUE)
    w_opt <- pmax(lb, pmin(w_opt, ub))

    data$iptw[ag_idx] <- w_opt / mean(w_opt)
  }

  # normalize weights
  data$iptw <- as.numeric(data$iptw)
  data$iptw <- ifelse(is.na(data$iptw) | data$iptw <= 0, 1.0, data$iptw)
  data$iptw <- data$iptw / mean(data$iptw)

  ess <- (sum(data$iptw, na.rm = TRUE)^2) / sum(data$iptw^2, na.rm = TRUE)
  attr(data, "ESS") <- ess
  data
}

# =========================================================================
# 1) Main Logic Control Function (Preprocessing)  --- UNCHANGED
# =========================================================================
survival_CTx_analysis2_logic_control <- function() {
  analysis_env <- new.env()
  clear_reactive_data(OUTPUT_DATA)
  gc()

  with(analysis_env, {
    withProgress(message = sample(nietzsche)[1], {
      Data_case_target = Data_case()

      if (!is.null(input$gene_group_analysis) &&
          input$gene_group_analysis %in% c("Only cases with mutations in the gene set are analyzed",
                                           "Only cases without mutations in the gene set are analyzed")) {
        Data_case_target = Data_case_target %>%
          dplyr::filter(C.CAT調査結果.基本項目.ハッシュID %in% Data_report()$Tumor_Sample_Barcode)
      }

      Data_case_target = Data_case_target %>%
        dplyr::select(
          C.CAT調査結果.基本項目.ハッシュID,
          症例.基本情報.年齢,
          Lymph_met, Brain_met, Lung_met, Bone_met, Liver_met, Other_met,
          EP_option, EP_treat, YoungOld,
          症例.基本情報.性別.名称.,
          症例.基本情報.がん種.OncoTree.,
          症例.基本情報.がん種.OncoTree..名称.,
          症例.基本情報.がん種.OncoTree.LEVEL1.,
          症例.検体情報.パネル.名称.,
          症例.背景情報.ECOG.PS.名称.,
          症例.背景情報.喫煙歴有無.名称.,
          症例.背景情報.アルコール多飲有無.名称.,
          症例.背景情報.重複がん有無.異なる臓器..名称.,
          症例.背景情報.多発がん有無.同一臓器..名称.,
          HER2_IHC, MSI_PCR, MMR_IHC,
          final_observe,
          censor,
          CTx_lines_before_CGP,
          pre_CGP_best_RECIST,
          treat_group, treat_group_2, treat_group_3,
          time_enroll_final,
          time_palliative_final,
          time_palliative_enroll,
          time_2L_final,
          time_diagnosis_enroll,
          time_diagnosis_final,
          time_2L_enroll
        )

      if (input$HER2 == "No") {
        Data_case_target = Data_case_target %>% dplyr::select(-HER2_IHC)
      } else {
        Data_case_target = Data_case_target %>% dplyr::filter(HER2_IHC != "Unknown")
      }
      if (input$MSI == "No") {
        Data_case_target = Data_case_target %>% dplyr::select(-MSI_PCR)
      } else {
        Data_case_target = Data_case_target %>% dplyr::filter(MSI_PCR != "Unknown")
      }
      if (input$MMR == "No") {
        Data_case_target = Data_case_target %>% dplyr::select(-MMR_IHC)
      } else {
        Data_case_target = Data_case_target %>% dplyr::filter(MMR_IHC != "Unknown")
      }

      Data_case_target$Cancers = Data_case_target$症例.基本情報.がん種.OncoTree.
      incProgress(1 / 13)

      Data_drug = Data_drug_raw()
      OUTPUT_DATA$figure_surv_CTx_Data_drug_control = Data_drug

      Data_case_target = Data_case_target %>%
        dplyr::distinct(.keep_all = TRUE, C.CAT調査結果.基本項目.ハッシュID)

      Data_MAF = Data_report() %>%
        dplyr::filter(
          !stringr::str_detect(Hugo_Symbol, ",") &
            Hugo_Symbol != "" &
            Evidence_level %in% c("", "A", "B", "C", "D", "E", "F") &
            Variant_Classification != "expression"
        ) %>%
        dplyr::arrange(desc(Evidence_level)) %>%
        dplyr::distinct(Tumor_Sample_Barcode, Hugo_Symbol, Start_Position, .keep_all = TRUE)

      if (length(Data_MAF[Data_MAF$TMB > 30, ]$TMB) > 0) {
        Data_MAF[Data_MAF$TMB > 30, ]$TMB = 30
      }

      if (nrow(Data_case_target) > 0) {
        Data_report_TMB = Data_MAF %>%
          dplyr::filter(Tumor_Sample_Barcode %in% Data_case_target$C.CAT調査結果.基本項目.ハッシュID) %>%
          dplyr::select(Tumor_Sample_Barcode, TMB) %>%
          dplyr::distinct(Tumor_Sample_Barcode, .keep_all = TRUE)

        colnames(Data_report_TMB) = c("C.CAT調査結果.基本項目.ハッシュID", "TMB")
        Data_case_target = dplyr::left_join(Data_case_target, Data_report_TMB, by = "C.CAT調査結果.基本項目.ハッシュID")

        Data_MAF_target = Data_MAF %>%
          dplyr::filter(Tumor_Sample_Barcode %in% Data_case_target$C.CAT調査結果.基本項目.ハッシュID)

        if (input$patho == "Only pathogenic muts") {
          Data_MAF_target = Data_MAF_target %>% dplyr::filter(Evidence_level == "F")
        }

        Data_MAF_target = Data_MAF_target %>%
          dplyr::distinct(Tumor_Sample_Barcode, Hugo_Symbol, .keep_all = TRUE)

        if (length(Data_case_target$censor[is.na(Data_case_target$censor)]) > 0) {
          Data_case_target$censor[is.na(Data_case_target$censor)] = 0
        }
        incProgress(1 / 13)

        Data_case_target$time_pre = Data_case_target$time_diagnosis_enroll
        Data_case_target$time_all = Data_case_target$time_diagnosis_final
        Data_case_target = Data_case_target %>% dplyr::filter(time_pre > 0)

        adjustment = TRUE
        OUTPUT_DATA$figure_surv_CTx_adjustment_control = adjustment

        Data_case_target = Data_case_target %>% dplyr::filter(
          !is.na(time_enroll_final) & is.finite(time_enroll_final) & time_enroll_final > 0 &
            !is.na(time_pre) & is.finite(time_pre) &
            !is.na(censor) & is.finite(censor)
        )

        Data_survival = Data_case_target

        Data_MAF_target = Data_MAF_target %>%
          dplyr::filter(Tumor_Sample_Barcode %in% Data_survival$C.CAT調査結果.基本項目.ハッシュID)

        OUTPUT_DATA$figure_surv_CTx_Data_MAF_target_control = Data_MAF_target
        incProgress(1 / 13)

        Data_survival_interactive = Data_survival
        OUTPUT_DATA$figure_surv_CTx_Data_survival_interactive_control = Data_survival_interactive

        candidate_genes = sort(unique(c(
          Data_MAF_target$Hugo_Symbol,
          paste0(input$special_gene, "_", input$special_gene_mutation_1_name),
          paste0(input$special_gene, "_", input$special_gene_mutation_2_name),
          paste0(input$special_gene, "_NOS")
        )))
        candidate_genes = candidate_genes[!candidate_genes %in% c("", "_", "_NOS", paste0(input$special_gene, "_"))]
        candidate_genes = candidate_genes[!is.na(candidate_genes)]

        Top_gene = unique(names(sort(table(Data_MAF_target$Hugo_Symbol), decreasing = TRUE)))
        Top_gene = Top_gene[Top_gene %in% candidate_genes]

        candidate_drugs = sort(unique(c(Data_drug$Drug)))

        OUTPUT_DATA$figure_surv_interactive_Top_gene_control = Top_gene
        OUTPUT_DATA$figure_surv_interactive_candidate_genes_control = candidate_genes
        OUTPUT_DATA$figure_surv_interactive_candidate_drugs_control = candidate_drugs[!is.na(candidate_drugs)]
        OUTPUT_DATA$figure_surv_interactive_candidate_Age_control = sort(unique(Data_survival_interactive$YoungOld))
        OUTPUT_DATA$figure_surv_interactive_candidate_Sex_control = sort(unique(Data_survival_interactive$症例.基本情報.性別.名称.))
        OUTPUT_DATA$figure_surv_interactive_candidate_Histology_control = sort(unique(Data_survival_interactive$Cancers))
      }

      incProgress(1 / 13)
    })
  })

  rm(analysis_env)
  gc()
}

# =========================================================================
# 2) Survival Plot: Replaced with Proposed Age-stratified matching via weights
# =========================================================================
output$figure_survival_CTx_interactive_1_control = renderPlot({
  req(OUTPUT_DATA$figure_surv_CTx_Data_survival_interactive_control,
      OUTPUT_DATA$figure_surv_CTx_Data_MAF_target_control,
      OUTPUT_DATA$figure_surv_CTx_Data_drug_control)

  # Force reactivity
  lapply(1:2, function(i) {
    prefix <- paste0("gene_survival_interactive_", i, "_")
    list(
      input[[paste0(prefix, "P_1_control")]],
      input[[paste0(prefix, "P_2_control")]],
      input[[paste0(prefix, "W_control")]],
      input[[paste0(prefix, "A_control")]],
      input[[paste0(prefix, "S_control")]],
      input[[paste0(prefix, "H_control")]],
      input[[paste0(prefix, "D_control")]],
      input[[paste0(prefix, "ND_control")]]
    )
  })
  input$registry_cancer_type
  input$survival_data_source

  Data_survival_interactive = OUTPUT_DATA$figure_surv_CTx_Data_survival_interactive_control
  Data_MAF_target = OUTPUT_DATA$figure_surv_CTx_Data_MAF_target_control
  Data_drug = OUTPUT_DATA$figure_surv_CTx_Data_drug_control

  # --- Age parsing (robust) + age_class
  Data_survival_interactive <- Data_survival_interactive %>%
    dplyr::mutate(
      age_num = extract_age_num(症例.基本情報.年齢),
      age_class = make_age_class_6grp(age_num)
    )

  # --- Build macro reference survival
  ref_surv_list <- build_ref_surv_list_realdata(input)

  # Ensure keys exist for the age classes we generate
  for (ag in unique(Data_survival_interactive$age_class)) {
    if (!ag %in% names(ref_surv_list)) {
      if ("全年齢" %in% names(ref_surv_list)) ref_surv_list[[ag]] <- ref_surv_list[["全年齢"]]
      else ref_surv_list[[ag]] <- ref_surv_list[[1]]
    }
  }

  # --- Calculate Matching Weights (Proposed Pipeline)
  Data_survival_interactive <- calculate_obs_matching_weights_real(Data_survival_interactive, ref_surv_list)

  # ========================================================
  # Extract Group IDs
  # ========================================================
  extract_group_ids <- function(group_num) {
    IDs <- unique(Data_survival_interactive$C.CAT調査結果.基本項目.ハッシュID)
    input_prefix <- paste0("gene_survival_interactive_", group_num, "_")

    p1_input <- input[[paste0(input_prefix, "P_1_control")]]
    if (!all(is.null(p1_input))) {
      IDs <- intersect(IDs, (Data_MAF_target %>% dplyr::filter(Hugo_Symbol %in% p1_input))$Tumor_Sample_Barcode)
    }

    p2_input <- input[[paste0(input_prefix, "P_2_control")]]
    if (!all(is.null(p2_input))) {
      IDs <- intersect(IDs, (Data_MAF_target %>% dplyr::filter(Hugo_Symbol %in% p2_input))$Tumor_Sample_Barcode)
    }

    w_input <- input[[paste0(input_prefix, "W_control")]]
    if (!all(is.null(w_input))) {
      IDs <- setdiff(IDs, (Data_MAF_target %>% dplyr::filter(Hugo_Symbol %in% w_input))$Tumor_Sample_Barcode)
    }

    clinical_filters <- list(A_control = "YoungOld", S_control = "症例.基本情報.性別.名称.", H_control = "Cancers")
    for (filter_key in names(clinical_filters)) {
      filter_input <- input[[paste0(input_prefix, filter_key)]]
      if (!all(is.null(filter_input))) {
        column_name <- clinical_filters[[filter_key]]
        filter_expr <- paste0(column_name, " %in% filter_input")
        IDs <- intersect(IDs, (Data_survival_interactive %>% dplyr::filter(!!rlang::parse_expr(filter_expr)))$C.CAT調査結果.基本項目.ハッシュID)
      }
    }

    d_input <- input[[paste0(input_prefix, "D_control")]]
    if (!all(is.null(d_input))) IDs <- intersect(IDs, (Data_drug %>% dplyr::filter(Drug %in% d_input))$ID)

    nd_input <- input[[paste0(input_prefix, "ND_control")]]
    if (!all(is.null(nd_input))) IDs <- setdiff(IDs, (Data_drug %>% dplyr::filter(Drug %in% nd_input))$ID)

    IDs
  }

  ID_1 <- extract_group_ids(1)
  ID_2 <- extract_group_ids(2)

  Data_survival_1 <- Data_survival_interactive %>%
    dplyr::filter(C.CAT調査結果.基本項目.ハッシュID %in% ID_1) %>%
    dplyr::mutate(Group = "1")

  Data_survival_2 <- Data_survival_interactive %>%
    dplyr::filter(C.CAT調査結果.基本項目.ハッシュID %in% ID_2) %>%
    dplyr::mutate(Group = "2")

  Data_EP0 <- Data_survival_interactive %>%
    dplyr::filter(EP_treat == 0) %>%
    dplyr::mutate(Group = "EP0")

  Data_model <- dplyr::bind_rows(Data_EP0, Data_survival_1, Data_survival_2) %>%
    dplyr::mutate(
      Sex = as.factor(`症例.基本情報.性別.名称.`),
      Cancers = as.factor(Cancers)
    ) %>%
    dplyr::filter(time_pre > 0, time_all > time_pre)

  shiny::validate(shiny::need(nrow(Data_model) > 0, "No valid data to fit the model."))
  Data_model$Group <- factor(Data_model$Group, levels = c("EP0", "1", "2"))

  shiny::validate(shiny::need(nrow(Data_model %>% dplyr::filter(Group == "EP0")) >= 5,
                              "Not enough EP_treat == 0 patients to establish the reference baseline."))

  # ========================================================
  # Plot using survival_compare_and_plot_CTx directly with weights
  # ========================================================
  survival_compare_and_plot_CTx(
    data = Data_model,
    time_var1 = "time_pre",
    time_var2 = "time_all",
    status_var = "censor",     # censor==1 means death event
    group_var = "Group",
    plot_title = "Standardized OS (Proposed: Age-stratified matching via weights)",
    adjustment = FALSE,
    color_var_surv_CTx_1 = "diagnosis",
    weights_var = "iptw"
  )
})

# =========================================================================
# 3) Forest Plot: Replaced with Proposed Weighted AFT on Observed OS
# =========================================================================
output$forest_plot_multivariate = renderPlot({
  req(OUTPUT_DATA$figure_surv_CTx_Data_survival_interactive_control,
      OUTPUT_DATA$figure_surv_CTx_Data_MAF_target_control)

  shiny::validate(shiny::need(requireNamespace("patchwork", quietly = TRUE),
                              "Please install the 'patchwork' package. Run: install.packages('patchwork')"))

  Data_survival_interactive = OUTPUT_DATA$figure_surv_CTx_Data_survival_interactive_control
  Data_MAF_target = OUTPUT_DATA$figure_surv_CTx_Data_MAF_target_control
  selected_genes <- input$gene_survival_interactive_1_P_1_control_forest

  input$registry_cancer_type
  input$survival_data_source

  # --- Age parsing (robust) + age_class
  Data_survival_interactive <- Data_survival_interactive %>%
    dplyr::mutate(
      age_num = extract_age_num(症例.基本情報.年齢),
      age_class = make_age_class_6grp(age_num),
      Sex = as.factor(`症例.基本情報.性別.名称.`),
      Cancers = as.character(Cancers)
    )

  # --- Macro reference + Calculate Matching Weights (Proposed Pipeline)
  ref_surv_list <- build_ref_surv_list_realdata(input)
  for (ag in unique(Data_survival_interactive$age_class)) {
    if (!ag %in% names(ref_surv_list)) {
      if ("全年齢" %in% names(ref_surv_list)) ref_surv_list[[ag]] <- ref_surv_list[["全年齢"]]
      else ref_surv_list[[ag]] <- ref_surv_list[[1]]
    }
  }
  Data_survival_interactive <- calculate_obs_matching_weights_real(Data_survival_interactive, ref_surv_list)

  # --- Prepare modeling cohort
  Data_forest <- Data_survival_interactive %>%
    dplyr::filter(is.finite(time_pre), is.finite(time_all), is.finite(censor), time_all > time_pre) %>%
    dplyr::mutate(
      Cancers = as.character(Cancers)
    )

  shiny::validate(shiny::need(nrow(Data_forest) > 50, "Not enough valid cases for multivariate analysis."))
  total_cohort_N <- nrow(Data_forest)

  # Sex reference = most frequent
  if (length(levels(Data_forest$Sex)) > 1) {
    top_sex <- names(sort(table(Data_forest$Sex), decreasing = TRUE))[1]
    Data_forest$Sex <- relevel(Data_forest$Sex, ref = top_sex)
  }

  # Lump rare cancers (<50) into Others, set ref = most frequent
  cancer_counts <- table(Data_forest$Cancers)
  rare_cancers <- names(cancer_counts)[cancer_counts < 50]
  Data_forest$Cancers[Data_forest$Cancers %in% rare_cancers] <- "Others"
  Data_forest$Cancers <- as.factor(Data_forest$Cancers)
  top_cancer <- names(sort(table(Data_forest$Cancers), decreasing = TRUE))[1]
  Data_forest$Cancers <- relevel(Data_forest$Cancers, ref = top_cancer)

  # --- Add selected genes as covariates (Mutated vs WT)
  valid_gene_covariates <- c()
  if (length(selected_genes) > 0) {
    for (gene in selected_genes) {
      mutated_IDs <- (Data_MAF_target %>% dplyr::filter(Hugo_Symbol == gene))$Tumor_Sample_Barcode
      col_name <- paste0("Gene_", make.names(gene))

      Data_forest[[col_name]] <- factor(
        ifelse(Data_forest$C.CAT調査結果.基本項目.ハッシュID %in% mutated_IDs, "Mutated", "WT"),
        levels = c("WT", "Mutated")
      )
      if (length(unique(Data_forest[[col_name]])) > 1) valid_gene_covariates <- c(valid_gene_covariates, col_name)
    }
  }

  # --- Build covariates
  base_covariates <- c("age_num", "Sex", "Cancers")
  all_covariates <- c()
  for (v in base_covariates) {
    if (v %in% names(Data_forest) && length(unique(stats::na.omit(Data_forest[[v]]))) > 1) all_covariates <- c(all_covariates, v)
  }
  all_covariates <- c(all_covariates, valid_gene_covariates)
  shiny::validate(shiny::need(length(all_covariates) > 0, "No valid covariates found."))

  # ========================================================
  # Proposed (Weighted AFT) Modeling:
  # Fit single synchronized model on observed OS (time_all)
  # ========================================================
  form_final <- as.formula(paste("survival::Surv(time_all, censor) ~", paste(all_covariates, collapse = " + ")))
  fit_forest <- tryCatch(
    flexsurv::flexsurvreg(form_final, data = Data_forest, weights = iptw, dist = "llogis"),
    error = function(e) NULL
  )

  shiny::validate(shiny::need(!is.null(fit_forest), "Final synchronized multivariate model failed to converge."))

  # ========================================================
  # Extract results for forest plot (TR = exp(beta))
  # ========================================================
  res <- fit_forest$res
  cov_names <- rownames(res)[!rownames(res) %in% c("shape", "scale")]

  plot_data <- data.frame(
    Variable_Raw = cov_names,
    Estimate = exp(res[cov_names, "est"]),
    Lower = exp(res[cov_names, "L95%"]),
    Upper = exp(res[cov_names, "U95%"]),
    Category = NA, Label = NA, Positive_N = NA,
    Total_N = total_cohort_N, Text_CI = NA
  )

  # Age per +10 years
  for (i in seq_len(nrow(plot_data))) {
    var <- plot_data$Variable_Raw[i]
    if (grepl("^age_num", var)) {
      plot_data$Estimate[i] <- exp(res[var, "est"] * 10)
      plot_data$Lower[i]  <- exp(res[var, "L95%"] * 10)
      plot_data$Upper[i]  <- exp(res[var, "U95%"] * 10)
    }
  }
  plot_data$Text_CI <- sprintf("%.2f (%.2f-%.2f)", plot_data$Estimate, plot_data$Lower, plot_data$Upper)

  for (i in seq_len(nrow(plot_data))) {
    var <- plot_data$Variable_Raw[i]

    if (grepl("^age_num", var)) {
      plot_data$Label[i] <- "Age (per +10 years)"
      plot_data$Category[i] <- "Demographics"

    } else if (grepl("^Sex", var)) {
      level_name <- sub("^Sex", "", var)
      plot_data$Label[i] <- paste0("Sex: ", level_name, " (vs ", levels(Data_forest$Sex)[1], ")")
      plot_data$Category[i] <- "Demographics"
      plot_data$Positive_N[i] <- sum(Data_forest$Sex == level_name, na.rm = TRUE)

    } else if (grepl("^Cancers", var)) {
      cancer_type <- sub("^Cancers", "", var)
      plot_data$Label[i] <- paste0("Histology: ", cancer_type, " (vs ", levels(Data_forest$Cancers)[1], ")")
      plot_data$Category[i] <- "Histology"
      plot_data$Positive_N[i] <- sum(Data_forest$Cancers == cancer_type, na.rm = TRUE)

    } else if (grepl("^Gene_", var)) {
      gene_name <- sub("^Gene_", "", var)
      gene_name <- sub("Mutated$", "", gene_name)
      col_name <- paste0("Gene_", gene_name)
      plot_data$Label[i] <- paste0("Gene: ", gene_name, " (Mutated vs WT)")
      plot_data$Category[i] <- "Genomics"
      if (col_name %in% names(Data_forest)) {
        plot_data$Positive_N[i] <- sum(Data_forest[[col_name]] == "Mutated", na.rm = TRUE)
      } else {
        plot_data$Positive_N[i] <- NA
      }

    } else {
      plot_data$Label[i] <- var
      plot_data$Category[i] <- "Others"
    }
  }

  plot_data$Category <- factor(plot_data$Category, levels = c("Genomics", "Demographics", "Histology", "Others"))
  plot_data <- plot_data[order(plot_data$Category, plot_data$Estimate), ]
  plot_data$Label <- factor(plot_data$Label, levels = plot_data$Label)
  plot_data$Text_PosN <- ifelse(is.na(plot_data$Positive_N), "-", as.character(plot_data$Positive_N))

  # ========================================================
  # Plot
  # ========================================================
  p_left <- ggplot(plot_data, aes(x = Estimate, y = Label, color = Category)) +
    geom_errorbarh(aes(xmin = Lower, xmax = Upper), height = 0.3, linewidth = 0.8) +
    geom_point(size = 4) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "#7f8c8d", linewidth = 1) +
    scale_x_log10(breaks = c(0.1, 0.2, 0.5, 1, 2, 5, 10)) +
    scale_color_manual(values = c("Genomics" = "#e74c3c", "Demographics" = "#2980b9", "Histology" = "#27ae60", "Others" = "#8e44ad")) +
    theme_minimal(base_size = 14) +
    labs(x = "Time Ratio (Log Scale, 95% CI)", y = "") +
    theme(
      panel.grid.minor = element_blank(),
      axis.text.y = element_text(face = "bold", size = 12),
      axis.title.x = element_text(margin = margin(t = 10)),
      legend.position = "none"
    )

  p_right <- ggplot(plot_data, aes(y = Label)) +
    geom_text(aes(x = 0, label = Total_N), hjust = 0.5, size = 4.5, color = "#2c3e50") +
    geom_text(aes(x = 1, label = Text_PosN), hjust = 0.5, size = 4.5, color = "#2c3e50") +
    geom_text(aes(x = 2, label = Text_CI), hjust = 0.5, size = 4.5, color = "#2c3e50") +
    scale_x_continuous(limits = c(-0.5, 2.5), breaks = c(0, 1, 2),
                       labels = c("Total N\n(Real)", "Positive N\n(Real)", "TR\n(95% CI)")) +
    theme_minimal(base_size = 14) +
    labs(x = "", y = "") +
    theme(
      panel.grid = element_blank(),
      axis.text.y = element_blank(),
      axis.text.x = element_text(face = "bold", color = "#2c3e50"),
      axis.ticks = element_blank()
    )

  final_plot <- p_left + p_right + plot_layout(widths = c(2, 1.2)) +
    plot_annotation(
      title = "Independent Prognostic Impact (Standardized AFT Model)",
      subtitle = "Computed via Age-stratified External OS Matching correcting Immortal Time & Dependent Truncation\nTR > 1: Prolonged Survival | TR < 1: Shortened Survival",
      theme = theme(
        plot.title = element_text(face = "bold", size = 16),
        plot.subtitle = element_text(size = 13, color = "#34495e")
      )
    )

  final_plot
}, height = function() {
  Data_whole <- OUTPUT_DATA$figure_surv_CTx_Data_survival_interactive_control
  if (is.null(Data_whole)) return(400)

  n_genes <- length(input$gene_survival_interactive_1_P_1_control_forest)
  cancer_counts <- table(Data_whole$Cancers)
  n_cancers <- sum(cancer_counts >= 50)
  if (n_cancers <= 1) n_cancers <- 1
  estimated_items <- 1 + 1 + (n_cancers - 1) + n_genes
  max(450, estimated_items * 35 + 200)
})
