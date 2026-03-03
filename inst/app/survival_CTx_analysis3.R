# =========================================================================
# Tamura & Ikegami Model Ver 2.3.2 (Real-Data Port)
#  - censor == 1 means EVENT (death)
#  - Keep output variable names: time_pre, time_all, censor, Group
#  - Replaced downstream analysis/plotting with Proposed Simulation Pipeline:
#    Age-stratified external OS curve matching via weights on T_obs (time_all),
#    followed by weighted AFT modeling on observed data.
# =========================================================================
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
# Tamura & Ikegami Model Ver 2.4.0 (Real-Data Port)
#  - censor == 1 means EVENT (death)
#  - Keep output variable names: time_pre, time_all, censor, Group
#  - Proposed pipeline (same spirit as simulation):
#      (1) Age-stratified external OS curve matching via weights on time_all
#      (2) Weighted AFT (llogis) for conditional Time Ratio (TR=exp(beta))
#  - Outputs:
#      * Survival curve (external target + unweighted KM + weighted KM)
#      * Forest plot (weighted AFT)
# =========================================================================

# =========================================================================
# 0) Helpers
# =========================================================================
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

# ----------------------------
# External reference survival list (1-5y, %)
# ----------------------------
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
      for (ag in age_groups_all) ref_surv_list[[ag]] <- as.numeric(vals)

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

  # hard fallback
  if (length(ref_surv_list) == 0) {
    for (ag in age_groups_all) ref_surv_list[[ag]] <- c(50, 30, 20, 15, 10)
  } else {
    for (ag in age_groups_all) {
      if (is.null(ref_surv_list[[ag]]) || length(ref_surv_list[[ag]]) != 5) {
        ref_surv_list[[ag]] <- ref_surv_list[[1]]
      }
    }
  }

  # sanitize 0-100
  for (ag in names(ref_surv_list)) {
    ref_surv_list[[ag]] <- pmax(0, pmin(100, as.numeric(ref_surv_list[[ag]])))
  }

  ref_surv_list
}

# =========================================================================
# 1) External survival smoothing (llogis) same as simulation
# =========================================================================
fit_llogis_from_S15 <- function(S_y_percent) {
  t_years <- 1:5
  S <- pmax(pmin(as.numeric(S_y_percent) / 100, 0.999), 0.001)

  obj <- function(par) {
    shape <- exp(par[1])
    scale <- exp(par[2])
    S_hat <- 1 / (1 + (t_years / scale)^shape)
    sum((S_hat - S)^2)
  }

  init <- c(log(1.5), log(3))
  opt <- tryCatch(
    optim(init, obj, method = "Nelder-Mead", control = list(maxit = 2000)),
    error = function(e) NULL
  )
  if (is.null(opt) || any(!is.finite(opt$par))) return(list(shape = 1.5, scale = 3.0))
  list(shape = exp(opt$par[1]), scale = exp(opt$par[2]))
}

llogis_surv <- function(t_years, shape, scale) {
  t_years <- pmax(t_years, 1e-8)
  1 / (1 + (t_years / scale)^shape)
}

# =========================================================================
# 2) Proposed weights: Age-stratified curve matching on observed OS (time_all)
#   - This is the real-data analogue of calculate_obs_matching_weights()
# =========================================================================
calculate_obs_matching_weights_real <- function(data, ref_surv_list) {
  if (!("age_class" %in% names(data))) stop("age_class missing")
  if (!all(c("time_all", "censor") %in% names(data))) stop("time_all/censor missing")

  data$iptw <- 1.0

  # 0.5-year grid up to 5 years
  t_points_years <- seq(0.5, 5.0, by = 0.5)
  t_points_days  <- t_points_years * 365.25

  for (ag in unique(data$age_class)) {
    if (!(ag %in% names(ref_surv_list))) next
    ag_idx <- which(data$age_class == ag)
    if (length(ag_idx) < 10) next

    ag_data <- data[ag_idx, , drop = FALSE]

    pars <- fit_llogis_from_S15(ref_surv_list[[ag]])
    S_target <- llogis_surv(t_points_years, pars$shape, pars$scale)
    S_target <- pmax(pmin(S_target, 0.999), 0.001)

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

      S_est <- approx(km$time, km$surv, xout = t_points_days,
                      method = "constant", f = 0, rule = 2)$y
      S_est <- pmax(pmin(ifelse(is.finite(S_est), S_est, 0), 0.999), 0.001)

      sum((S_est - S_target)^2) + 0.001 * sum(theta^2)
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

  data$iptw <- as.numeric(data$iptw)
  data$iptw <- ifelse(is.na(data$iptw) | data$iptw <= 0, 1.0, data$iptw)
  data$iptw <- data$iptw / mean(data$iptw, na.rm = TRUE)

  ess <- (sum(data$iptw, na.rm = TRUE)^2) / sum(data$iptw^2, na.rm = TRUE)
  attr(data, "ESS") <- ess
  data
}

# =========================================================================
# 3) Plot helpers: KM -> df, plus external target curve line
# =========================================================================
km_to_df <- function(km_fit, model_label) {
  s <- summary(km_fit)
  data.frame(
    time_years = s$time / 365.25,
    surv = s$surv,
    strata = if (!is.null(s$strata)) as.character(s$strata) else "All",
    Model = model_label,
    stringsAsFactors = FALSE
  )
}

# =========================================================================
# PATCH: EP0 removed + ggsurvplot + title stats (median + HR with CI)
# =========================================================================

# ---- helper: external target curve (llogis) as data.frame ----
make_external_target_df <- function(ref_surv_list, age_key = "全年齢") {
  if (!(age_key %in% names(ref_surv_list))) age_key <- names(ref_surv_list)[1]
  pars <- fit_llogis_from_S15(ref_surv_list[[age_key]])
  t <- seq(0, 5, by = 0.05)
  data.frame(
    time_years = t,
    surv = llogis_surv(t, pars$shape, pars$scale),
    label = paste0("External target (", age_key, ", llogis)"),
    stringsAsFactors = FALSE
  )
}

# ---- helper: weighted KM median (+95%CI) per group (months) ----
calc_weighted_median_ci_months_by_group <- function(data, group_var = "Group",
                                                    time_var = "time_all", status_var = "censor",
                                                    weight_var = "iptw") {
  shiny::validate(shiny::need(requireNamespace("survminer", quietly = TRUE),
                              "Please install the 'survminer' package. Run: install.packages('survminer')"))

  data <- data %>%
    dplyr::filter(is.finite(.data[[time_var]]), .data[[time_var]] > 0,
                  is.finite(.data[[status_var]]), !is.na(.data[[group_var]]))

  sf <- tryCatch(
    survival::survfit(
      as.formula(paste0("survival::Surv(", time_var, ", ", status_var, ") ~ ", group_var)),
      data = data,
      weights = data[[weight_var]],
      conf.type = "log"   # 安定（好みで "plain" など）
    ),
    error = function(e) NULL
  )
  if (is.null(sf)) return(NULL)

  md <- tryCatch(survminer::surv_median(sf), error = function(e) NULL)
  if (is.null(md) || nrow(md) == 0) return(NULL)

  # strata の "Group=1" みたいな表記を "1" に
  md$Group <- ifelse(grepl("=", md$strata), sub(".*=", "", md$strata), md$strata)

  data.frame(
    Group = md$Group,
    Med_month = md$median / 365.25 * 12,
    L_month   = md$lower  / 365.25 * 12,
    U_month   = md$upper  / 365.25 * 12,
    stringsAsFactors = FALSE
  )
}

# ---- helper: weighted Cox HR (robust SE) between 2 groups ----
calc_weighted_hr <- function(data, group_var = "Group",
                             time_var = "time_all", status_var = "censor",
                             weight_var = "iptw", min_events_per_group = 5) {
  data <- data %>%
    dplyr::filter(is.finite(.data[[time_var]]), .data[[time_var]] > 0,
                  is.finite(.data[[status_var]]),
                  is.finite(.data[[weight_var]]), .data[[weight_var]] > 0,
                  !is.na(.data[[group_var]])) %>%
    dplyr::mutate(
      Group2 = droplevels(as.factor(.data[[group_var]]))
    )

  lv <- levels(data$Group2)
  if (length(lv) != 2) return(NULL)

  # 参照を固定（例： "1" をrefにしたいなら levels順をここで設定）
  # data$Group2 <- relevel(data$Group2, ref = lv[1])

  # イベント数チェック（censor==1がイベント）
  ev_by_g <- tapply(data[[status_var]], data$Group2, function(x) sum(x == 1, na.rm = TRUE))
  if (any(is.na(ev_by_g)) || any(ev_by_g < min_events_per_group)) {
    message(sprintf("[HR] Too few events: %s", paste(names(ev_by_g), ev_by_g, sep="=", collapse=" | ")))
    return(NULL)
  }

  fit <- tryCatch(
    survival::coxph(
      survival::Surv(data[[time_var]], data[[status_var]]) ~ Group2,
      data = data,
      weights = data[[weight_var]],
      robust = TRUE,
      ties = "efron"
    ),
    error = function(e) { message("[HR] coxph error: ", e$message); NULL }
  )
  if (is.null(fit)) return(NULL)

  sm <- summary(fit)
  if (any(!is.finite(sm$coefficients[1, c("exp(coef)")])) ||
      any(!is.finite(sm$conf.int[1, c("lower .95","upper .95")]))) {
    message("[HR] Non-finite HR/CI produced (separation or numerical issue).")
    return(NULL)
  }

  list(
    ref  = lv[1],
    comp = lv[2],
    HR = unname(sm$coefficients[1, "exp(coef)"]),
    L  = unname(sm$conf.int[1, "lower .95"]),
    U  = unname(sm$conf.int[1, "upper .95"])
  )
}

# ---- helper: build multi-line title (wrapped) ----
build_surv_title <- function(med_df, hr_obj, ess, max_width = 56) {
  # medians with CI
  if (is.null(med_df) || nrow(med_df) == 0) {
    med_txt <- "Median OS (months): N/A"
  } else {
    med_df <- med_df %>% dplyr::arrange(Group)
    med_pairs <- paste0(
      med_df$Group, "=",
      sprintf("%.1f", med_df$Med_month),
      " (", sprintf("%.1f", med_df$L_month), "–", sprintf("%.1f", med_df$U_month), ")"
    )
    med_txt <- paste0("Median OS (months, 95%CI): ", paste(med_pairs, collapse = " | "))
  }

  hr_txt <- "HR: N/A"
  if (!is.null(hr_obj)) {
    hr_txt <- sprintf("HR (%s vs %s)=%.2f (95%%CI %.2f–%.2f)",
                      hr_obj$comp, hr_obj$ref, hr_obj$HR, hr_obj$L, hr_obj$U)
  }

  ess_txt <- if (is.finite(ess)) sprintf("ESS=%.1f", ess) else "ESS=N/A"

  title1 <- paste(strwrap("Standardized OS (calibrated weighted KM; ggsurvplot)", width = max_width), collapse = "\n")
  title2 <- paste(strwrap(paste(med_txt, hr_txt, ess_txt, sep = "  |  "), width = max_width), collapse = "\n")
  paste(title1, title2, sep = "\n")
}

# =========================================================================
# Survival curve output (EP0 removed, ggsurvplot, title stats)
# =========================================================================
output$figure_survival_CTx_interactive_1_control <- renderPlot({
  req(OUTPUT_DATA$figure_surv_CTx_Data_survival_interactive_control,
      OUTPUT_DATA$figure_surv_CTx_Data_MAF_target_control,
      OUTPUT_DATA$figure_surv_CTx_Data_drug_control)

  # Touch inputs for reactivity (your original)
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

  Data_survival_interactive <- OUTPUT_DATA$figure_surv_CTx_Data_survival_interactive_control
  Data_MAF_target <- OUTPUT_DATA$figure_surv_CTx_Data_MAF_target_control
  Data_drug <- OUTPUT_DATA$figure_surv_CTx_Data_drug_control

  # --- Age parsing + age_class
  Data_survival_interactive <- Data_survival_interactive %>%
    dplyr::mutate(
      age_num = extract_age_num(症例.基本情報.年齢),
      age_class = make_age_class_6grp(age_num)
    )

  # --- External reference + ensure keys
  ref_surv_list <- build_ref_surv_list_realdata(input)
  for (ag in unique(Data_survival_interactive$age_class)) {
    if (!ag %in% names(ref_surv_list)) {
      if ("全年齢" %in% names(ref_surv_list)) ref_surv_list[[ag]] <- ref_surv_list[["全年齢"]]
      else ref_surv_list[[ag]] <- ref_surv_list[[1]]
    }
  }

  # --- Calibrated weights on observed OS (time_all)
  Data_survival_interactive <- calculate_obs_matching_weights_real(Data_survival_interactive, ref_surv_list)
  ESS <- attr(Data_survival_interactive, "ESS")

  # ========================================================
  # Extract Group IDs (same as your code)
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

  # ---- EP0 removed: model cohort = Group 1 & 2 only ----
  Data_model <- dplyr::bind_rows(Data_survival_1, Data_survival_2) %>%
    dplyr::mutate(
      Sex = as.factor(`症例.基本情報.性別.名称.`),
      Cancers = as.factor(Cancers)
    ) %>%
    dplyr::filter(is.finite(time_all), is.finite(censor), time_all > 0)

  shiny::validate(shiny::need(nrow(Data_model) > 30, "Not enough valid data to plot."))
  Data_model$Group <- factor(Data_model$Group, levels = c("1", "2"))
  shiny::validate(shiny::need(nlevels(Data_model$Group) == 2, "Need exactly 2 groups (1 and 2) for HR display."))

  # --- KM fits (unweighted + calibrated weighted)
  km_unw <- survival::survfit(survival::Surv(time_all, censor) ~ Group, data = Data_model)
  km_w   <- survival::survfit(survival::Surv(time_all, censor) ~ Group, data = Data_model, weights = Data_model$iptw)

  # --- stats for title
  med_df <- calc_weighted_median_ci_months_by_group(Data_model, weight_var = "iptw")
  hr_obj <- calc_weighted_hr(Data_model, weight_var = "iptw")
  title_txt <- build_surv_title(med_df, hr_obj, ESS, max_width = 56)

  gp <- survminer::ggsurvplot(
    fit = km_w,
    data = Data_model,
    conf.int = FALSE,
    risk.table = TRUE,
    risk.table.height = 0.28,
    risk.table.y.text = FALSE,
    censor = TRUE,
    xlim = c(0, 5 * 365.25),
    break.time.by = 365.25,
    xlab = "Time from diagnosis (years)",
    ylab = "Overall survival",
    xscale = "d_y",
    ggtheme = ggplot2::theme_minimal(base_size = 15),
    legend = "bottom",
    legend.title = ""
  )

  ext_df <- make_external_mixture_target_df(
    ref_surv_list = ref_surv_list,
    age_class_vec = Data_model$age_class,   # <- ここが重要（描画対象の年齢構成）
    t_max_years = 5,
    dt_years = 0.05,
    label_prefix = "External mixture target"
  )
  # ---- overlay external target line with unambiguous legend ----
  gp$plot <- gp$plot +
    ggplot2::geom_line(
      data = ext_df,
      ggplot2::aes(x = time_years * 365.25, y = surv, linetype = label),
      linewidth = 1.1,
      inherit.aes = FALSE
    ) +
    ggplot2::scale_linetype_manual(values = c("External Cancer Registry Data" = "dashed")) +
    ggplot2::ggtitle(title_txt) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 14, lineheight = 1.05),
      plot.margin = ggplot2::margin(t = 10, r = 15, b = 5, l = 10),
      legend.position = "bottom",
      legend.box = "vertical"
    )

  gp$table <- gp$table +
    ggplot2::theme(
      plot.margin = ggplot2::margin(t = 0, r = 15, b = 10, l = 10)
    )

  # ---- IMPORTANT FIX: arrange_ggsurvplots needs a LIST of ggsurvplot objects ----
  arranged <- survminer::arrange_ggsurvplots(
    list(gp),
    ncol = 1, nrow = 1,
    heights = c(0.72, 0.28)
  )

  # In renderPlot, explicitly print the arranged grob
  print(arranged)
})
# =========================================================================
# 7) Forest plot output: weighted AFT on observed OS (time_all)
# =========================================================================
output$forest_plot_multivariate <- renderPlot({
  req(OUTPUT_DATA$figure_surv_CTx_Data_survival_interactive_control,
      OUTPUT_DATA$figure_surv_CTx_Data_MAF_target_control)

  Data_survival_interactive <- OUTPUT_DATA$figure_surv_CTx_Data_survival_interactive_control
  Data_MAF_target <- OUTPUT_DATA$figure_surv_CTx_Data_MAF_target_control
  selected_genes <- input$gene_survival_interactive_1_P_1_control_forest

  input$registry_cancer_type
  input$survival_data_source

  # --- Age parsing + age_class + factors
  Data_survival_interactive <- Data_survival_interactive %>%
    dplyr::mutate(
      age_num = extract_age_num(症例.基本情報.年齢),
      age_class = make_age_class_6grp(age_num),
      Sex = as.factor(`症例.基本情報.性別.名称.`),
      Cancers = as.character(Cancers)
    )

  # --- External reference + weights
  ref_surv_list <- build_ref_surv_list_realdata(input)
  for (ag in unique(Data_survival_interactive$age_class)) {
    if (!ag %in% names(ref_surv_list)) {
      if ("全年齢" %in% names(ref_surv_list)) ref_surv_list[[ag]] <- ref_surv_list[["全年齢"]]
      else ref_surv_list[[ag]] <- ref_surv_list[[1]]
    }
  }
  Data_survival_interactive <- calculate_obs_matching_weights_real(Data_survival_interactive, ref_surv_list)

  # --- Modeling cohort
  Data_forest <- Data_survival_interactive %>%
    dplyr::filter(is.finite(time_pre), is.finite(time_all), is.finite(censor),
                  time_all > time_pre, time_all > 0) %>%
    dplyr::mutate(Cancers = as.character(Cancers))

  shiny::validate(shiny::need(nrow(Data_forest) > 50, "Not enough valid cases for multivariate analysis."))

  # Sex reference = most frequent
  if (length(levels(Data_forest$Sex)) > 1) {
    top_sex <- names(sort(table(Data_forest$Sex), decreasing = TRUE))[1]
    Data_forest$Sex <- relevel(Data_forest$Sex, ref = top_sex)
  }

  # Lump rare cancers (<50) into Others, ref = most frequent
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

  # --- Covariates
  base_covariates <- c("age_num", "Sex", "Cancers")
  covariates <- c()
  for (v in base_covariates) {
    if (v %in% names(Data_forest) && length(unique(stats::na.omit(Data_forest[[v]]))) > 1) covariates <- c(covariates, v)
  }
  covariates <- c(covariates, valid_gene_covariates)
  shiny::validate(shiny::need(length(covariates) > 0, "No valid covariates found."))

  # --- Forest plot (weighted AFT)
  build_forest_plot_weighted_aft(
    data_forest = Data_forest,
    covariates = covariates,
    time_var = "time_all",
    status_var = "censor",
    weight_var = "iptw"
  )
}, height = function() {
  Data_whole <- OUTPUT_DATA$figure_surv_CTx_Data_survival_interactive_control
  if (is.null(Data_whole)) return(450)

  n_genes <- length(input$gene_survival_interactive_1_P_1_control_forest)
  cancer_counts <- table(Data_whole$Cancers)
  n_cancers <- sum(cancer_counts >= 50)
  if (n_cancers <= 1) n_cancers <- 1
  estimated_items <- 1 + 1 + (n_cancers - 1) + n_genes  # Age + Sex + cancers + genes
  max(500, estimated_items * 35 + 220)
})

build_forest_plot_weighted_aft <- function(data_forest, covariates,
                                           time_var = "time_all", status_var = "censor",
                                           weight_var = "iptw") {
  shiny::validate(shiny::need(requireNamespace("patchwork", quietly = TRUE),
                              "Please install the 'patchwork' package. Run: install.packages('patchwork')"))

  form_final <- as.formula(
    paste0("survival::Surv(", time_var, ", ", status_var, ") ~ ", paste(covariates, collapse = " + "))
  )

  fit <- tryCatch(
    flexsurv::flexsurvreg(form_final, data = data_forest, weights = data_forest[[weight_var]], dist = "llogis"),
    error = function(e) NULL
  )
  shiny::validate(shiny::need(!is.null(fit), "Weighted AFT (llogis) failed to converge."))

  res <- fit$res
  cov_names <- rownames(res)[!rownames(res) %in% c("shape", "scale")]

  total_N <- nrow(data_forest)

  plot_data <- data.frame(
    Variable_Raw = cov_names,
    Estimate = exp(res[cov_names, "est"]),
    Lower = exp(res[cov_names, "L95%"]),
    Upper = exp(res[cov_names, "U95%"]),
    Category = "Others",
    Label = cov_names,
    Positive_N = NA_real_,
    Total_N = total_N,
    stringsAsFactors = FALSE
  )

  # Age per +10 years if present
  for (i in seq_len(nrow(plot_data))) {
    var <- plot_data$Variable_Raw[i]
    if (grepl("^age_num", var)) {
      plot_data$Estimate[i] <- exp(res[var, "est"] * 10)
      plot_data$Lower[i]    <- exp(res[var, "L95%"] * 10)
      plot_data$Upper[i]    <- exp(res[var, "U95%"] * 10)
      plot_data$Label[i]    <- "Age (per +10 years)"
      plot_data$Category[i] <- "Demographics"
    }
    if (grepl("^Sex", var)) plot_data$Category[i] <- "Demographics"
    if (grepl("^Cancers", var)) plot_data$Category[i] <- "Histology"
    if (grepl("^Gene_", var)) plot_data$Category[i] <- "Genomics"
  }

  # Pretty labels + Positive_N
  if ("Sex" %in% covariates && "Sex" %in% names(data_forest)) {
    ref_sex <- levels(data_forest$Sex)[1]
    for (i in seq_len(nrow(plot_data))) {
      var <- plot_data$Variable_Raw[i]
      if (grepl("^Sex", var)) {
        lv <- sub("^Sex", "", var)
        plot_data$Label[i] <- paste0("Sex: ", lv, " (vs ", ref_sex, ")")
        plot_data$Positive_N[i] <- sum(data_forest$Sex == lv, na.rm = TRUE)
      }
    }
  }

  if ("Cancers" %in% covariates && "Cancers" %in% names(data_forest)) {
    ref_can <- levels(data_forest$Cancers)[1]
    for (i in seq_len(nrow(plot_data))) {
      var <- plot_data$Variable_Raw[i]
      if (grepl("^Cancers", var)) {
        lv <- sub("^Cancers", "", var)
        plot_data$Label[i] <- paste0("Histology: ", lv, " (vs ", ref_can, ")")
        plot_data$Positive_N[i] <- sum(data_forest$Cancers == lv, na.rm = TRUE)
      }
    }
  }

  for (i in seq_len(nrow(plot_data))) {
    var <- plot_data$Variable_Raw[i]
    if (grepl("^Gene_", var)) {
      gene_nm <- sub("^Gene_", "", var)
      gene_nm <- sub("Mutated$", "", gene_nm)
      col_name <- paste0("Gene_", gene_nm)
      plot_data$Label[i] <- paste0("Gene: ", gene_nm, " (Mutated vs WT)")
      if (col_name %in% names(data_forest)) {
        plot_data$Positive_N[i] <- sum(data_forest[[col_name]] == "Mutated", na.rm = TRUE)
      }
    }
  }

  plot_data$Text_CI <- sprintf("%.2f (%.2f-%.2f)", plot_data$Estimate, plot_data$Lower, plot_data$Upper)
  plot_data$Text_PosN <- ifelse(is.na(plot_data$Positive_N), "-", as.character(plot_data$Positive_N))

  plot_data$Category <- factor(plot_data$Category, levels = c("Genomics", "Demographics", "Histology", "Others"))
  plot_data <- plot_data[order(plot_data$Category, plot_data$Estimate), ]
  plot_data$Label <- factor(plot_data$Label, levels = plot_data$Label)

  p_left <- ggplot2::ggplot(plot_data, ggplot2::aes(x = Estimate, y = Label, color = Category)) +
    ggplot2::geom_errorbarh(ggplot2::aes(xmin = Lower, xmax = Upper), height = 0.3, linewidth = 0.8) +
    ggplot2::geom_point(size = 3.8) +
    ggplot2::geom_vline(xintercept = 1, linetype = "dashed", linewidth = 1) +
    ggplot2::scale_x_log10(breaks = c(0.1, 0.2, 0.5, 1, 2, 5, 10)) +
    ggplot2::theme_minimal(base_size = 14) +
    ggplot2::labs(x = "Time Ratio (log scale, 95% CI)", y = "") +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_text(face = "bold", size = 11),
      legend.position = "none"
    )

  p_right <- ggplot2::ggplot(plot_data, ggplot2::aes(y = Label)) +
    ggplot2::geom_text(ggplot2::aes(x = 0, label = Total_N), hjust = 0.5, size = 4.2) +
    ggplot2::geom_text(ggplot2::aes(x = 1, label = Text_PosN), hjust = 0.5, size = 4.2) +
    ggplot2::geom_text(ggplot2::aes(x = 2, label = Text_CI), hjust = 0.5, size = 4.2) +
    ggplot2::scale_x_continuous(limits = c(-0.5, 2.5), breaks = c(0, 1, 2),
                                labels = c("Total N\n(Real)", "Positive N\n(Real)", "TR\n(95% CI)")) +
    ggplot2::theme_minimal(base_size = 14) +
    ggplot2::labs(x = "", y = "") +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank()
    )

  p_left + p_right + patchwork::plot_layout(widths = c(2.0, 1.25)) +
    patchwork::plot_annotation(
      title = "Independent Prognostic Impact (Calibrated weighted AFT; llogis)",
      subtitle = "TR > 1: prolonged survival | TR < 1: shortened survival\nWeights are learned by age-stratified matching of observed OS to external target OS",
      theme = ggplot2::theme(
        plot.title = ggplot2::element_text(face = "bold", size = 16),
        plot.subtitle = ggplot2::element_text(size = 12)
      )
    )
}

make_external_mixture_target_df <- function(ref_surv_list,
                                            age_class_vec,
                                            t_max_years = 5,
                                            dt_years = 0.05,
                                            label_prefix = "External mixture target") {

  age_class_vec <- as.character(age_class_vec)
  age_class_vec <- age_class_vec[!is.na(age_class_vec) & age_class_vec != "全年齢"]

  # age_classが空なら fallback
  if (length(age_class_vec) == 0) {
    return(make_external_target_df(ref_surv_list, age_key = "全年齢"))
  }

  # mixture weights (proportions in the plotting cohort)
  wtab <- table(age_class_vec)
  pi <- as.numeric(wtab) / sum(wtab)
  names(pi) <- names(wtab)

  # make sure each component exists in ref_surv_list
  valid_keys <- intersect(names(pi), names(ref_surv_list))
  if (length(valid_keys) == 0) {
    return(make_external_target_df(ref_surv_list, age_key = "全年齢"))
  }
  pi <- pi[valid_keys]
  pi <- pi / sum(pi)

  # time grid
  t <- seq(0, t_max_years, by = dt_years)

  # compute mixture survival: S_mix(t) = Σ pi_a * S_a(t)
  S_mix <- rep(0, length(t))
  for (ag in names(pi)) {
    pars <- fit_llogis_from_S15(ref_surv_list[[ag]])
    S_ag <- llogis_surv(t, pars$shape, pars$scale)
    S_mix <- S_mix + pi[ag] * S_ag
  }

  # label includes mixture composition (optional, short)
  mix_txt <- paste0(names(pi), "=", sprintf("%.2f", pi), collapse = ",")
  label <- paste0(label_prefix, " (", mix_txt, ")")

  data.frame(
    time_years = t,
    surv = pmax(pmin(S_mix, 0.999), 0.001),
    label = label,
    stringsAsFactors = FALSE
  )
}
