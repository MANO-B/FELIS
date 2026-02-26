# =========================================================================
# Server Logic for Copula Model (Robust Extraction Version)
# =========================================================================

copula_results <- eventReactive(input$run_copula, {

  shiny::validate(shiny::need(requireNamespace("depend.truncation", quietly = TRUE), "Please install 'depend.truncation' package."))
  shiny::validate(shiny::need(requireNamespace("copula", quietly = TRUE), "Please install 'copula' package."))

  id <- showNotification("Data sampling... (Max N=150 for NPMLE computation)", duration = NULL, type = "message")
  on.exit(removeNotification(id), add = TRUE)

  library(copula)
  library(depend.truncation)
  library(survival)

  N_macro <- input$cop_n
  tau_target <- input$cop_tau
  cens_rate <- input$cop_cens / 100

  # =========================================================================
  # Generate Data
  # =========================================================================
  theta_frank <- tryCatch({ iTau(frankCopula(), tau_target) }, error = function(e) { 5.0 })
  cop <- frankCopula(param = theta_frank, dim = 2)

  uv <- rCopula(N_macro, cop)
  u <- pmax(pmin(uv[, 1], 0.999), 0.001)
  v <- pmax(pmin(uv[, 2], 0.999), 0.001)

  T_true <- qweibull(u, shape = 1.5, scale = 3 * 365.25)
  T1_latent <- qweibull(v, shape = 1.5, scale = 1.5 * 365.25)

  idx_cgp <- which(T_true > T1_latent)
  T_true_obs <- T_true[idx_cgp]
  T1_obs <- T1_latent[idx_cgp]

  T2_true <- T_true_obs - T1_obs
  if (cens_rate > 0) {
    rate_c <- cens_rate / mean(T2_true, na.rm = TRUE)
    C2_obs <- rexp(length(T_true_obs), rate = rate_c)
  } else {
    C2_obs <- rep(Inf, length(T_true_obs))
  }

  T_obs <- T1_obs + pmin(T2_true, C2_obs)
  Event <- ifelse(T2_true <= C2_obs, 1, 0)

  # Jitter to prevent ties crashing NPMLE
  set.seed(Sys.time())
  T1_obs <- T1_obs + runif(length(T1_obs), 0, 0.05)
  T_obs <- T_obs + runif(length(T_obs), 0, 0.05)

  Data_cgp <- data.frame(T1 = T1_obs, T_obs = T_obs, Event = Event)
  shiny::validate(shiny::need(sum(Data_cgp$Event) > 20, "Not enough events generated. Increase N or reduce censoring."))

  max_n_for_copula <- 150
  if (nrow(Data_cgp) > max_n_for_copula) {
    Data_cgp <- Data_cgp[sample(1:nrow(Data_cgp), max_n_for_copula), ]
  }

  # =========================================================================
  # Fit Models
  # =========================================================================
  fit_true <- survfit(Surv(T_true, rep(1, N_macro)) ~ 1)
  fit_naive <- survfit(Surv(T_obs, Event) ~ 1, data = Data_cgp)
  fit_lb <- survfit(Surv(T1, T_obs, Event) ~ 1, data = Data_cgp)

  showNotification("Fitting NPMLE.Frank (Optimizing Likelihood)... Expected time: 3-10 seconds.", duration = 10, type = "warning")

  start_time <- Sys.time()
  fit_copula <- tryCatch({
    depend.truncation::NPMLE.Frank(Data_cgp$T1, Data_cgp$T_obs, Data_cgp$Event)
  }, error = function(e) { e$message })
  end_time <- Sys.time()

  shiny::validate(shiny::need(is.list(fit_copula), paste("Copula model failed. Error:", fit_copula)))
  showNotification(paste("Copula fit completed in", round(as.numeric(difftime(end_time, start_time, units="secs")), 1), "seconds!"), type = "message")

  # =========================================================================
  # Dynamic & Robust Extraction of NPMLE Output
  # =========================================================================
  out_names <- names(fit_copula)

  # Extract Tau safely
  copula_tau <- if ("tau" %in% out_names) fit_copula$tau else if ("Tau" %in% out_names) fit_copula$Tau else NA

  # Extract Survival Vector (Sy, S.y, S_y, S.Y, etc.)
  copula_surv <- NULL
  if ("Sy" %in% out_names) copula_surv <- fit_copula$Sy
  else if ("S.y" %in% out_names) copula_surv <- fit_copula$S.y
  else if ("S_y" %in% out_names) copula_surv <- fit_copula$S_y
  else if ("S.Y" %in% out_names) copula_surv <- fit_copula$S.Y
  else if ("S_Y" %in% out_names) copula_surv <- fit_copula$S_Y

  # Extract Time Vector (y, Y, time, etc.)
  copula_time <- NULL
  if ("y" %in% out_names) copula_time <- fit_copula$y
  else if ("Y" %in% out_names) copula_time <- fit_copula$Y
  else if ("time" %in% out_names) copula_time <- fit_copula$time

  # If time vector is missing but survival exists, assign sorted unique observed times
  if (is.null(copula_time) && !is.null(copula_surv)) {
    unique_T <- sort(unique(Data_cgp$T_obs[Data_cgp$Event == 1]))
    if (length(copula_surv) == length(unique_T)) {
      copula_time <- unique_T
    } else {
      copula_time <- seq(min(Data_cgp$T_obs), max(Data_cgp$T_obs), length.out = length(copula_surv))
    }
  }

  # ** Fail-safe UI Error Message **
  # If extraction still fails, it will print the EXACT variable names the package returned.
  error_msg <- paste("Extraction failed. The package returned these variables:", paste(out_names, collapse = ", "))
  shiny::validate(shiny::need(!is.null(copula_surv) && !is.null(copula_time), error_msg))

  list(
    fit_true = fit_true, fit_naive = fit_naive, fit_lb = fit_lb,
    copula_time = copula_time, copula_surv = copula_surv,
    copula_tau = copula_tau
  )
})

# =========================================================================
# テーブル描画
# =========================================================================
output$copula_result_table <- renderTable({
  res <- copula_results()

  get_median <- function(time_vec, surv_vec) {
    idx <- which(surv_vec <= 0.5)
    if (length(idx) == 0) return(NA)
    return(time_vec[min(idx)])
  }

  med_true <- summary(res$fit_true)$table["median"] / 365.25
  med_naive <- summary(res$fit_naive)$table["median"] / 365.25
  med_lb <- summary(res$fit_lb)$table["median"] / 365.25
  med_copula <- get_median(res$copula_time, res$copula_surv) / 365.25

  data.frame(
    Model = c("1. True Marginal Survival",
              "2. Naive KM (Ignores Truncation)",
              "3. Standard LT-KM (Lynden-Bell)",
              "4. Copula Model (Frank NPMLE)"),
    `Median Survival (Years)` = c(sprintf("%.2f", med_true),
                                  sprintf("%.2f", med_naive),
                                  sprintf("%.2f", med_lb),
                                  ifelse(is.na(med_copula), "Not Reached", sprintf("%.2f", med_copula))),
    `Bias (Years)` = c("-",
                       sprintf("%+.2f", med_naive - med_true),
                       sprintf("%+.2f", med_lb - med_true),
                       ifelse(is.na(med_copula), "-", sprintf("%+.2f", med_copula - med_true)))
  )
}, bordered = TRUE, striped = TRUE, hover = TRUE, align = "c")

# =========================================================================
# プロット描画
# =========================================================================
output$copula_survival_plot <- renderPlot({
  res <- copula_results()
  library(ggplot2)

  # 原点 (0, 1) を追加してカーブがY軸から描画されるようにする関数
  add_origin <- function(df) {
    if (nrow(df) > 0 && df$Time[1] > 0) {
      rbind(data.frame(Time = 0, Surv = 1, Method = df$Method[1]), df)
    } else df
  }

  df_true <- data.frame(Time = res$fit_true$time / 365.25, Surv = res$fit_true$surv, Method = "1. True Survival")
  df_naive <- data.frame(Time = res$fit_naive$time / 365.25, Surv = res$fit_naive$surv, Method = "2. Naive KM")
  df_lb <- data.frame(Time = res$fit_lb$time / 365.25, Surv = res$fit_lb$surv, Method = "3. Standard LT-KM (Lynden-Bell)")
  df_copula <- data.frame(Time = res$copula_time / 365.25, Surv = res$copula_surv, Method = "4. Copula Model (Frank NPMLE)")

  plot_df <- rbind(add_origin(df_true), add_origin(df_naive), add_origin(df_lb), add_origin(df_copula))

  tau_str <- if(is.na(res$copula_tau)) "N/A" else sprintf("%.2f", res$copula_tau)

  ggplot(plot_df, aes(x = Time, y = Surv, color = Method, linetype = Method)) +
    geom_step(size = 1.2) +
    scale_color_manual(values = c("1. True Survival" = "black",
                                  "2. Naive KM" = "#e74c3c",
                                  "3. Standard LT-KM (Lynden-Bell)" = "#f39c12",
                                  "4. Copula Model (Frank NPMLE)" = "#27ae60")) +
    scale_linetype_manual(values = c("1. True Survival" = "solid",
                                     "2. Naive KM" = "dotted",
                                     "3. Standard LT-KM (Lynden-Bell)" = "dashed",
                                     "4. Copula Model (Frank NPMLE)" = "solid")) +
    scale_y_continuous(limits = c(0, 1), labels = scales::percent_format()) +
    scale_x_continuous(limits = c(0, 8)) +
    theme_minimal(base_size = 15) +
    labs(
      title = "Marginal Survival Recovery under Dependent Truncation",
      subtitle = paste0("Copula Model Estimated Kendall's Tau: ", tau_str),
      x = "Time from Diagnosis (Years)",
      y = "Overall Survival Probability"
    ) +
    theme(
      plot.title = element_text(face = "bold"),
      legend.position = "bottom",
      legend.title = element_blank(),
      legend.direction = "vertical"
    )
})
