get_cpu_architecture <- function() {
  arch <- R.version$arch
  if (arch == "aarch64") {
    return("arm64")
  } else if (arch == "x86_64") {
    return("amd64")
  } else {
    return("unknown")
  }
}

QS_READ <- function(nthreads, file, ...) {
  cpu_arch <- get_cpu_architecture()
  if (cpu_arch == "arm64") {
    # Arm64の場合、nthreadsをデフォルト設定
    qs2::qs_read(file=file, nthreads=nthreads, ...)
  } else if(!file.exists("/.dockerenv")){
    # その他のアーキテクチャの場合、nthreadsを1に設定
    qs2::qs_read(file=file, nthreads = 1, ...)
  } else {
    qs2::qs_read(file=file, nthreads = 1, ...)
  }
}
QS_SAVE <- function(nthreads, object, file, ...) {
  cpu_arch <- get_cpu_architecture()
  if (cpu_arch == "arm64") {
    # Arm64の場合、nthreadsをデフォルト設定
    qs2::qs_save(nthreads=nthreads, object=object, file=file, ...)
  } else if(!file.exists("/.dockerenv")){
    # その他のアーキテクチャの場合、nthreadsを1に設定
    qs2::qs_save(nthreads = 1, object=object, file=file, ...)
  } else {
    qs2::qs_save(nthreads = 1, object=object, file=file, ...)
  }
}
Abs = function(x){
  x = as.numeric(x)
  x[is.na(x)]=1
  return(x)
}
get_type_optimized = function(x) {
  if (is.na(x) || x == "") return(character(0))
  unlist(strsplit(x, ";", fixed = TRUE))
}

format_p <- function(x, digits = 1, scientific = TRUE) {
  # ゼロの場合
  zero_case <- sprintf(paste0("0.%0", digits, "d"), 0)

  # 通常表示の場合
  normal_case <- sprintf(paste0("%.", digits, "f"), x)

  # 指数表示の場合
  e <- floor(log10(abs(x)))
  m <- x / (10^e)
  sci_case <- paste0(sprintf(paste0("%.", digits, "f"), m), "×10^", e)

  # 条件に基づいて選択
  result <- ifelse(x == 0, zero_case,
                   ifelse(abs(x) >= 10^(-digits) & abs(x) < 1000, normal_case, sci_case))

  return(result)
}

apply_custom_format <- function(tbl, columns = c("estimate", "conf.low", "conf.high", "p.value", "q.value")) {
  tbl |>
    modify_table_styling(
      columns = columns,
      fmt_fun = function(x) {
        sapply(x, function(val) {
          if (is.na(val)) return(NA_character_)
          format_p(val, digits = 2, scientific = TRUE)
        })
      }
    )
}

fun_zero <- function(a, b){
  if(length(a) != length(b)){
    if(length(a) / length(b) == as.integer(length(a) / length(b))){
      b = rep(b,length(a) / length(b))
    }
  }
  return(ifelse(b == 0, 0, a / b))
}
convert_date <- function(x) {
  case_when(
    x < "1924-02-19" & x > "1901-01-01" ~ paste0("20", str_sub(x, 3, 10)),
    TRUE ~ x
  )
}

gg_empty = function(){
  g = ggplot()
  g = g + geom_blank()
  g = g + ggtitle("")
  g = g + theme_void()
  return(g)
}

ggsurvplot_empty = function(){
  g = ggsurvplot(survfit(Surv(time, status) ~ 1, data = lung),
                 title = "",
                 palette = "white",
                 ggtheme=theme_void())
  g$plot <- g$plot + theme(axis.text.x = element_blank(),
                           axis.text.y = element_blank(),
                           plot.title = element_blank(),
                           plot.subtitle = element_blank(),
                           legend.title=element_blank(),
                           legend.text=element_blank(),
                           axis.ticks=element_blank(),
                           axis.ticks.x = element_blank())
  return(g)
}

gg_drug_plot = function(data, x_name, y_name, x_lab, y_lab, Total_pts, x_max_patients){
  Cor = format_p(cor.test(data[,colnames(data) == x_name], data[,colnames(data) == y_name])$estimate,digits = 3)
  P = format_p(cor.test(data[,colnames(data) == x_name], data[,colnames(data) == y_name])$p.value,digits = 3)
  cor = data.frame(Cor, P)
  Max_x = max(data[,colnames(data) == x_name],na.rm = T) * 1.05
  Max_y = max(data[,colnames(data) == y_name],na.rm = T)
  g <- ggplot(data, aes_(x = as.name(x_name), y = as.name(y_name), color = as.name("Cancers")))
  g <- g + geom_point(size = Total_pts/x_max_patients*10, alpha = .5)
  g <- g + scale_fill_nejm()
  g <- g + coord_cartesian(xlim = c(0, Max_x), ylim = c(0, Max_y*1.05))
  g <- g + theme_bw()
  g <- g + labs(x = x_lab,y = y_lab, color = "Cancer type")
  g <- g + geom_smooth(method = lm, se = FALSE, na.rm = T, colour = "red")
  g <- g + geom_text(data=cor, colour = "red", mapping = aes(x = 0.3*Max_x, y = Max_y*1.03, label = paste("r =",Cor, ", p =",P)))
}

shannon.entropy = function(x, x_length){
  x = x[(!is.na(x))]
  x = x[(x != 0)]
  if(sum(x) <= 0) {
    return(log(x_length, base=2))
  } else {
    p = x/sum(x)
    score = -sum(p * log(p, base=2))
    return(score)
  }
}

odds.ratio <- function(a, b, c, d, correct=FALSE){
  cl <- function(x){
    or*exp(c(1,-1)*qnorm(x)*sqrt(1/a+1/b+1/c+1/d))
  }
  if (correct || a*b*c*d==0) {
    a <- a+0.5
    b <- b+0.5
    c <- c+0.5
    d <- d+0.5
  }
  or <- a*d/(b*c)
  conf <- rbind(cl90=cl(0.05), cl95=cl(0.025), cl99=cl(0.005), cl999=cl(0.0005))
  conf <- data.frame(conf)
  colnames(conf) <- paste(c("lower","upper"), " limit" , sep="")
  rownames(conf) <- paste(c(90, 95, 99, 99.9), "%CI" , sep="")
  list(or=or, conf=conf)
}

surv_curv_entry <- function(fit, data, title, legend_, diff_0, diff_1, diff_2=NULL){

  g = ggsurvplot(
    fit = fit,
    combine = TRUE,
    data = data,
    xlab = "Time from enrollment (months)",
    ylab = "Survival Probability",
    censor = TRUE,
    conf.int = FALSE,
    surv.scale = "percent",
    font.title = 8,
    font.subtitle = 8,
    font.main = 8,
    font.submain = 8,
    font.caption = 8,
    font.legend = 8,
    pval = FALSE,
    surv.median.line = "v",
    palette = "Dark2",
    risk.table = TRUE,
    risk.table.y.text = FALSE,
    cumevents = FALSE,
    cumcensor = FALSE,
    tables.theme = clean_theme(),
    legend = c(0.8,0.8),
    xlim = c(0, min(365.25*10, max(data$time_enroll_final, na.rm = TRUE)) * 1.05),
    xscale = "d_m",
    break.x.by = max(6 * 365.25 / 12, ceiling(max(data$time_enroll_final) / 365.25 / 6) * 365.25 / 2),
    legend.labs = legend_
  )

  # --- Median OS（群ごとに必ず作る） ---
  tab <- as.data.frame(summary(fit)$table)

  get_col <- function(df, cand) {
    nm <- names(df)
    hit <- cand[cand %in% nm]
    if (length(hit) >= 1) return(df[[hit[1]]])
    # 互換：X0.95LCLなど
    if ("0.95LCL" %in% cand) {
      h <- grep("0\\.95LCL", nm, fixed = FALSE, value = TRUE)
      if (length(h) >= 1) return(df[[h[1]]])
    }
    if ("0.95UCL" %in% cand) {
      h <- grep("0\\.95UCL", nm, fixed = FALSE, value = TRUE)
      if (length(h) >= 1) return(df[[h[1]]])
    }
    return(rep(NA_real_, nrow(df)))
  }

  med <- get_col(tab, c("median"))
  lcl <- get_col(tab, c("0.95LCL","X0.95LCL"))
  ucl <- get_col(tab, c("0.95UCL","X0.95UCL"))

  # medianが到達していないと NA になることがあるので、その場合の表示も整える
  fmt_median <- function(m, lo, up) {
    if (!is.finite(m)) return("NR")  # Not Reached
    paste0(
      format_p(digits = 1, m / 365.25 * 12),
      " (",
      ifelse(is.finite(lo), format_p(digits = 1, lo / 365.25 * 12), "NA"),
      "-",
      ifelse(is.finite(up), format_p(digits = 1, up / 365.25 * 12), "NA"),
      ")"
    )
  }

  legends <- vapply(seq_len(nrow(tab)), function(i) fmt_median(med[i], lcl[i], ucl[i]), character(1))
  legends_txt <- paste(legends, collapse = ", ")

  # --- HR 表示（従来通り） ---
  title_HR <- NULL
  if(!is.null(diff_2)){
    data_tmp <- tidy(diff_2, exponentiate=TRUE, conf.int=TRUE)
    if(nrow(data_tmp) > 0){
      for(i in 1:nrow(data_tmp)){
        title_HR <- c(title_HR,
                      paste0(format_p(data_tmp$estimate[i],digits=2), " (",
                             format_p(data_tmp$conf.low[i],digits=2), "-",
                             format_p(data_tmp$conf.high[i],digits=2), ") ",
                             "p=", format_p(data_tmp$p.value[i],digits=2)))
      }
    }
  }
  hr_txt <- if (!is.null(title_HR)) paste(title_HR, collapse = "/") else NULL

  # --- タイトル/サブタイトル ---
  if (is.null(diff_0)) {
    # p値無し（以前はここで1群扱いになってたのを修正）
    g$plot <- g$plot +
      labs(
        title = if (!is.null(hr_txt)) paste0(title, ": \nhazard ratio (vs 1st group)=", hr_txt) else title,
        subtitle = paste0("Median OS, ", legends_txt, " months")
      )
  } else {
    # p値あり
    p_lr <- 1 - pchisq(diff_0$chisq, length(diff_0$n)-1, lower.tail = TRUE)
    p_wx <- 1 - pchisq(diff_1$chisq, length(diff_1$n)-1, lower.tail = TRUE)

    g$plot <- g$plot +
      labs(
        title = paste0(title, ": log-rank, p=", format_p(p_lr, digits=3),
                       "/Wilcoxon, p=", format_p(p_wx, digits=3),
                       if (!is.null(hr_txt)) paste0("\nhazard ratio (vs 1st group)=", hr_txt) else ""),
        subtitle = paste0("Median OS, ", legends_txt, " months")
      )
  }

  g$table <- g$table + theme(plot.title = element_blank(),
                             plot.subtitle = element_blank())
  return(g)
}

weighted_survdiff_like <- function(time, status, group, weights, rho = 0) {
  # group: 2群想定（factor/character ok）
  # weights: IPTW等（非負、finite）
  ok <- is.finite(time) & is.finite(status) & !is.na(group) & is.finite(weights) & weights >= 0
  time <- time[ok]; status <- status[ok]; group <- group[ok]; weights <- weights[ok]

  g <- as.factor(group)
  if (nlevels(g) != 2) stop("weighted_survdiff_like() は2群のみ対応です。")

  g01 <- as.integer(g == levels(g)[2])  # 2nd level を1群扱い

  # イベント時刻（status==1）で計算
  event_times <- sort(unique(time[status == 1]))
  if (length(event_times) == 0) {
    return(list(chisq = NA_real_, n = as.numeric(table(g))))
  }

  # Wilcoxon(rho=1)等に必要：pooled weighted KMで S(t-) を作る
  # rho=0なら全部1になるので計算は軽い
  if (rho != 0) {
    sf_all <- survival::survfit(survival::Surv(time, status) ~ 1, weights = weights)
    # sf_all$time はイベント/検閲の変化点、sf_all$surv は直後の生存
    # S(t-) を得るため、各 event_time に対し「直前のsurv」を拾う
    # (t=0は1)
    tt <- c(0, sf_all$time)
    ss <- c(1, sf_all$surv)
    S_tminus <- function(t) {
      idx <- max(which(tt < t))
      ss[idx]
    }
    w_rho <- vapply(event_times, function(t) S_tminus(t)^rho, numeric(1))
  } else {
    w_rho <- rep(1, length(event_times))
  }

  U <- 0
  V <- 0

  for (k in seq_along(event_times)) {
    t0 <- event_times[k]
    wk <- w_rho[k]

    at_risk <- time >= t0
    # weighted risk sets
    Y1 <- sum(weights[at_risk & (g01 == 1)])
    Y0 <- sum(weights[at_risk & (g01 == 0)])
    Y  <- Y0 + Y1
    if (!is.finite(Y) || Y <= 0) next

    # weighted events at t0
    at_event <- (time == t0) & (status == 1)
    d1 <- sum(weights[at_event & (g01 == 1)])
    d0 <- sum(weights[at_event & (g01 == 0)])
    d  <- d0 + d1
    if (!is.finite(d) || d <= 0) next

    exp_d1 <- d * (Y1 / Y)

    U <- U + wk * (d1 - exp_d1)

    # 分散（連続重みでの近似）：wk^2 * d * p * (1-p)
    p1 <- Y1 / Y
    V <- V + (wk^2) * d * p1 * (1 - p1)
  }

  chisq <- if (is.finite(V) && V > 0) (U^2 / V) else NA_real_
  list(chisq = chisq, n = as.numeric(table(g)))
}

surv_curv_CTx <- function(fit, data, title, legend_, diff_0, Xlab = "Time from CTx start, risk-set adjusted (months)"){
  g = ggsurvplot(
    fit = fit,
    combine = TRUE,
    data = data,
    xlab = Xlab,
    ylab = "Survival Probability",
    censor = TRUE,
    conf.int = FALSE,
    surv.scale = "percent",
    font.title = 8,
    font.subtitle = 8,
    font.main = 8,
    font.submain = 8,
    font.caption = 8,
    font.legend = 8,
    pval = FALSE,
    surv.median.line = "v",
    palette = "Dark2",
    risk.table = TRUE,
    risk.table.y.text = FALSE,
    cumevents = FALSE,
    cumcensor = FALSE,
    tables.theme = clean_theme(),
    legend = c(0.8,0.8),
    xlim = c(0, min(365.25*10, max(data$time_all, na.rm = T)) * 1.05),
    xscale = "d_m",
    # break.x.by = 12 * 365.25 / 12,
    break.x.by = max(6 * 365.25 / 12, ceiling(max(data$time_all) / 365.25 / 6) * 365.25 / 2),
    legend.labs = legend_
  )
  if(is.null(diff_0)){
    tmp = summary(fit)$table
    legends = paste0(
      format_p(digits = 1, tmp[[7]] / 365.25 * 12),
      " (",
      format_p(digits = 1, tmp[[8]] / 365.25 * 12),
      "-",
      format_p(digits = 1, tmp[[9]] / 365.25 * 12),
      ")"
    )
    g$plot <- g$plot +
      labs(title = title,
           subtitle = paste0("Median OS, ", legends, " months"))
  } else{
    tmp = data.frame(summary(fit)$table)
    legends = paste0(
      format_p(digits = 1, tmp$median[1] / 365.25 * 12),
      " (",
      format_p(digits = 1, tmp$X0.95LCL[1] / 365.25 * 12),
      "-",
      format_p(digits = 1, tmp$X0.95UCL[1] / 365.25 * 12),
      ")"
    )
    for(i in 2:nrow(tmp)){
      legends = paste(
        legends,
        paste0(
          format_p(digits = 1, tmp$median[i] / 365.25 * 12),
          " (",
          format_p(digits = 1, tmp$X0.95LCL[i] / 365.25 * 12),
          "-",
          format_p(digits = 1, tmp$X0.95UCL[i] / 365.25 * 12),
          ")"
        ),
        sep = ", "
      )
    }
    data_tmp = tidy(diff_0, exponentiate=TRUE, conf.int=TRUE)
    title_HR = NULL
    if(nrow(data_tmp) > 0){
      for(i in 1:nrow(data_tmp)){
        title_HR = c(title_HR,
                     paste0(format_p(data_tmp$estimate[i],digits=2), " (",
                            format_p(data_tmp$conf.low[i],digits=2), "-",
                            format_p(data_tmp$conf.high[i],digits=2), ") ",
                            "p=", format_p(data_tmp$p.value[i],digits=2)))
      }
    }
    g$plot <- g$plot +
      labs(title = paste0(title, ": \nhazard ratio (vs 1st group)=", paste(title_HR, collapse = "/")),
           subtitle = paste0("Median OS, ", legends, " months"))
  }
  g$table <- g$table + theme(plot.title = element_blank(),
                             plot.subtitle = element_blank())
  return(g)
}

surv_curv_drug <- function(fit, data, title, diff_0, diff_1, diff_2=NULL){
  g = ggsurvplot(
    fit = fit,
    combine = TRUE,
    data = data,
    xlab = "Time from drug initiation (months)",
    ylab = "Survival Probability",
    censor = TRUE,
    surv.scale = "percent",
    conf.int = FALSE,
    font.title = 8,
    font.subtitle = 8,
    font.main = 8,
    font.submain = 8,
    font.caption = 8,
    font.legend = 8,
    pval = FALSE,
    surv.median.line = "v",
    palette = "Dark2",
    risk.table = TRUE,
    risk.table.y.text = FALSE,
    cumevents = FALSE,
    cumcensor = FALSE,
    tables.theme = clean_theme(),
    legend = c(0.8,0.8),
    xscale = "d_m",
    xlim = c(0, 3*365.35),
    break.x.by = 3 * 365.25 / 12,
    # break.x.by = ceiling(max(data$time_enroll_final) / 365.25 / 6,
  )
  if(is.null(diff_0)){
    if(is.null(diff_2)){
      tmp = summary(fit)$table
      legends = paste0(
        format_p(digits = 1, tmp[[7]] / 365.25 * 12),
        " (",
        format_p(digits = 1, tmp[[8]] / 365.25 * 12),
        "-",
        format_p(digits = 1, tmp[[9]] / 365.25 * 12),
        ")"
      )
      g$plot <- g$plot +
        labs(title = title,
             subtitle = paste0("Median OS, ", legends, " months"))
      g$table <- g$table + theme(plot.title = element_blank(),
                                 plot.subtitle = element_blank())
    } else {
      tmp = summary(fit)$table
      legends = paste0(
        format_p(digits = 1, tmp[[7]] / 365.25 * 12),
        " (",
        format_p(digits = 1, tmp[[8]] / 365.25 * 12),
        "-",
        format_p(digits = 1, tmp[[9]] / 365.25 * 12),
        ")"
      )
      data_tmp = tidy(diff_2, exponentiate=TRUE, conf.int=TRUE)
      title_HR = NULL
      if(nrow(data_tmp) > 0){
        for(i in 1:nrow(data_tmp)){
          title_HR = c(title_HR,
                       paste0(format_p(data_tmp$estimate[i],digits=2), " (",
                              format_p(data_tmp$conf.low[i],digits=2), "-",
                              format_p(data_tmp$conf.high[i],digits=2), ") ",
                              "p=", format_p(data_tmp$p.value[i],digits=3)))
        }
      }
      g$plot <- g$plot +
        labs(title = paste0(title, ": hazard ratio=", paste(title_HR, collapse = "/")),
             subtitle = paste0("Median OS, ", legends, " months"))
      g$table <- g$table + theme(plot.title = element_blank(),
                                 plot.subtitle = element_blank())
    }
  } else{
    tmp = data.frame(summary(fit)$table)
    legends = paste0(
      format_p(digits = 1, tmp$median[1] / 365.25 * 12),
      " (",
      format_p(digits = 1, tmp$X0.95LCL[1] / 365.25 * 12),
      "-",
      format_p(digits = 1, tmp$X0.95UCL[1] / 365.25 * 12),
      ")"
    )
    for(i in 2:nrow(tmp)){
      legends = paste(
        legends,
        paste0(
          format_p(digits = 1, tmp$median[i] / 365.25 * 12),
          " (",
          format_p(digits = 1, tmp$X0.95LCL[i] / 365.25 * 12),
          "-",
          format_p(digits = 1, tmp$X0.95UCL[i] / 365.25 * 12),
          ")"
        ),
        sep = ", "
      )
    }
    data_tmp = tidy(diff_2, exponentiate=TRUE, conf.int=TRUE)
    title_HR = NULL
    if(nrow(data_tmp) > 0){
      for(i in 1:nrow(data_tmp)){
        title_HR = c(title_HR,
                     paste0(format_p(data_tmp$estimate[i],digits=2), " (",
                            format_p(data_tmp$conf.low[i],digits=2), "-",
                            format_p(data_tmp$conf.high[i],digits=2), ") ",
                            "p=", format_p(data_tmp$p.value[i],digits=3)))
      }
    }
    # g$plot <- g$plot +
    #   labs(title = paste0(title, ": hazard ratio=", paste(title_HR, collapse = "/")),
    #        subtitle = paste0("Median OS, ", legends, " months"))
    g$plot <- g$plot +
      labs(title = paste0(title, ": log-rank, p=", format_p(
        1 - pchisq(diff_0$chisq, length(diff_0$n)-1, lower.tail = TRUE),digits=3),
        "/Wilcoxon, p=", format_p(
          1 - pchisq(diff_1$chisq, length(diff_1$n)-1, lower.tail = TRUE),digits=3)),
        subtitle = paste0("Median OS, ", legends," months, hazard ratio (vs 1st group)=", paste(title_HR, collapse = "/")))
    g$table <- g$table + theme(plot.title = element_blank(),
                               plot.subtitle = element_blank())
  }
  return(g)
}

surv_curv_drug_IPTW <- function(fit, data, title, diff_2=NULL){
  g = ggsurvplot(
    fit = fit,
    combine = TRUE,
    data = data,
    xlab = "Time from drug initiation (months)",
    ylab = "Survival Probability, IPTW corrected",
    censor = TRUE,
    surv.scale = "percent",
    conf.int = FALSE,
    font.title = 8,
    font.subtitle = 8,
    font.main = 8,
    font.submain = 8,
    font.caption = 8,
    font.legend = 8,
    pval = FALSE,
    surv.median.line = "v",
    palette = "Dark2",
    risk.table = TRUE,
    risk.table.y.text = FALSE,
    cumevents = FALSE,
    cumcensor = FALSE,
    tables.theme = clean_theme(),
    legend = c(0.8,0.8),
    xscale = "d_m",
    xlim = c(0, 3*365.35),
    break.x.by = 3 * 365.25 / 12,
    # break.x.by = ceiling(max(data$time_enroll_final) / 365.25 / 6,
  )
  tmp = summary(fit)$table
  legends = paste0(
    format_p(digits = 1, tmp[[7]] / 365.25 * 12),
    " (",
    format_p(digits = 1, tmp[[8]] / 365.25 * 12),
    "-",
    format_p(digits = 1, tmp[[9]] / 365.25 * 12),
    ")"
  )
  tmp = data.frame(summary(fit)$table)
  legends = paste0(
    format_p(digits = 1, tmp$median[1] / 365.25 * 12),
    " (",
    format_p(digits = 1, tmp$X0.95LCL[1] / 365.25 * 12),
    "-",
    format_p(digits = 1, tmp$X0.95UCL[1] / 365.25 * 12),
    ")"
  )
  for(i in 2:nrow(tmp)){
    legends = paste(
      legends,
      paste0(
        format_p(digits = 1, tmp$median[i] / 365.25 * 12),
        " (",
        format_p(digits = 1, tmp$X0.95LCL[i] / 365.25 * 12),
        "-",
        format_p(digits = 1, tmp$X0.95UCL[i] / 365.25 * 12),
        ")"
      ),
      sep = ", "
    )
  }
  data_tmp = tidy(diff_2, exponentiate=TRUE, conf.int=TRUE)
  title_HR = NULL
  if(nrow(data_tmp) > 0){
    for(i in 1:nrow(data_tmp)){
      title_HR = c(title_HR,
                   paste0(format_p(data_tmp$estimate[i],digits=2), " (",
                          format_p(data_tmp$conf.low[i],digits=2), "-",
                          format_p(data_tmp$conf.high[i],digits=2), ") ",
                          "p=", format_p(data_tmp$p.value[i],digits=3)))
    }
  }
  g$plot <- g$plot +
    labs(title = paste0(title, ": hazard ratio=", paste(title_HR, collapse = "/")),
         subtitle = paste0("Median OS, ", legends, " months. Risk table is not valid due to IPTW correction."))
  g$table <- g$table + theme(plot.title = element_blank(),
                             plot.subtitle = element_blank())
  return(g)
}

survival_compare_and_plot <- function(data,
                                      time_var = "time_enroll_final",
                                      status_var = "censor",
                                      group_var = "treatment",
                                      plot_title = "Survival curve",
                                      input_rmst_cgp = 2,
                                      group_labels = NULL,
                                      weight_var = NULL) {
  surv_formula <- as.formula(paste0("Surv(", time_var, ", ", status_var, ") ~ ", group_var))
  surv_fit <- eval(substitute(
    survfit(formula = FORMULA, data = data, conf.type = "log-log"),
    list(FORMULA = surv_formula)
  ))
  if(group_var == "1"){
    surv_curv_entry(surv_fit, data, plot_title, NULL, NULL, NULL)
  } else{
    group_vals <- unique(na.omit(data[[group_var]]))
    n_group <- length(group_vals)
    if (n_group == 1) {
      surv_curv_entry(surv_fit, data, plot_title, NULL, NULL, NULL)
    } else if (n_group == 2) {
      surv_obj <- with(data, Surv(get(time_var), get(status_var)))
      # tau: min of max times in each group or input
      data <- data %>%
        mutate(treatment_numeric = ifelse(.data[[group_var]] == group_vals[1], 0,
                                          ifelse(.data[[group_var]] == group_vals[2], 1, NA)))

      tau0 <- max(data[data[[group_var]] == group_vals[1], ] %>% pull(!!sym(time_var)), na.rm = TRUE) / 365.25
      tau1 <- max(data[data[[group_var]] == group_vals[2], ] %>% pull(!!sym(time_var)), na.rm = TRUE) / 365.25
      tau <- floor(min(tau0, tau1, input_rmst_cgp) * 10) / 10

      # rmst
      rmst_result <- rmst2(
        time = data[[time_var]],
        status = data[[status_var]],
        arm = data$treatment_numeric,
        tau = tau * 365.25
      )

      # 差と信頼区間
      rmst_diff <- format_p(rmst_result$unadjusted.result[1], digits = 1)
      rmst_ll <- format_p(rmst_result$unadjusted.result[4], digits = 1)
      rmst_ul <- format_p(rmst_result$unadjusted.result[7], digits = 1)

      # log-rank
      diff_0 <- survdiff(surv_obj ~ data[[group_var]], rho = 0)
      diff_1 <- survdiff(surv_obj ~ data[[group_var]], rho = 1)
      diff_2 <- coxph(surv_obj ~ data[[group_var]])
      # title更新
      plot_title_full <- paste0(plot_title, ", ", tau, "-year RMST diff.: ", rmst_diff,
                                " (", rmst_ll, "-", rmst_ul, ") days")

      surv_curv_entry(surv_fit, data, plot_title_full, group_labels, diff_0, diff_1, diff_2)

    } else if (n_group > 2) {
      surv_obj <- with(data, Surv(get(time_var), get(status_var)))
      # 通常の log-rankのみ
      diff_0 <- survdiff(surv_obj ~ data[[group_var]], rho = 0)
      diff_1 <- survdiff(surv_obj ~ data[[group_var]], rho = 1)
      diff_2 <- coxph(surv_obj ~ data[[group_var]])

      surv_curv_entry(surv_fit, data, plot_title, group_labels, diff_0, diff_1, diff_2)
    }
  }
}

survival_compare_and_plot_match <- function(data,
                                            time_var = "time_enroll_final",
                                            status_var = "censor",
                                            group_var = "treatment",
                                            plot_title = "Survival curve",
                                            input_rmst_cgp = 2,
                                            group_labels = NULL,
                                            pair_var = NULL,
                                            n_boot = 1000,
                                            seed = 1,
                                            weight_var = NULL) {

  # --- weights 取得（無ければNULL） ---
  w <- NULL
  if (!is.null(weight_var) && weight_var %in% names(data)) {
    w <- data[[weight_var]]
    # 念のため
    w[!is.finite(w)] <- NA
  }

  surv_formula <- as.formula(paste0("Surv(", time_var, ", ", status_var, ") ~ ", group_var))

  # ★ survfit に weights を渡す（IPTWカーブ用）
  if (is.null(w)) {
    surv_fit <- eval(substitute(
      survfit(formula = FORMULA, data = data, conf.type = "log-log"),
      list(FORMULA = surv_formula)
    ))
  } else {
    surv_fit <- eval(substitute(
      survfit(formula = FORMULA, data = data, weights = W, conf.type = "log-log"),
      list(FORMULA = surv_formula, W = w)
    ))
  }

  if (group_var == "1") {
    return(surv_curv_entry(surv_fit, data, plot_title, NULL, NULL, NULL, NULL))
  }

  group_vals <- unique(na.omit(data[[group_var]]))
  n_group <- length(group_vals)

  if (n_group == 1) {
    return(surv_curv_entry(surv_fit, data, plot_title, NULL, NULL, NULL, NULL))
  }

  # surv_obj は共通
  surv_obj <- with(data, Surv(get(time_var), get(status_var)))

  if (n_group == 2) {

    data <- data %>%
      mutate(treatment_numeric = ifelse(.data[[group_var]] == group_vals[1], 0,
                                        ifelse(.data[[group_var]] == group_vals[2], 1, NA)))

    tau0 <- max(data[data[[group_var]] == group_vals[1], ] %>% pull(!!sym(time_var)), na.rm = TRUE) / 365.25
    tau1 <- max(data[data[[group_var]] == group_vals[2], ] %>% pull(!!sym(time_var)), na.rm = TRUE) / 365.25
    tau  <- floor(min(tau0, tau1, input_rmst_cgp) * 10) / 10

    # --- RMST ---
    if (is.null(w)) {
      # 従来通り（rmst2）
      rmst_result <- rmst2(
        time   = data[[time_var]],
        status = data[[status_var]],
        arm    = data$treatment_numeric,
        tau    = tau * 365.25
      )
      rmst_diff_point <- rmst_result$unadjusted.result[1]

      # CI（pair bootstrap or rmst2 CI）
      if (!is.null(pair_var) && pair_var %in% names(data)) {
        set.seed(seed)
        pairs <- unique(data[[pair_var]])

        boot_rmst_diff <- function(d) {
          if (length(unique(d$treatment_numeric)) < 2) return(NA_real_)
          out <- rmst2(
            time   = d[[time_var]],
            status = d[[status_var]],
            arm    = d$treatment_numeric,
            tau    = tau * 365.25
          )
          out$unadjusted.result[1]
        }

        boot_vals <- replicate(n_boot, {
          sampled_pairs <- sample(pairs, size = length(pairs), replace = TRUE)
          d_boot <- data[data[[pair_var]] %in% sampled_pairs, , drop = FALSE]
          boot_rmst_diff(d_boot)
        })

        boot_vals <- boot_vals[is.finite(boot_vals)]
        rmst_ll <- unname(quantile(boot_vals, 0.025, na.rm = TRUE))
        rmst_ul <- unname(quantile(boot_vals, 0.975, na.rm = TRUE))

      } else {
        rmst_ll <- rmst_result$unadjusted.result[4]
        rmst_ul <- rmst_result$unadjusted.result[7]
      }

    } else {
      # ★ IPTW時：重み付きKM（step関数）を tau_days まで厳密積分してRMST差(days)を推定し、
      #           bootstrap CI も days スケールで整合させる

      tau_years <- floor(min(tau0, tau1, input_rmst_cgp) * 10) / 10
      tau_days  <- tau_years * 365.25

      rmst_ci <- rmst_diff_iptw_boot_ci(
        data = data,
        time_var = time_var,
        status_var = status_var,
        group_var = group_var,
        weight_var = weight_var,   # ★列名を渡す
        tau_days = tau_days,
        pair_var = pair_var,
        n_boot = n_boot,
        seed = seed
      )

      rmst_diff_point <- rmst_ci$point
      rmst_ll <- rmst_ci$ll
      rmst_ul <- rmst_ci$ul
    }
    rmst_diff <- format_p(rmst_diff_point, digits = 1)
    rmst_ll_f <- format_p(rmst_ll, digits = 1)
    rmst_ul_f <- format_p(rmst_ul, digits = 1)

    # --- log-rank / Wilcoxon ---
    if (is.null(w)) {
      diff_0 <- survdiff(surv_obj ~ data[[group_var]], rho = 0)
      diff_1 <- survdiff(surv_obj ~ data[[group_var]], rho = 1)
    } else {
      # ★重み付き log-rank / Wilcoxon 相当
      wd0 <- weighted_survdiff_like(
        time = data[[time_var]],
        status = data[[status_var]],
        group = data[[group_var]],
        weights = w,
        rho = 0
      )
      wd1 <- weighted_survdiff_like(
        time = data[[time_var]],
        status = data[[status_var]],
        group = data[[group_var]],
        weights = w,
        rho = 1
      )
      # survdiff と同じ形っぽくして既存コードを活かす
      diff_0 <- list(chisq = wd0$chisq, n = wd0$n)
      diff_1 <- list(chisq = wd1$chisq, n = wd1$n)
    }

    # --- Cox ---
    if (!is.null(pair_var) && pair_var %in% names(data)) {
      cox_formula <- as.formula(
        paste0("Surv(", time_var, ", ", status_var, ") ~ ",
               group_var, " + strata(", pair_var, ")")
      )
    } else {
      cox_formula <- as.formula(
        paste0("Surv(", time_var, ", ", status_var, ") ~ ", group_var)
      )
    }

    if (is.null(w)) {
      diff_2 <- coxph(cox_formula, data = data)
    } else {
      diff_2 <- coxph(cox_formula, data = data, weights = w, robust = TRUE)
    }

    plot_title_full <- paste0(
      plot_title, ", ", tau, "-year RMST diff.: ",
      rmst_diff, " (", rmst_ll_f, "-", rmst_ul_f, ") days"
    )

    return(surv_curv_entry(surv_fit, data, plot_title_full,
                           group_labels, diff_0, diff_1, diff_2))

  } else {  # n_group > 2

    if (is.null(w)) {
      diff_0 <- survdiff(surv_obj ~ data[[group_var]], rho = 0)
      diff_1 <- survdiff(surv_obj ~ data[[group_var]], rho = 1)
    } else {
      diff_0 <- NULL
      diff_1 <- NULL
    }

    if (!is.null(pair_var) && pair_var %in% names(data)) {
      cox_formula <- as.formula(
        paste0("Surv(", time_var, ", ", status_var, ") ~ ",
               group_var, " + strata(", pair_var, ")")
      )
    } else {
      cox_formula <- as.formula(
        paste0("Surv(", time_var, ", ", status_var, ") ~ ", group_var)
      )
    }

    if (is.null(w)) {
      diff_2 <- coxph(cox_formula, data = data)
    } else {
      diff_2 <- coxph(cox_formula, data = data, weights = w, robust = TRUE)
    }

    return(surv_curv_entry(surv_fit, data, plot_title,
                           group_labels, diff_0, diff_1, diff_2))
  }
}
# survival_compare_and_plot_match <- function(data,
#                                             time_var = "time_enroll_final",
#                                             status_var = "censor",
#                                             group_var = "treatment",
#                                             plot_title = "Survival curve",
#                                             input_rmst_cgp = 2,
#                                             group_labels = NULL,
#                                             pair_var = NULL,     # ペアID列名（例 "pair_id" or "subclass"）
#                                             n_boot = 1000,
#                                             seed = 1) {
#
#   surv_formula <- as.formula(paste0("Surv(", time_var, ", ", status_var, ") ~ ", group_var))
#   surv_fit <- eval(substitute(
#     survfit(formula = FORMULA, data = data, conf.type = "log-log"),
#     list(FORMULA = surv_formula)
#   ))
#
#   if (group_var == "1") {
#     return(surv_curv_entry(surv_fit, data, plot_title, NULL, NULL, NULL))
#   }
#
#   group_vals <- unique(na.omit(data[[group_var]]))
#   n_group <- length(group_vals)
#
#   if (n_group == 1) {
#     return(surv_curv_entry(surv_fit, data, plot_title, NULL, NULL, NULL))
#   }
#
#   # surv_obj は共通で使う
#   surv_obj <- with(data, Surv(get(time_var), get(status_var)))
#
#   if (n_group == 2) {
#
#     data <- data %>%
#       mutate(treatment_numeric = ifelse(.data[[group_var]] == group_vals[1], 0,
#                                         ifelse(.data[[group_var]] == group_vals[2], 1, NA)))
#
#     tau0 <- max(data[data[[group_var]] == group_vals[1], ] %>% pull(!!sym(time_var)), na.rm = TRUE) / 365.25
#     tau1 <- max(data[data[[group_var]] == group_vals[2], ] %>% pull(!!sym(time_var)), na.rm = TRUE) / 365.25
#     tau  <- floor(min(tau0, tau1, input_rmst_cgp) * 10) / 10
#
#     # RMST point estimate
#     rmst_result <- rmst2(
#       time   = data[[time_var]],
#       status = data[[status_var]],
#       arm    = data$treatment_numeric,
#       tau    = tau * 365.25
#     )
#     rmst_diff_point <- rmst_result$unadjusted.result[1]
#
#     # CI: pair bootstrap if pair_var is available
#     if (!is.null(pair_var) && pair_var %in% names(data)) {
#
#       set.seed(seed)
#       pairs <- unique(data[[pair_var]])
#
#       boot_rmst_diff <- function(d){
#         if (length(unique(d$treatment_numeric)) < 2) return(NA_real_)
#         out <- rmst2(
#           time   = d[[time_var]],
#           status = d[[status_var]],
#           arm    = d$treatment_numeric,
#           tau    = tau * 365.25
#         )
#         out$unadjusted.result[1]
#       }
#
#       boot_vals <- replicate(n_boot, {
#         sampled_pairs <- sample(pairs, size = length(pairs), replace = TRUE)
#         d_boot <- data[data[[pair_var]] %in% sampled_pairs, , drop = FALSE]
#         boot_rmst_diff(d_boot)
#       })
#
#       boot_vals <- boot_vals[is.finite(boot_vals)]
#       rmst_ll <- unname(quantile(boot_vals, 0.025, na.rm=TRUE))
#       rmst_ul <- unname(quantile(boot_vals, 0.975, na.rm=TRUE))
#
#     } else {
#       rmst_ll <- rmst_result$unadjusted.result[4]
#       rmst_ul <- rmst_result$unadjusted.result[7]
#     }
#
#     rmst_diff <- format_p(rmst_diff_point, digits = 1)
#     rmst_ll_f <- format_p(rmst_ll, digits = 1)
#     rmst_ul_f <- format_p(rmst_ul, digits = 1)
#
#     # log-rank / Wilcoxon
#     diff_0 <- survdiff(surv_obj ~ data[[group_var]], rho = 0)
#     diff_1 <- survdiff(surv_obj ~ data[[group_var]], rho = 1)
#
#     # ★ Cox: strata(pair_var) を条件付きで入れる
#     if (!is.null(pair_var) && pair_var %in% names(data)) {
#       cox_formula <- as.formula(
#         paste0("Surv(", time_var, ", ", status_var, ") ~ ",
#                group_var, " + strata(", pair_var, ")")
#       )
#       diff_2 <- coxph(cox_formula, data = data)
#     } else {
#       cox_formula <- as.formula(
#         paste0("Surv(", time_var, ", ", status_var, ") ~ ", group_var)
#       )
#       diff_2 <- coxph(cox_formula, data = data)
#     }
#
#     plot_title_full <- paste0(
#       plot_title, ", ", tau, "-year RMST diff.: ",
#       rmst_diff, " (", rmst_ll_f, "-", rmst_ul_f, ") days"
#     )
#
#     return(surv_curv_entry(surv_fit, data, plot_title_full,
#                            group_labels, diff_0, diff_1, diff_2))
#
#   } else {  # n_group > 2
#
#     diff_0 <- survdiff(surv_obj ~ data[[group_var]], rho = 0)
#     diff_1 <- survdiff(surv_obj ~ data[[group_var]], rho = 1)
#
#     # ★ Cox: strata(pair_var) 条件付き
#     if (!is.null(pair_var) && pair_var %in% names(data)) {
#       cox_formula <- as.formula(
#         paste0("Surv(", time_var, ", ", status_var, ") ~ ",
#                group_var, " + strata(", pair_var, ")")
#       )
#       diff_2 <- coxph(cox_formula, data = data)
#     } else {
#       cox_formula <- as.formula(
#         paste0("Surv(", time_var, ", ", status_var, ") ~ ", group_var)
#       )
#       diff_2 <- coxph(cox_formula, data = data)
#     }
#
#     return(surv_curv_entry(surv_fit, data, plot_title,
#                            group_labels, diff_0, diff_1, diff_2))
#   }
# }

# --- KM step関数を tau_days まで「厳密に」積分して RMST(days) を返す ---
rmst_from_survfit_step_exact <- function(sf, tau_days) {
  if (!is.finite(tau_days) || tau_days <= 0) return(NA_real_)

  # survfit の time/surv は t=0 を含まないので補う
  tt <- c(0, sf$time)
  ss <- c(1, sf$surv)

  # 変化点が無ければ S(t)=1 として tau まで
  if (length(tt) == 1) return(tau_days)

  # tau が最初の変化点より小さいなら S(t)=1 のまま
  if (tau_days <= tt[2]) return(tau_days)

  # tau 未満の変化点を取り出し、最後に tau を追加（Sは直前値で延長）
  keep <- tt < tau_days
  tt2 <- tt[keep]
  ss2 <- ss[keep]
  tt2 <- c(tt2, tau_days)
  ss2 <- c(ss2, tail(ss2, 1))

  # step関数の厳密積分：区間 [t_i, t_{i+1}) の値は ss2[i]
  sum(diff(tt2) * head(ss2, -1))
}

# --- 2群のIPTW RMST差（treated - control）を days で返す ---
rmst_diff_iptw_exact <- function(data, time_var, status_var, group_var, weight_var, tau_days) {
  # 前処理
  g <- as.factor(data[[group_var]])
  if (nlevels(g) != 2) return(NA_real_)

  w <- data[[weight_var]]
  ok <- is.finite(data[[time_var]]) & is.finite(data[[status_var]]) & !is.na(g) &
    is.finite(w) & w >= 0
  d <- data[ok, , drop = FALSE]
  g <- droplevels(as.factor(d[[group_var]]))
  if (nlevels(g) != 2) return(NA_real_)
  w <- d[[weight_var]]

  lev <- levels(g)
  d0 <- d[g == lev[1], , drop = FALSE]
  d1 <- d[g == lev[2], , drop = FALSE]
  w0 <- w[g == lev[1]]
  w1 <- w[g == lev[2]]

  # group別に重み付きKM
  sf0 <- survival::survfit(survival::Surv(d0[[time_var]], d0[[status_var]]) ~ 1,
                           data = d0, weights = w0)
  sf1 <- survival::survfit(survival::Surv(d1[[time_var]], d1[[status_var]]) ~ 1,
                           data = d1, weights = w1)

  rmst0 <- rmst_from_survfit_step_exact(sf0, tau_days)
  rmst1 <- rmst_from_survfit_step_exact(sf1, tau_days)

  rmst1 - rmst0  # days
}

# --- IPTW RMST差のCI（bootstrap, percentile）---
# pair_var があればペア単位、無ければ weights に比例したPPS bootstrap
rmst_diff_iptw_boot_ci <- function(data, time_var, status_var, group_var,
                                   weight_var, tau_days,
                                   pair_var = NULL,
                                   n_boot = 2000, seed = 1) {
  set.seed(seed)

  # 点推定（ここは元のまま正解）
  point <- rmst_diff_iptw_exact(data, time_var, status_var, group_var, weight_var, tau_days)

  boot_once <- function() {
    if (!is.null(pair_var) && pair_var %in% names(data)) {
      # ペアがある場合（PSM）: ペア単位でリサンプリング
      pairs <- unique(data[[pair_var]])
      sp <- sample(pairs, length(pairs), replace = TRUE)
      # ここはOK
      idx <- unlist(lapply(sp, function(p) which(data[[pair_var]] == p)))
      d_boot <- data[idx, , drop = FALSE]
    } else {
      # --- 修正箇所: IPWの場合 ---
      w <- data[[weight_var]]
      ok <- is.finite(w) & w >= 0
      if (!any(ok)) return(NA_real_)

      idx_ok <- which(ok)

      # 【修正】prob = p を削除。単純な復元抽出にする。
      # 重み付き解析なので、データは公平に選び、解析関数内で weights=w を効かせるのが正しい。
      idx <- sample(idx_ok, size = length(idx_ok), replace = TRUE)

      d_boot <- data[idx, , drop = FALSE]
    }

    # ここで weights=w が使われるので、抽出自体はランダムで良い
    rmst_diff_iptw_exact(d_boot, time_var, status_var, group_var, weight_var, tau_days)
  }

  boots <- replicate(n_boot, boot_once())
  boots <- boots[is.finite(boots)]

  # 点推定値とCIのズレを確認するため、mean/medianも見ておくと良い
  # print(summary(boots))

  ll <- unname(stats::quantile(boots, 0.025, na.rm = TRUE))
  ul <- unname(stats::quantile(boots, 0.975, na.rm = TRUE))

  list(point = point, ll = ll, ul = ul, boots = boots)
}

survival_compare_and_plot_CTx <- function(data,
                                          time_var1 = "time_pre",
                                          time_var2 = "time_all",
                                          status_var = "censor",
                                          group_var = "treatment",
                                          plot_title = "Survival curve",
                                          adjustment = TRUE,
                                          color_var_surv_CTx_1 = "CTx",
                                          group_labels = NULL,
                                          weights_var = NULL) {

  data$time_pre = data[[time_var1]]
  data$time_all = data[[time_var2]]

  # Check if valid weights are provided
  use_weights <- !is.null(weights_var) && (weights_var %in% colnames(data))

  if (use_weights) {
    # =========================================================================
    # FIX: Use Left Truncation (time_pre, time_all) combined with IPTW weights.
    # This prevents the Immortal Time Bias (which caused the 5-year median OS),
    # while the IPTW weights prevent the hazard explosion of the standard model.
    # =========================================================================
    surv_formula <- as.formula(paste0("Surv(", time_var1, ", ", time_var2, ", ", status_var, ") ~ ", group_var))
    Xlab <- paste0("Time from ", color_var_surv_CTx_1, ", IPTW adjusted (months)")
    weight_vector <- data[[weights_var]]
  } else {
    # Traditional behavior (supports older function calls)
    if (adjustment) {
      surv_formula <- as.formula(paste0("Surv(", time_var1, ", ", time_var2, ", ", status_var, ") ~ ", group_var))
      Xlab <- paste0("Time from ", color_var_surv_CTx_1, ", risk-set adjusted (months)")
    } else if (plot_title == "Unbiased OS Simulation (IPTW + Left-Truncated Llogis + Rank Match)") {
      surv_formula <- as.formula(paste0("Surv(", time_var2, ", ", status_var, ") ~ ", group_var))
      Xlab <- paste0("Time from ", color_var_surv_CTx_1, ", simulation (months)")
    } else {
      surv_formula <- as.formula(paste0("Surv(", time_var2, ", ", status_var, ") ~ ", group_var))
      Xlab <- paste0("Time from ", color_var_surv_CTx_1, ", bias not adj (months)")
    }
    weight_vector <- rep(1, nrow(data))
  }

  # Fit model
  surv_fit <- eval(substitute(
    survfit(formula = FORMULA, data = data, weights = WEIGHTS, conf.type = "log-log"),
    list(FORMULA = surv_formula, WEIGHTS = weight_vector)
  ))

  if(group_var == "1"){
    surv_curv_CTx(surv_fit, data, plot_title, NULL, NULL, Xlab)
  } else {
    group_vals <- unique(na.omit(data[[group_var]]))
    n_group <- length(group_vals)
    if (n_group == 1) {
      surv_curv_CTx(surv_fit, data, plot_title, NULL, NULL, Xlab)
    } else {
      # Use robust=TRUE when weights are applied
      if (use_weights) {
        diff_0 <- coxph(surv_formula, data = data, weights = weight_vector, robust = TRUE)
      } else {
        diff_0 <- coxph(surv_formula, data = data)
      }
      surv_curv_CTx(surv_fit, data, plot_title, group_labels, diff_0, Xlab)
    }
  }
}

rename_factors_survival_CGP <- function(names_vec) {
  name_map <- c(
    "Lymph_met" = "Lymphatic metastasis",
    "Liver_met" = "Liver metastasis",
    "Lung_met" = "Lung metastasis",
    "Bone_met" = "Bone metastasis",
    "Brain_met" = "Brain metastasis",
    "EP_option" = "Treatment recommended",
    "PS" = "Performance status",
    "Lines" = "CTx lines before CGP",
    "Best_effect" = "Best CTx effect before CGP"
  )
  return(ifelse(names_vec %in% names(name_map), name_map[names_vec], names_vec))
}
create_gt_table <- function(Data_forest_tmp, Factor_names, Factor_names_univariant) {
  withProgress(message = sample(nietzsche)[1], {
    Data_forest_tmp_table <- Data_forest_tmp %>%
      dplyr::mutate(EP_option = case_when(
        EP_option == 1 ~ "Yes",
        TRUE ~ "No"
      ))
    # --- Histology level reduction (top 30, others -> "OTHER") ---
    if (("Histology" %in% Factor_names_univariant) || ("Histology" %in% Factor_names)) {
      if ("Histology" %in% colnames(Data_forest_tmp_table)) {
        n_uniq_hist <- dplyr::n_distinct(Data_forest_tmp_table$Histology, na.rm = TRUE)
        if (n_uniq_hist >= 30) {
          top_hist <- Data_forest_tmp_table %>%
            dplyr::filter(!is.na(Histology)) %>%
            dplyr::count(Histology, sort = TRUE) %>%
            dplyr::slice_head(n = 30) %>%
            dplyr::pull(Histology)

          Data_forest_tmp_table <- Data_forest_tmp_table %>%
            dplyr::mutate(
              Histology = dplyr::if_else(
                is.na(Histology) | Histology %in% top_hist,
                as.character(Histology),
                "OTHER"
              )
            )
        }
      }
    }
    # --- end Histology reduction ---
    Data_forest_tmp_table$Histology <- factor(Data_forest_tmp_table$Histology)
    colnames(Data_forest_tmp_table) <- rename_factors_survival_CGP(colnames(Data_forest_tmp_table))
    Factor_names <- rename_factors_survival_CGP(Factor_names)
    Factor_names_univariant <- rename_factors_survival_CGP(Factor_names_univariant)
    Factor_names <- Factor_names[Factor_names != "CTx lines before CGP"]
    formula_str <- paste0("Surv(time_enroll_final, censor) ~ ", paste(paste0("`", Factor_names, "`"), collapse = " + "))
    linelistsurv_cox <- coxph(
      formula = as.formula(formula_str),
      data = Data_forest_tmp_table,
      control = coxph.control(iter.max = 50)
    )
    incProgress(1/4)
    univ_tab <- Data_forest_tmp_table %>%
      tbl_uvregression(
        method = coxph,
        y = Surv(time = time_enroll_final, event = censor),
        include = Factor_names_univariant,
        exponentiate = TRUE
      ) |>
      add_global_p() |>
      add_n(location = "level") |>
      add_nevent(location = "level") |>
      add_q() |>
      bold_p() |>
      apply_custom_format(columns = c("estimate", "conf.low", "conf.high", "p.value", "q.value")) |>
      bold_labels()
    incProgress(1/4)
    final_mv_reg <- linelistsurv_cox %>%
      stats::step(direction = "backward", trace = FALSE)
    incProgress(1/4)
  })
  if (!is.null(final_mv_reg$xlevels)) {
    mv_tab <- final_mv_reg |>
      tbl_regression(exponentiate = TRUE) |>
      add_global_p() |>
      add_q() |>
      bold_p() |>
      apply_custom_format(columns = c("estimate", "conf.low", "conf.high", "p.value", "q.value")) |>
      bold_labels()
    tbl_merge(
      tbls = list(univ_tab, mv_tab),
      tab_spanner = c("**Univariate**", "**Multivariable**")
    ) |>
      modify_caption("Hazard ratio for death after CGP (Only factors with >2 events observed)") |> as_gt()
  } else {
    univ_tab |>
      modify_caption("Hazard ratio for death after CGP, no significant factor in multivariable analysis (Only factors with >2 events observed)") |> as_gt()
  }
}

manual_one_hot <- function(data, target_cols) {
  result_data <- data
  for(col in target_cols) {
    if(col %in% names(data)) {
      unique_vals <- unique(data[[col]])
      unique_vals <- unique_vals[!is.na(unique_vals)]
      for(val in unique_vals) {
        new_col_name <- paste0(col, "_", val)
        result_data[[new_col_name]] <- as.numeric(data[[col]] == val)
      }
      result_data[[col]] <- NULL
    }
  }
  return(result_data)
}


format_numeric_columns <- function(df, digits = 3) {
  is_integer_column <- function(x) {
    is.integer(x) || (is.numeric(x) && all(x %% 1 == 0, na.rm = TRUE))
  }

  df[] <- lapply(df, function(col) {
    if (is.numeric(col) && !is_integer_column(col)) {
      formatted <- ifelse(
        is.na(col), NA_character_,
        ifelse(col < 0.001, "<0.001",
               ifelse(col > 1000, ">1000",
                      format(round(col, 3), nsmall = 3,  scientific = FALSE)))
      )
      return(formatted)
    } else {
      return(col)  # integer や文字列はそのまま
    }
  })
  return(df)
}


# create_datatable_with_confirm <- function(data,
#                                           page_length = 100,
#                                           scroll_y = "1000px",
#                                           buttons = c("csv", "excel"),
#                                           messages = list(
#                                             csv = "I will use the downloaded csv file in accordance with the terms and conditions.",
#                                             excel = "I will use the downloaded excel file in accordance with the terms and conditions."
#                                           )) {
#
#   # ボタン設定を動的に作成
#   button_configs <- list()
#
#   if ("csv" %in% buttons) {
#     button_configs <- append(button_configs, list(list(
#       extend = 'csv',
#       text = 'CSV',
#       action = DT::JS(paste0("
#         function(e, dt, node, config) {
#           if (confirm('", messages$csv, "')) {
#             $.fn.dataTable.ext.buttons.csvHtml5.action.call(this, e, dt, node, config);
#           }
#         }
#       "))
#     )))
#   }
#
#   if ("excel" %in% buttons) {
#     button_configs <- append(button_configs, list(list(
#       extend = 'excel',
#       text = 'Excel',
#       action = DT::JS(paste0("
#         function(e, dt, node, config) {
#           if (confirm('", messages$excel, "')) {
#             $.fn.dataTable.ext.buttons.excelHtml5.action.call(this, e, dt, node, config);
#           }
#         }
#       "))
#     )))
#   }
#
#   if ("copy" %in% buttons) {
#     button_configs <- append(button_configs, list(list(
#       extend = 'copy',
#       text = 'Copy',
#       action = DT::JS(paste0("
#         function(e, dt, node, config) {
#           if (confirm('", messages$copy, "')) {
#             $.fn.dataTable.ext.buttons.copyHtml5.action.call(this, e, dt, node, config);
#           }
#         }
#       "))
#     )))
#   }
#
#   DT::datatable(data,
#                 filter = 'top',
#                 extensions = c('Buttons'),
#                 options = list(
#                   pageLength = page_length,
#                   lengthMenu = c(100, 500, 1000, nrow(data)),
#                   scrollX = TRUE,
#                   scrollY = scroll_y,
#                   scrollCollapse = TRUE,
#                   dom = "Blfrtip",
#                   buttons = button_configs
#                 ))
# }

create_datatable_with_confirm <- function(data,
                                          message = NULL,
                                          page_length = 100,
                                          scroll_y = "1000px",
                                          buttons = c("csv", "excel", "copy"),
                                          messages = list(
                                            csv = "Please manage according to the personal information handling rules\\n=================================\\nPlease comply with the contents of the Agreement on Utilization of C-CAT Data and the service specification conformity disclosure of/for C-CAT Research-Use Portal site, and handle it properly.\\n=================================\\nBe careful not to leave the downloaded file on your computer",
                                            excel = "Please manage according to the personal information handling rules\\n=================================\\nPlease comply with the contents of the Agreement on Utilization of C-CAT Data and the service specification conformity disclosure of/for C-CAT Research-Use Portal site, and handle it properly.\\n=================================\\nBe careful not to leave the downloaded file on your computer",
                                            copy = "Please manage according to the personal information handling rules\\n=================================\\nPlease comply with the contents of the Agreement on Utilization of C-CAT Data and the service specification conformity disclosure of/for C-CAT Research-Use Portal site, and handle it properly.\\n=================================\\nBe careful not to leave the downloaded file on your computer"
                                          )) {

  # message をJSに安全に渡す（JSONとしてエスケープされた文字列）
  # 例: "FELIS CSV download about case summary"
  desc_json <- jsonlite::toJSON(
    if (is.null(message) || !nzchar(message)) NA_character_ else message,
    auto_unbox = TRUE
  )
  # --- 共通で使う fetch の JavaScript ロジックを定義 ---
  # btnText 変数（'csv', 'excel', 'copy' のいずれか）を使って通知する
  fetch_logic <- paste0("
    fetch(window.location.origin + '/@@/api/writelog', {
      method: 'POST',
      headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify({
        type: btnText,
        time: new Date().toISOString(),
        message: ", desc_json, "
      })
    }).catch(err => console.error('Download log fetch error:', err));
  ")

  # ボタン設定を動的に作成
  button_configs <- list()

  if ("csv" %in% buttons) {
    button_configs <- append(button_configs, list(list(
      extend = 'csv',
      text = 'CSV',
      action = DT::JS(paste0("
        function(e, dt, node, config) {
          // 1. 確認ダイアログを表示
          if (confirm('", messages$csv, "')) {
            // 2. [追加] 確認OKなら、通知を送信
            var btnText = 'csv';
            ", fetch_logic, "

            // 3. 本来のダウンロードアクションを実行
            $.fn.dataTable.ext.buttons.csvHtml5.action.call(this, e, dt, node, config);
          }
        }
      "))
    )))
  }

  if ("excel" %in% buttons) {
    button_configs <- append(button_configs, list(list(
      extend = 'excel',
      text = 'Excel',
      action = DT::JS(paste0("
        function(e, dt, node, config) {
          // 1. 確認ダイアログを表示
          if (confirm('", messages$excel, "')) {
            // 2. [追加] 確認OKなら、通知を送信
            var btnText = 'excel';
            ", fetch_logic, "

            // 3. 本来のダウンロードアクションを実行
            $.fn.dataTable.ext.buttons.excelHtml5.action.call(this, e, dt, node, config);
          }
        }
      "))
    )))
  }

  if ("copy" %in% buttons) {
    button_configs <- append(button_configs, list(list(
      extend = 'copy',
      text = 'Copy',
      action = DT::JS(paste0("
        function(e, dt, node, config) {
          // 1. 確認ダイアログを表示
          if (confirm('", messages$copy, "')) {
            // 2. [追加] 確認OKなら、通知を送信
            var btnText = 'copy';
            ", fetch_logic, "

            // 3. 本来のコピーアクションを実行
            $.fn.dataTable.ext.buttons.copyHtml5.action.call(this, e, dt, node, config);
          }
        }
      "))
    )))
  }

  DT::datatable(data,
                filter = 'top',
                extensions = c('Buttons'),
                options = list(
                  pageLength = page_length,
                  lengthMenu = c(100, 500, 1000, nrow(data)),
                  scrollX = TRUE,
                  scrollY = scroll_y,
                  scrollCollapse = TRUE,
                  dom = "Blfrtip",
                  buttons = button_configs
                  # initComplete は使いません
                ))
}


clear_reactive_data <- function(reactive_obj, verbose = FALSE) {

  tryCatch({
    isolate({
      # 全要素を取得
      all_names <- names(reactive_obj)

      if(verbose) cat("Found", length(all_names), "elements to clear\n")

      if(length(all_names) > 0) {
        # 各要素をNULLに設定
        for(nm in all_names) {
          tryCatch({
            reactive_obj[[nm]] <- NULL
            if(verbose) cat("Cleared:", nm, "\n")
          }, error = function(e) {
            if(verbose) cat("Failed to clear:", nm, "-", e$message, "\n")
          })
        }

        # 内部環境の安全なチェックと削除
        tryCatch({
          # reactiveValuesの内部構造を安全に確認
          if(is.environment(reactive_obj) && ".values" %in% names(reactive_obj)) {
            internal_env <- reactive_obj$.values

            if(is.environment(internal_env)) {
              internal_names <- ls(envir = internal_env)

              if(length(internal_names) > 0) {
                if(verbose) cat("Found", length(internal_names), "internal objects\n")

                # 内部オブジェクトを削除
                rm(list = internal_names, envir = internal_env)
                if(verbose) cat("Removed", length(internal_names), "internal objects\n")
              }
            }
          }
        }, error = function(e) {
          if(verbose) cat("Could not access internal environment:", e$message, "\n")
          # 内部環境にアクセスできなくても、API層のクリアは成功しているので続行
        })
      } else {
        if(verbose) cat("No elements to clear\n")
      }
    })

    # 強制ガベージコレクション
    gc(verbose = FALSE)

    if(verbose) cat("Clear operation completed successfully\n")
    return(TRUE)

  }, error = function(e) {
    if(verbose) cat("Error in clear_reactive_data:", e$message, "\n")
    gc(verbose = FALSE)
    return(FALSE)
  })
}

htmlOutputWithPopover <- function(outputId,
                                  popoverTitle,
                                  popoverContent,
                                  helpId = NULL,
                                  placement = "top",
                                  trigger = "hover") {

  if (is.null(helpId)) {
    helpId <- paste0("help_", outputId)
  }

  # HTMLコンテンツをJavaScript用にエスケープ
  escaped_content <- gsub('"', '\\"', gsub("'", "\\'", as.character(popoverContent)))
  escaped_title <- gsub('"', '\\"', gsub("'", "\\'", popoverTitle))

  list(
    div(style = "display: flex; align-items: flex-start; width: 100%; min-height: 25px;",
        div(style = "flex: 1; min-width: 0; padding-right: 5px;",
            htmlOutput(outputId)
        ),
        div(style = "flex-shrink: 0; width: 20px; height: 20px; margin-top: 2px;",
            actionButton(helpId, "",
                         icon = icon("question-circle"),
                         style = "border: none; background: transparent; color: #3c8dbc; width: 18px; height: 18px; padding: 0; font-size: 13px;")
        )
    ),

    # JavaScriptでHTMLポップオーバーを初期化
    tags$script(HTML(paste0("
      $(document).ready(function(){
        $('#", helpId, "').popover({
          title: '", escaped_title, "',
          content: '", escaped_content, "',
          html: true,
          placement: '", placement, "',
          trigger: '", trigger, "',
          container: 'body',
          template: '<div class=\"popover\" role=\"tooltip\"><div class=\"arrow\"></div><h3 class=\"popover-title\"></h3><div class=\"popover-content\"></div></div>'
        });
      });
    ")))
  )
}

BootNoSet = function(Data){
  return(as.integer(min(500, max(25, as.integer(50000 / nrow(Data))*50))/BootNoVar))
}

run_crossval <- function(cross_validation_samples, formula_treat, Penal, Toler) {
  try(silent = FALSE,
      cross_validation_samples %>%
        dplyr::mutate(
          glm_analysis = lapply(splits, function(split) {
            lrm(formula = formula_treat,
                data = rsample::analysis(split),
                x = TRUE, y = TRUE,
                penalty = Penal, tol = Toler)
          })
        )
  )
}

optimize_data_datatable <- function(Data_forest_tmp_, input) {
  # EP_treat == 1の行を事前に計算
  ep_treat_mask <- Data_forest_tmp_$EP_treat == 1

  # 削除する列を収集するベクトル
  cols_to_remove <- character()

  # 各列の処理を関数化
  check_and_process <- function(col_name, value_check, alt_value = NULL) {
    if (!(col_name %in% names(Data_forest_tmp_))) return(NULL)

    col_data <- Data_forest_tmp_[[col_name]]
    if (is.null(alt_value)) {
      # 単純な削除チェック
      count_match <- sum(col_data == value_check & ep_treat_mask, na.rm = TRUE)
      count_not_match <- sum(col_data != value_check & ep_treat_mask, na.rm = TRUE)

      if (count_match < 3 || count_not_match < 3) {
        return(col_name)
      }
    }
    return(NULL)
  }

  # PS列の処理
  if ("PS" %in% names(Data_forest_tmp_)) {
    ps_data <- Data_forest_tmp_$PS
    count_0 <- sum(ps_data == "0" & ep_treat_mask)
    count_not_0 <- sum(ps_data != "0" & ep_treat_mask)
    count_2_4 <- sum(ps_data == "2_4" & ep_treat_mask)

    if (count_0 < 3 || count_not_0 < 3) {
      cols_to_remove <- c(cols_to_remove, "PS")
    } else if (count_2_4 < 3) {
      Data_forest_tmp_$PS[ps_data %in% c("1", "2_4")] <- "1_4"
    }
  }

  # 単純な二値変数の処理
  binary_cols <- list(
    c("Smoking_history", "Yes"),
    c("Alcoholic_history", "Yes"),
    c("time_diagnosis_enroll", ">1-year"),
    c("Lymph_met", "Yes"),
    c("Lung_met", "Yes"),
    c("Brain_met", "Yes"),
    c("Bone_met", "Yes"),
    c("Liver_met", "Yes")
  )

  for (col_info in binary_cols) {
    col_to_remove <- check_and_process(col_info[1], col_info[2])
    if (!is.null(col_to_remove)) {
      cols_to_remove <- c(cols_to_remove, col_to_remove)
    }
  }

  # Panel列の処理
  if ("Panel" %in% names(Data_forest_tmp_)) {
    panel_data <- Data_forest_tmp_$Panel
    solid_panels <- c("NCC OncoPanel", "FoundationOne CDx", "GenMineTOP")

    if (any(sapply(solid_panels, function(p) sum(panel_data == p & ep_treat_mask) < 3))) {
      Data_forest_tmp_$Panel[panel_data %in% solid_panels] <- "Solid"

      if (length(unique(Data_forest_tmp_$Panel)) > 1) {
        panel_table <- table(Data_forest_tmp_$Panel)
        ref_level <- names(sort(panel_table, decreasing = TRUE))[1]
        Data_forest_tmp_$Panel <- relevel(factor(Data_forest_tmp_$Panel), ref = ref_level)
      }

      count_solid <- sum(Data_forest_tmp_$Panel == "Solid" & ep_treat_mask)
      count_not_solid <- sum(Data_forest_tmp_$Panel != "Solid" & ep_treat_mask)

      if (count_solid < 3 || count_not_solid < 3) {
        cols_to_remove <- c(cols_to_remove, "Panel")
      }
    }
  }

  # Enroll_date列の処理
  if ("Enroll_date" %in% names(Data_forest_tmp_)) {
    enroll_data <- Data_forest_tmp_$Enroll_date
    count_2019 <- sum(enroll_data == "2019" & ep_treat_mask, na.rm = TRUE)

    if (count_2019 < 3) {
      Data_forest_tmp_$Enroll_date[enroll_data %in% c("2019", "2020")] <- "2019/20"

      count_2019_20 <- sum(Data_forest_tmp_$Enroll_date == "2019/20" & ep_treat_mask, na.rm = TRUE)

      if (count_2019_20 < 3) {
        Data_forest_tmp_$Enroll_date[Data_forest_tmp_$Enroll_date %in% c("2019/20", "2021")] <- "2019-21"

        count_2019_21 <- sum(Data_forest_tmp_$Enroll_date == "2019-21" & ep_treat_mask, na.rm = TRUE)
        count_not_2019_21 <- sum(Data_forest_tmp_$Enroll_date != "2019-21" & ep_treat_mask, na.rm = TRUE)

        if (count_2019_21 < 3 || count_not_2019_21 < 3) {
          cols_to_remove <- c(cols_to_remove, "Enroll_date")
        }
      } else {
        count_not_2019_20 <- sum(Data_forest_tmp_$Enroll_date != "2019/20" & ep_treat_mask, na.rm = TRUE)
        if (count_not_2019_20 < 3) {
          cols_to_remove <- c(cols_to_remove, "Enroll_date")
        }
      }
    } else {
      count_not_2019 <- sum(enroll_data != "2019" & ep_treat_mask, na.rm = TRUE)
      if (count_not_2019 < 3) {
        cols_to_remove <- c(cols_to_remove, "Enroll_date")
      }
    }
  }

  # Age列の処理
  if ("Age" %in% names(Data_forest_tmp_)) {
    age_data <- Data_forest_tmp_$Age
    count_younger <- sum(age_data == "Younger" & ep_treat_mask)
    count_older <- sum(age_data == "Older" & ep_treat_mask)

    if (count_younger < 3 || count_older < 3) {
      cols_to_remove <- c(cols_to_remove, "Age")
    }
  }

  # Best_effect列の処理（注：元のコードにバグがある可能性）
  if ("Best_effect" %in% names(Data_forest_tmp_)) {
    best_data <- Data_forest_tmp_$Best_effect
    count_sd <- sum(best_data == "SD" & ep_treat_mask, na.rm = TRUE)
    count_not_sd <- sum(best_data != "SD" & ep_treat_mask, na.rm = TRUE)

    if (count_sd < 3) {
      cols_to_remove <- c(cols_to_remove, "Best_effect")
    } else if (count_not_sd < 3) {
      # 元のコードはselect(Best_effect)となっているが、これは他の列を削除する意味？
      # バグの可能性があるため、コメントとして残す
      # Data_forest_tmp_ <- Data_forest_tmp_[, "Best_effect", drop = FALSE]
      cols_to_remove <- c(cols_to_remove, setdiff(names(Data_forest_tmp_), c("Best_effect", "EP_treat")))
    }
  }

  # Sex列の処理
  if ("Sex" %in% names(Data_forest_tmp_)) {
    sex_data <- Data_forest_tmp_$Sex
    count_male <- sum(sex_data == "Male" & ep_treat_mask, na.rm = TRUE)
    count_not_male <- sum(sex_data != "Male" & ep_treat_mask, na.rm = TRUE)

    if (count_male < 3 || count_not_male < 3) {
      cols_to_remove <- c(cols_to_remove, "Sex")
    }
  }

  # Lines列の処理
  if ("Lines" %in% names(Data_forest_tmp_)) {
    lines_data <- Data_forest_tmp_$Lines
    count_1 <- sum(lines_data == "1" & ep_treat_mask)
    count_not_1 <- sum(lines_data != "1" & ep_treat_mask)
    count_0 <- sum(lines_data == "0" & ep_treat_mask)
    count_2_plus <- sum(lines_data == "2~" & ep_treat_mask)

    if (count_1 < 3 || count_not_1 < 3) {
      cols_to_remove <- c(cols_to_remove, "Lines")
    } else if (count_0 < 3) {
      Data_forest_tmp_$Lines[lines_data %in% c("0", "1")] <- "0~1"
    } else if (count_2_plus < 3) {
      Data_forest_tmp_$Lines[lines_data %in% c("2~", "1")] <- "1~"
    }
  }

  # 条件付き列の処理
  if (!is.null(input$HER2) && input$HER2 != "No" && "HER2_IHC" %in% names(Data_forest_tmp_)) {
    her2_data <- Data_forest_tmp_$HER2_IHC
    count_pos <- sum(her2_data == "Positive" & ep_treat_mask)
    count_not_pos <- sum(her2_data != "Positive" & ep_treat_mask)

    if (count_pos < 3 || count_not_pos < 3) {
      cols_to_remove <- c(cols_to_remove, "HER2_IHC")
    }
  }

  if (!is.null(input$MSI) && input$MSI != "No" && "MSI_PCR" %in% names(Data_forest_tmp_)) {
    msi_data <- Data_forest_tmp_$MSI_PCR
    count_pos <- sum(msi_data == "Positive" & ep_treat_mask)
    count_not_pos <- sum(msi_data != "Positive" & ep_treat_mask)

    if (count_pos < 3 || count_not_pos < 3) {
      cols_to_remove <- c(cols_to_remove, "MSI_PCR")
    }
  }

  if (!is.null(input$MMR) && input$MMR != "No" && "MMR_IHC" %in% names(Data_forest_tmp_)) {
    mmr_data <- Data_forest_tmp_$MMR_IHC
    count_dmmr <- sum(mmr_data == "dMMR" & ep_treat_mask)
    count_not_dmmr <- sum(mmr_data != "dMMR" & ep_treat_mask)

    if (count_dmmr < 3 || count_not_dmmr < 3) {
      cols_to_remove <- c(cols_to_remove, "MMR_IHC")
    }
  }

  # 一度にすべての列を削除
  if (length(cols_to_remove) > 0) {
    cols_to_keep <- setdiff(names(Data_forest_tmp_), cols_to_remove)
    Data_forest_tmp_ <- Data_forest_tmp_[,  ..cols_to_keep, drop = FALSE]
  }

  return(Data_forest_tmp_)
}

# Helper function to calculate age-stratified IPTW
calculate_iptw_age <- function(data, ref_surv_list, time_var = "time_pre", age_var = "症例.基本情報.年齢") {
  init_pop <- 10000
  max_years <- 10
  bin_width <- 0.5 # 6-month window
  breaks <- seq(0, max_years, by = bin_width)
  n_bins <- length(breaks) - 1

  # Categorize age and assign to 6-month time bins
  data <- data %>%
    dplyr::mutate(
      age_num = as.numeric(!!sym(age_var)),
      age_class = dplyr::case_when(
        age_num < 40 ~ "40未満",
        age_num < 50 ~ "40代",
        age_num < 60 ~ "50代",
        age_num < 70 ~ "60代",
        age_num < 80 ~ "70代",
        age_num >= 80 ~ "80以上",
        TRUE ~ "Unknown"
      ),
      time_years = !!sym(time_var) / 365.25,
      time_bin = ceiling(time_years / bin_width),
      time_bin = ifelse(time_bin > n_bins, n_bins, time_bin),
      time_bin = ifelse(time_bin == 0, 1, time_bin) # Safety for time=0
    )

  # Build Person-Time (PT) reference table
  pt_table <- expand.grid(age_class = names(ref_surv_list), time_bin = 1:n_bins, stringsAsFactors = FALSE)
  pt_table$pt_ref <- 0

  t_points <- 1:5 # We have data for years 1, 2, 3, 4, 5

  for(ag in names(ref_surv_list)) {
    surv_rates <- ref_surv_list[[ag]]

    if(length(surv_rates) >= 5) {
      # Convert % to probabilities and avoid log(0) bounds
      S_t <- surv_rates[1:5] / 100
      S_t <- pmax(pmin(S_t, 0.999), 0.001)

      # Log-logistic linearization: log(1/S(t) - 1) = p * log(lambda) + p * log(t)
      y <- log(1/S_t - 1)
      x <- log(t_points)

      # Fit linear model to find parameters
      fit <- lm(y ~ x)
      p <- coef(fit)[2]
      p_log_lambda <- coef(fit)[1]
      lambda <- exp(p_log_lambda / p)

      # Define smooth Log-logistic survival function
      S_fit <- function(t_y) {
        1 / (1 + (lambda * t_y)^p)
      }

      # Calculate Expected PT for each 6-month window using trapezoidal rule
      pt_bins <- numeric(n_bins)
      for(i in 1:n_bins) {
        t_start <- breaks[i]
        t_end <- breaks[i+1]
        pt_bins[i] <- init_pop * (S_fit(t_start) + S_fit(t_end)) / 2 * bin_width
      }
      pt_table[pt_table$age_class == ag, "pt_ref"] <- pt_bins
    }
  }

  # Count actual N in CGP data per age class and 6-month bin
  bin_counts <- data %>% dplyr::count(age_class, time_bin, name = "N_cgp")

  # Calculate IPTW (Weight = PT / N_cgp)
  data <- data %>%
    dplyr::left_join(pt_table, by = c("age_class", "time_bin")) %>%
    dplyr::left_join(bin_counts, by = c("age_class", "time_bin")) %>%
    dplyr::mutate(
      raw_weight = ifelse(!is.na(N_cgp) & N_cgp > 0 & !is.na(pt_ref), pt_ref / N_cgp, 0)
    )

  # =========================================================================
  # 極端な重み（外れ値）による分散の爆発を防ぐための2.5%〜97.5%トリミング
  # =========================================================================
  if (any(data$raw_weight > 0, na.rm = TRUE)) {
    lower_bound <- quantile(data$raw_weight[data$raw_weight > 0], 0.025, na.rm = TRUE)
    upper_bound <- quantile(data$raw_weight[data$raw_weight > 0], 0.975, na.rm = TRUE)

    data <- data %>%
      dplyr::mutate(
        raw_weight = ifelse(raw_weight > 0 & raw_weight < lower_bound, lower_bound, raw_weight),
        raw_weight = ifelse(raw_weight > upper_bound, upper_bound, raw_weight)
      )
  }
  # Stabilize weights (mean = 1)
  mean_w <- mean(data$raw_weight[data$raw_weight > 0], na.rm = TRUE)
  data$iptw <- ifelse(data$raw_weight > 0, data$raw_weight / mean_w, 1.0)

  return(data)
}
