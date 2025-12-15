observeEvent(input$start_prediction, {
  withProgress(message = "Prediction model loading", {
    incProgress(1 / 13)
    if(file.exists(paste0("source/m2", input$prediction_var, ".qs"))){
      OUTPUT_DATA$m2 = QS_READ(nthreads = max(1, parallel::detectCores() - 1, na.rm = TRUE), file=paste0("source/m2", input$prediction_var, ".qs"))
    }
    if(file.exists(paste0("source/lgb_m2", input$prediction_var, ".qs"))){
      OUTPUT_DATA$lgb_m2 = QS_READ(nthreads = max(1, parallel::detectCores() - 1, na.rm = TRUE), file=paste0("source/lgb_m2", input$prediction_var, ".qs"))
    }
    if(file.exists(paste0("source/rf_m2", input$prediction_var, ".qs"))){
      OUTPUT_DATA$rf_m2 = QS_READ(nthreads = max(1, parallel::detectCores() - 1, na.rm = TRUE), file=paste0("source/rf_m2", input$prediction_var, ".qs"))
    }
  })
  req(OUTPUT_DATA$m2, input$predict_age, input$predict_organ, input$predict_Panel)
  withProgress(message = sample(nietzsche)[1], {
    incProgress(1 / 13)
    predict_age = ifelse(input$predict_age <= input$mid_age, "Younger", "Older")
    
    OUTPUT_DATA$predict_input_data = data.frame(
      Age = predict_age,
      Lymph_met = input$predict_Lymph_met,
      Smoking_history = input$predict_Smoking_history,
      Alcoholic_history = input$predict_Alcoholic_history,
      Enroll_date = input$predict_Enroll_date,
      time_diagnosis_enroll = input$predict_time_diagnosis_enroll,
      Brain_met = input$predict_Brain_met,
      Lung_met = input$predict_Lung_met,
      Bone_met = input$predict_Bone_met,
      Liver_met = input$predict_Liver_met,
      Sex = input$predict_Sex,
      Histology = input$predict_organ,
      PS = input$predict_PS,
      Lines = input$predict_Lines,
      Best_effect = input$predict_Best_effect,
      Panel = input$predict_Panel
    )
    OUTPUT_DATA$predict_input_data_ML = data.frame(
      Age = if_else(input$predict_age <= input$mid_age, "Younger", "Older"),
      Lymph_met = input$predict_Lymph_met,
      Smoking_history = input$predict_Smoking_history,
      Alcoholic_history = input$predict_Alcoholic_history,
      Enroll_date = input$predict_Enroll_date,
      time_diagnosis_enroll = input$predict_time_diagnosis_enroll,
      Brain_met = input$predict_Brain_met,
      Lung_met = input$predict_Lung_met,
      Bone_met = input$predict_Bone_met,
      Liver_met = input$predict_Liver_met,
      Sex = input$predict_Sex,
      Histology = input$predict_organ_raw,
      PS = input$predict_PS_raw,
      Lines = input$predict_Lines_raw,
      Best_effect = input$predict_Best_effect,
      Panel = input$predict_Panel
    )
    prediction_var <- if(exists("input") && !is.null(input$prediction_var)) {
      input$prediction_var
    } else {
      NULL
    }
    
    # カテゴリ変数の定義
    category_values <- list(
      "PS" = c("0", "1", "2_4"),
      "Panel" = c("NCC OncoPanel", "FoundationOne CDx", "GenMineTOP", "Liquid"),
      "Best_effect" = c("CR", "PR", "SD", "PD", "Unknown", "CTx_naive"),
      "Lines" = c("0", "1", "2", "3~"),
      "Enroll_date" = c("2019", "2020", "2021", "2022", "2023", "2024", "2025")
    )
    
    # Histologyの動的読み込み
    if(!is.null(prediction_var)) {
      histology_file <- paste0("source/Organ_data", prediction_var, ".rda")
      if(file.exists(histology_file)) {
        load(histology_file)
        dynamic_data <- get(paste0("Organ_data", prediction_var))
        category_values[["Histology"]] <- dynamic_data$Organ_type_raw
      }
    }
    
    # データの正規化
    data <- OUTPUT_DATA$predict_input_data_ML
    
    # PS値の正規化
    if("PS" %in% names(data)) {
      data$PS <- ifelse(data$PS %in% c(2, 3, 4), "2_4", as.character(data$PS))
    }
    
    # Lines値の正規化
    if("Lines" %in% names(data)) {
      data$Lines <- ifelse(data$Lines >= 2, "2~", as.character(data$Lines))
    }
    
    # ワンホットエンコーディング
    for(var in names(category_values)) {
      if(var %in% names(data)) {
        for(val in category_values[[var]]) {
          col_name <- paste0(var, "_", gsub("[^A-Za-z0-9]", "_", val))
          col_name <- gsub("_+", "_", col_name)
          col_name <- gsub("^_|_$", "", col_name)
          data[[col_name]] <- as.numeric(data[[var]] == val)
        }
        data[[var]] <- NULL
      }
    }
    
    OUTPUT_DATA$predict_input_data_ML <- data
    if(!is.null(OUTPUT_DATA$m2)){
      OUTPUT_DATA$predicted_score = predict(OUTPUT_DATA$m2, OUTPUT_DATA$predict_input_data, se.fit=TRUE)
    }
    if(!is.null(OUTPUT_DATA$lgb_m2)){
      OUTPUT_DATA$predicted_score_lgb = predict(OUTPUT_DATA$lgb_m2, OUTPUT_DATA$predict_input_data_ML, type = "prob")
    }
    if(!is.null(OUTPUT_DATA$rf_m2)){
      OUTPUT_DATA$predicted_score_rf = predict(OUTPUT_DATA$rf_m2, OUTPUT_DATA$predict_input_data_ML, type = "prob")
      OUTPUT_DATA$predicted_score_rf_CI = predict(OUTPUT_DATA$rf_m2, OUTPUT_DATA$predict_input_data_ML, type = "conf_int")
    }
    incProgress(1 / 13)
  })
})

output$prediction_nomogram = renderUI({
  req(OUTPUT_DATA$predicted_score)
  if(!is.null(OUTPUT_DATA$predicted_score) &
     !is.null(OUTPUT_DATA$predicted_score_lgb) &
     !is.null(OUTPUT_DATA$predicted_score_rf)){
    HTML(paste0("Predicted rate", tags$br(),
           "Logistic regression: ", format_p(exp(OUTPUT_DATA$predicted_score[[1]][[1]])/(1+exp(OUTPUT_DATA$predicted_score[[1]][[1]]))*100, digits=2),
           "% (95%CI: ",
           format_p(exp(OUTPUT_DATA$predicted_score[[1]][[1]] - 1.96 * predict(OUTPUT_DATA$m2, OUTPUT_DATA$predict_input_data, se.fit=TRUE)[[2]][[1]])/(1+exp(OUTPUT_DATA$predicted_score[[1]][[1]] - 1.96 * predict(OUTPUT_DATA$m2, OUTPUT_DATA$predict_input_data, se.fit=TRUE)[[2]][[1]]))*100, digits=2), "-",
           format_p(exp(OUTPUT_DATA$predicted_score[[1]][[1]] + 1.96 * predict(OUTPUT_DATA$m2, OUTPUT_DATA$predict_input_data, se.fit=TRUE)[[2]][[1]])/(1+exp(OUTPUT_DATA$predicted_score[[1]][[1]] + 1.96 * predict(OUTPUT_DATA$m2, OUTPUT_DATA$predict_input_data, se.fit=TRUE)[[2]][[1]]))*100, digits=2),")", tags$br(),
           "LightGBM: ", format_p(OUTPUT_DATA$predicted_score_lgb[2]*100, digits=2), "%", tags$br(),
           "Random forest: ", format_p(predicted_score_rf[2]*100, digits=2), "% (95%CI: ",
           format_p(OUTPUT_DATA$predicted_score_rf_CI[[3]]*100, digits=2), "-",
           format_p(OUTPUT_DATA$predicted_score_rf_CI[[4]]*100, digits=2),")", tags$br()))
  } else if(!is.null(OUTPUT_DATA$predicted_score) &
            !is.null(OUTPUT_DATA$predicted_score_lgb)){
    HTML(paste0("Predicted rate", tags$br(),
           "Logistic regression: ", format_p(exp(OUTPUT_DATA$predicted_score[[1]][[1]])/(1+exp(OUTPUT_DATA$predicted_score[[1]][[1]]))*100, digits=2),
           "% (95%CI: ",
           format_p(exp(OUTPUT_DATA$predicted_score[[1]][[1]] - 1.96 * predict(OUTPUT_DATA$m2, OUTPUT_DATA$predict_input_data, se.fit=TRUE)[[2]][[1]])/(1+exp(OUTPUT_DATA$predicted_score[[1]][[1]] - 1.96 * predict(OUTPUT_DATA$m2, OUTPUT_DATA$predict_input_data, se.fit=TRUE)[[2]][[1]]))*100, digits=2), "-",
           format_p(exp(OUTPUT_DATA$predicted_score[[1]][[1]] + 1.96 * predict(OUTPUT_DATA$m2, OUTPUT_DATA$predict_input_data, se.fit=TRUE)[[2]][[1]])/(1+exp(OUTPUT_DATA$predicted_score[[1]][[1]] + 1.96 * predict(OUTPUT_DATA$m2, OUTPUT_DATA$predict_input_data, se.fit=TRUE)[[2]][[1]]))*100, digits=2),")", tags$br(),
           "LightGBM: ", format_p(OUTPUT_DATA$predicted_score_lgb[2]*100, digits=2), "%", tags$br()))
  } else if(!is.null(OUTPUT_DATA$predicted_score) &
            !is.null(OUTPUT_DATA$predicted_score_rf)){
    HTML(paste0("Predicted rate", tags$br(),
           "Logistic regression: ", format_p(exp(OUTPUT_DATA$predicted_score[[1]][[1]])/(1+exp(OUTPUT_DATA$predicted_score[[1]][[1]]))*100, digits=2),
           "% (95%CI: ",
           format_p(exp(OUTPUT_DATA$predicted_score[[1]][[1]] - 1.96 * predict(OUTPUT_DATA$m2, OUTPUT_DATA$predict_input_data, se.fit=TRUE)[[2]][[1]])/(1+exp(OUTPUT_DATA$predicted_score[[1]][[1]] - 1.96 * predict(OUTPUT_DATA$m2, OUTPUT_DATA$predict_input_data, se.fit=TRUE)[[2]][[1]]))*100, digits=2), "-",
           format_p(exp(OUTPUT_DATA$predicted_score[[1]][[1]] + 1.96 * predict(OUTPUT_DATA$m2, OUTPUT_DATA$predict_input_data, se.fit=TRUE)[[2]][[1]])/(1+exp(OUTPUT_DATA$predicted_score[[1]][[1]] + 1.96 * predict(OUTPUT_DATA$m2, OUTPUT_DATA$predict_input_data, se.fit=TRUE)[[2]][[1]]))*100, digits=2),")", tags$br(),
           "Random forest: ", format_p(OUTPUT_DATA$predicted_score_rf[2]*100, digits=2), "% (95%CI: ",
           format_p(OUTPUT_DATA$predicted_score_rf_CI[[3]]*100, digits=2), "-",
           format_p(OUTPUT_DATA$predicted_score_rf_CI[[4]]*100, digits=2),")", tags$br()))
  } else if(!is.null(OUTPUT_DATA$predicted_score)) {
    HTML(paste0("Predicted rate", tags$br(),
           "Logistic regression: ", format_p(exp(OUTPUT_DATA$predicted_score[[1]][[1]])/(1+exp(OUTPUT_DATA$predicted_score[[1]][[1]]))*100, digits=2),
           "% (95%CI: ",
           format_p(exp(OUTPUT_DATA$predicted_score[[1]][[1]] - 1.96 * predict(OUTPUT_DATA$m2, OUTPUT_DATA$predict_input_data, se.fit=TRUE)[[2]][[1]])/(1+exp(OUTPUT_DATA$predicted_score[[1]][[1]] - 1.96 * predict(OUTPUT_DATA$m2, OUTPUT_DATA$predict_input_data, se.fit=TRUE)[[2]][[1]]))*100, digits=2), "-",
           format_p(exp(OUTPUT_DATA$predicted_score[[1]][[1]] + 1.96 * predict(OUTPUT_DATA$m2, OUTPUT_DATA$predict_input_data, se.fit=TRUE)[[2]][[1]])/(1+exp(OUTPUT_DATA$predicted_score[[1]][[1]] + 1.96 * predict(OUTPUT_DATA$m2, OUTPUT_DATA$predict_input_data, se.fit=TRUE)[[2]][[1]]))*100, digits=2),")", tags$br()))
  }
})