observeEvent(input$start_prediction, {
  withProgress(message = "Prediction model loading", {
    incProgress(1 / 13)
    if(file.exists(paste0(tempdir(), "/m2", input$prediction_var, ".qs"))){
      OUTPUT_DATA$m2 = QS_READ(nthreads = max(1, parallel::detectCores() - 1, na.rm = TRUE), file=paste0(tempdir(), "/m2", input$prediction_var, ".qs"))
    }
    if(file.exists(paste0(tempdir(), "/lgb_m2", input$prediction_var, ".qs"))){
      OUTPUT_DATA$lgb_m2 = QS_READ(nthreads = max(1, parallel::detectCores() - 1, na.rm = TRUE), file=paste0(tempdir(), "/lgb_m2", input$prediction_var, ".qs"))
    }
    if(file.exists(paste0(tempdir(), "/rf_m2", input$prediction_var, ".qs"))){
      OUTPUT_DATA$rf_m2 = QS_READ(nthreads = max(1, parallel::detectCores() - 1, na.rm = TRUE), file=paste0(tempdir(), "/rf_m2", input$prediction_var, ".qs"))
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
    OUTPUT_DATA$predict_input_data_ML = OUTPUT_DATA$predict_input_data
    prediction_var <- if(exists("input") && !is.null(input$prediction_var)) {
      input$prediction_var
    } else {
      NULL
    }
    # カテゴリ変数の定義
    category_values <- list(
      "Best_effect" = c("CR", "PR", "SD", "PD", "Unknown", "CTx_naive"),
      "Age" = c("Older", "Younger"),
      "Lymph_met" = c("Yes", "No"),
      "Brain_met" = c("Yes", "No"),
      "Lung_met" = c("Yes", "No"),
      "Bone_met" = c("Yes", "No"),
      "Liver_met" = c("Yes", "No"),
      "Sex" = c("Female", "Male"),
      "Smoking_history" = c("Yes", "No"),
      "Alcoholic_history" = c("Yes", "No"),
      "time_diagnosis_enroll" = c("<=1-year", ">1-year"),

      "Histology" = "Not significant",
      "PS" = "Not significant",
      "Panel" = "Not significant",
      "Lines" = "Not significant"
    )

    # Histologyの動的読み込み
    if(!is.null(prediction_var)) {
      histology_file <- paste0(tempdir(), "/Organ_data", prediction_var, ".rda")
      if(file.exists(histology_file)) {
        load(histology_file)
        dynamic_data <- get(paste0("Organ_data", prediction_var))
        category_values[["Histology"]] <- dynamic_data$Organ_type_raw
      }
      PS_file <- paste0(tempdir(), "/PS_data", prediction_var, ".rda")
      if(file.exists(PS_file)) {
        load(PS_file)
        dynamic_data <- get(paste0("PS_data", prediction_var))
        category_values[["PS"]] <- dynamic_data$PS_type
      }
      Panel_file <- paste0(tempdir(), "/Panel_data", prediction_var, ".rda")
      if(file.exists(Panel_file)) {
        load(Panel_file)
        dynamic_data <- get(paste0("Panel_data", prediction_var))
        category_values[["Panel"]] <- dynamic_data$Panel_type
      }
      Lines_file <- paste0(tempdir(), "/Lines_data", prediction_var, ".rda")
      if(file.exists(Lines_file)) {
        load(Lines_file)
        dynamic_data <- get(paste0("Lines_data", prediction_var))
        category_values[["Lines"]] <- dynamic_data$Lines_type
      }
    }

    # データの正規化
    data <- OUTPUT_DATA$predict_input_data_ML

    for(var in names(category_values)) {
      if(var %in% names(data)) {
        for(val in category_values[[var]]) {
          # Preserve <=, >=, <, >, and spaces while replacing other special characters
          col_name <- paste0(var, "_", val)
          col_name <- gsub("[^A-Za-z0-9<>=_ ]", "_", col_name)  # Keep comparison operators and spaces
          col_name <- gsub("_+", "_", col_name)  # Remove multiple underscores
          col_name <- gsub("^_|_$", "", col_name)  # Remove leading/trailing underscores
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
      OUTPUT_DATA$predicted_score_lgb = predict(OUTPUT_DATA$lgb_m2, OUTPUT_DATA$predict_input_data, type = "prob")
    }
    if(!is.null(OUTPUT_DATA$rf_m2)){
      OUTPUT_DATA$predicted_score_rf = predict(OUTPUT_DATA$rf_m2, OUTPUT_DATA$predict_input_data, type = "prob")
      OUTPUT_DATA$predicted_score_rf_CI = predict(OUTPUT_DATA$rf_m2, OUTPUT_DATA$predict_input_data, type = "conf_int")
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
           "Random forest: ", format_p(OUTPUT_DATA$predicted_score_rf[2]*100, digits=2), "% (95%CI: ",
           format_p(OUTPUT_DATA$predicted_score_rf_CI[[3]][[1]]*100, digits=2), "-",
           format_p(OUTPUT_DATA$predicted_score_rf_CI[[4]][[1]]*100, digits=2),")", tags$br()))
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
           format_p(OUTPUT_DATA$predicted_score_rf_CI[[3]][[1]]*100, digits=2), "-",
           format_p(OUTPUT_DATA$predicted_score_rf_CI[[4]][[1]]*100, digits=2),")", tags$br()))
  } else if(!is.null(OUTPUT_DATA$predicted_score)) {
    HTML(paste0("Predicted rate", tags$br(),
           "Logistic regression: ", format_p(exp(OUTPUT_DATA$predicted_score[[1]][[1]])/(1+exp(OUTPUT_DATA$predicted_score[[1]][[1]]))*100, digits=2),
           "% (95%CI: ",
           format_p(exp(OUTPUT_DATA$predicted_score[[1]][[1]] - 1.96 * predict(OUTPUT_DATA$m2, OUTPUT_DATA$predict_input_data, se.fit=TRUE)[[2]][[1]])/(1+exp(OUTPUT_DATA$predicted_score[[1]][[1]] - 1.96 * predict(OUTPUT_DATA$m2, OUTPUT_DATA$predict_input_data, se.fit=TRUE)[[2]][[1]]))*100, digits=2), "-",
           format_p(exp(OUTPUT_DATA$predicted_score[[1]][[1]] + 1.96 * predict(OUTPUT_DATA$m2, OUTPUT_DATA$predict_input_data, se.fit=TRUE)[[2]][[1]])/(1+exp(OUTPUT_DATA$predicted_score[[1]][[1]] + 1.96 * predict(OUTPUT_DATA$m2, OUTPUT_DATA$predict_input_data, se.fit=TRUE)[[2]][[1]]))*100, digits=2),")", tags$br()))
  }
})
