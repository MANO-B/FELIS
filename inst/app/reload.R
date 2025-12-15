observeEvent(input$summary_base_reload, {
  clear_reactive_data(OUTPUT_DATA)
  gc()
  summary_base_logic()
})

observeEvent(input$CGP_benefit_analysis_reload, {
  clear_reactive_data(OUTPUT_DATA)
  gc()
  CGP_benefit_analysis_logic()
})

# 中間に反応制御用の reactive trigger を用意
oncoprint_trigger <- reactiveVal(0)
oncoprint_trigger_2 <- reactive(req((input$oncoprint_option)))
observeEvent(input$oncoprint_reload, {
  oncoprint_trigger(oncoprint_trigger() + 1)
})

# 中間に反応制御用の reactive trigger を用意
figure_survival_CGP_5_1_trigger <- reactiveVal(0)
figure_survival_CGP_5_1_trigger_2 <- reactive(req(input$figure_survival_CGP_5_1_max_samples))
observeEvent(input$figure_survival_CGP_5_1_reload, {
  figure_survival_CGP_5_1_trigger(figure_survival_CGP_5_1_trigger() + 1)
})

# 中間に反応制御用の reactive trigger を用意
figure_drug_trigger <- reactiveVal(0)
observeEvent(input$figure_survival_drug_reload, {
  figure_drug_trigger(figure_drug_trigger() + 1)
})

observeEvent(input$color_var_surv_CTx_reload, {
  survival_CTx_analysis2_logic()
})

observeEvent(input$figure_Bayes_reload, {
  survival_CTx_analysis_logic()
})

observeEvent(input$drug_table_reload, {
  drug_table_logic()
})

# figure_drug_load_trigger <- reactiveVal(0)

