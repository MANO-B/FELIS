library(shiny)
library(shinyWidgets)
library(shinydashboard)
library(shinyBS)
library(shinyjs)
library(shinyjqui)
library(ggplot2)
library(tidyr)
library(readr)
library(dplyr)
library(stringr)
library(RColorBrewer)
library(gt)
library(gtsummary)
library(flextable)
library(survival)
library(gridExtra)
library(survminer)
library(tranSurv)
library(DT)
library(ComplexHeatmap)
library(ggsci)
library(scales)
library(patchwork)
library(forcats)
library(PropCIs)
library(httr)
library(drawProteins)
library(ggrepel)
library(rms)
library(dcurves)
library(Matching)
library(survRM2)
library(pROC)
library(tidymodels)
library(ranger)
library(bonsai)
library(probably)
library(flexsurv)
library(data.table)
library(qs2)
library(mltools)
library(reshape2)
library(bigstep)
library(twang)
library(lightgbm)
library(ggbeeswarm)
library(ggiraph)
library(cmprsk)
library(zoo)
library(bayesplot)
library(tidybayes)
library(Rediscover)
library(rlang)
library(rstan)
library(matrixStats)
library(tools)
library(future)
library(furrr)
library(parallel)
library(aws.s3)
library(jsonlite)

APP_DIR <- system.file("app", package = "FELIS")
source_app <- function(fname) {
  f <- file.path(APP_DIR, fname)
  if (!file.exists(f)) stop("Missing file: ", f)
  source(f, local = parent.frame(), chdir = TRUE)
}

source_app("initialize.R")
source_app("RECIST.R")
source_app("survival_function.R")
source_app("ui.R")

# Define server logic required to draw a histogram
server <- function(input, output, session) {
  tmp_post <- reactiveValues()
  tmp_post$CCAT_FLAG = CCAT_FLAG
  analysis_env <- NULL

  OUTPUT_DATA = reactiveValues()
  source_app("initialize_each_run.R")
  source_app("case_load.R")
  source_app("drug_load.R")
  source_app("report_load.R")
  source_app("ui_helpers.R")
  source_app("cluster_load.R")
  source_app("reload.R")
  source_app("summary_base.R")
  source_app("oncoprint.R")
  source_app("mutually_exclusive.R")
  source_app("mut_subtype.R")
  source_app("custering_analysis.R")
  source_app("CGP_benefit_analysis.R")
  source_app("survival_simuration.R")
  source_app("survival_CGP_analysis.R")
  source_app("survival_CTx_analysis2.R")
  source_app("survival_CTx_analysis3.R")
  source_app("MCMC_function.R")
  source_app("survival_CTx_analysis.R")
  source_app("drug_table.R")
  source_app("volcano_function.R")
  source_app("drug_analysis.R")
  source_app("predict_nomogram.R")
  source_app("exit.R")
  source_app("memory.R")
  observeEvent(input$tabs, {
    if (input$tabs == "Settings") {
      analysis_env <- new.env()
      clear_reactive_data(OUTPUT_DATA)
      rm(analysis_env)
      gc()
    } else if (input$tabs == "Summarizedbymutationpattern" &&
        is.null(OUTPUT_DATA$table_summary_1_Data_summary)) {
      withProgress(message = "Data processing...", {
        summary_base_logic()
      })
    } else if (input$tabs %in% c("Oncoprint",
                          "Lolliplotfortheselectedgene",
                          "Tableofclinicalandmutationinformationperpatient") &&
        is.null(OUTPUT_DATA$oncoprint_Data_case_target_table)) {
      withProgress(message = "Data processing...", {
        oncoprint_logic()
      })
    } else if (input$tabs == "Mutualexclusivity" &&
               is.null(OUTPUT_DATA$figure_mutually_exclusive_Data_mut_ex)) {
      withProgress(message = "Data processing...", {
        mutually_exclusive_logic()
      })
    } else if (input$tabs == "Variantbyhistology" &&
               is.null(OUTPUT_DATA$figure_mut_subtype_plot_Data_tmp)) {
      withProgress(message = "Data processing...", {
        mut_subtype_logic()
      })
    } else if (input$tabs %in% c("Basicdata",
                                 "Frequencyofpatientswithtargetedtherapy",
                                 "Clusterandhistologyrelationship",
                                 "Heterogeneitywithinhistologictypes") &&
               is.null(OUTPUT_DATA$clustering_Data_mutation_cord)) {
      withProgress(message = "Data processing...", {
        custering_analysis_logic()
      })
    } else if (input$tabs %in% c("FactorsleadingtoTreatmentpre-CGPNomogram",
                                 "FactorsleadingtoTreatmentpre-CGPOddsratio",
                                 "FactorsleadingtoTreatmentdecisioncurve",
                                 "ROCnomogram",
                                 "Input_data") &&
               is.null(OUTPUT_DATA$table_prediction_Data_forest_tmp_8_GMT)) {
      withProgress(message = "Data processing...", {
        CGP_benefit_analysis_logic()
      })
    } else if (input$tabs %in% c("SurvivalafterCGP",
                                 "SurvivalandtreatmentafterCGP",
                                 "SurvivalafterCGPandmutationsforestplot",
                                 "HazardratioforsurvivalafterCGP_1",
                                 "Survival_treatment_reach") &&
               is.null(OUTPUT_DATA$figure_surv_CGP_Data_forest_tmp)) {
      withProgress(message = "Data processing...", {
        survival_CGP_analysis_logic()
      })
    } else if (input$tabs %in% c("SurvivalandtreatmentafterCTx",
                                 "SurvivalandtreatmentafterCTx2",
                                 "Geneticvariantsandsurvivalforestplot2",
                                 "Diagnosisandsurvivalforestplot2",
                                 "DiagnosisandsurvivalKM-curve2",
                                 "HazardratioforsurvivalafterCTx_1") &&
               is.null(OUTPUT_DATA$figure_surv_CTx_Data_survival_cancer)) {
      withProgress(message = "Data processing...", {
        survival_CTx_analysis2_logic()
      })
    } else if (input$tabs %in% c("ControlCustom"#,
                                 #"Simulation_Study"
                                 ) &&
               is.null(OUTPUT_DATA$figure_surv_CTx_Data_survival_cancer_control)) {
      withProgress(message = "Data processing...", {
        survival_CTx_analysis2_logic_control()
      })
    } else if (input$tabs %in% c("Survivalcorrectedforleft-truncationbias",
                                 "BayesCustom",
                                 "GeneticvariantsandsurvivalKM-curve",
                                 "DiagnosisandsurvivalKM-curve") &&
               is.null(OUTPUT_DATA$figure_surv_Bayes_Data_case_target)) {
      withProgress(message = "Data processing...", {
        survival_CTx_analysis_logic()
      })
    } else if (input$tabs %in% c("Drugusebylineoftreatment") &&
               is.null(OUTPUT_DATA$drug_table_Drug_summary)) {
      withProgress(message = "Data processing...", {
        drug_table_logic()
      })
    } else if (input$tabs %in% c("UseofdesignatedlinesanddrugswithToTinformationbymutationpattern",
                                 "Drugperpatient",
                                 "Timeontreatmentandpre-treatmentforthespecifiedtreatmentscatterplot",
                                 "ToT_interactive",
                                 "TimeontreatmentbygenemutationclusterKM-curve",
                                 "Timeontreatmentbymutatedgenesforestplot",
                                 "TimeontreatmentandmutationsofinterestKM-curve",
                                 "Hazardratioontimeontreatment_1",
                                 "Hazardratioontimeontreatment_2",
                                 "VolcanoPlot_ToT_1",
                                 "VolcanoPlotORR_1",
                                 "Oddsratioonobjectiveresponserate_1",
                                 "Oddsratioonobjectiveresponserate_2",
                                 "Oddsratioondiseasecontrolrate",
                                 "MutatedgenesandRECIST",
                                 "VolcanoPlotORR_1_AE",
                                 "Oddsratioonobjectiveresponserate_1_AE",
                                 "Oddsratioonobjectiveresponserate_2_AE",
                                 "Cumulative_incidence_AE",
                                 "Survival_drug"
                                 ) &&
               !is.null(OUTPUT_DATA$drug_table_Drug_summary) &&
               is.null(OUTPUT_DATA$Data_drug_Data_case_target)) {
      withProgress(message = "Data processing...", {
        drug_analysis_logic()
      })
    }
    if(input$tabs %in% c("VolcanoPlotORR_1_AE") &&
       is.null(OUTPUT_DATA$drug_analysis_g_volcano_AE) &&
       !is.null(OUTPUT_DATA$drug_analysis_Data_drug_RECIST) &&
       !is.null(input$target_line) &&
       !is.null(OUTPUT_DATA$Data_drug_Data_MAF_target_tmp) &&
       !is.null(OUTPUT_DATA$Data_drug_Data_case_target)) {
      volcano_AE_logic()
    }
    if(input$tabs %in% c("VolcanoPlotORR_1") &&
       is.null(OUTPUT_DATA$drug_analysis_g_volcano) &&
       !is.null(OUTPUT_DATA$drug_analysis_Data_drug_RECIST) &&
       !is.null(input$target_line) &&
       !is.null(OUTPUT_DATA$Data_drug_Data_MAF_target_tmp) &&
       !is.null(OUTPUT_DATA$Data_drug_Data_case_target)) {
      volcano_ORR_logic()
    }
    if(input$tabs %in% c("VolcanoPlot_ToT_1") &&
       is.null(OUTPUT_DATA$drug_analysis_g_volcano_ToT) &&
       !is.null(OUTPUT_DATA$Data_drug_Data_MAF_target_tmp) &&
       !is.null(OUTPUT_DATA$Data_drug_Data_drug_TTF)) {
      volcano_ToT_logic()
    }
  })
}

