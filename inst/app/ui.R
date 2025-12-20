# Define UI for application that draws a histogram
ui <- dashboardPage(
  skin = "black",
  title = paste("FELIS version ", SOFT_VERSION, ", data version ", DATA_VERSION),
  dashboardHeader(
    title = HTML(
      paste0(
        "<div style='display: flex; align-items: center; background-color: #FFFFFF;'>
           <a href='https://github.com/MANO-B/FELIS'>
             <img src='FELIS.png' height='30px' style='margin-right: 10px;'>
           </a>
           <h6 style='margin: 0;'>Data v", DATA_VERSION, ", ", ENV_, "</h6>
         </div>"
      )
    ),
    tags$li(class = "dropdown",
            style = "height: 30px; padding: 10px;",
            div(id = "memory-widget",
                style = "height: 28px; line-height: 28px; color: black; padding: 0px 10px;
                       padding: 0px 10px;
                       display: flex; align-items: center;",
                span(id = "memory-text", "Loading...")
            )
    ),
    tags$li(class = "dropdown",
            style = "height: 28px; padding: 0px;",
            actionButton(
              inputId = "logout_button",
              label = "Logout",
              icon = icon("sign-out"),
              style = "height: 28px; line-height: 28px; color: white; background-color: #d9534f; border-color: #ac2925;"
            )
    )
  ),
  dashboardSidebar(
    sidebarMenu(id = "tabs",
      h3("Settings"),
      menuItem("Settings", tabName = "Setting", icon = icon("dashboard")),
      hr(),
      h3("Results"),
      menuItem("Case summary", tabName = "Summarizedbymutationpattern", icon = icon("th")),
      menuItem("Oncoprint", icon = icon("th"),
               hr(),
               p("Figures"),
               menuSubItem("Oncoprint", tabName = "Oncoprint", icon = icon("angle-right")),
               menuSubItem("Lolliplot for the selected gene", tabName = "Lolliplotfortheselectedgene", icon = icon("angle-right")),
               p("Downloadable table"),
               menuSubItem("Table of clinical and mutation information per patient", tabName = "Tableofclinicalandmutationinformationperpatient", icon = icon("angle-right")),
               hr()
      ),
      menuItem("Mutual exclusivity", tabName = "Mutualexclusivity", icon = icon("th")),
      menuItem("Variat rate by histology", tabName = "Variationbyhistology", icon = icon("th")),
      menuItem("Mutation and treatment option", icon = icon("th"),
               hr(),
               p("Mutation-based clustering"),
               menuSubItem("Basic data", tabName = "Basicdata", icon = icon("angle-right")),
               menuSubItem("Frequency of patients with targeted therapy", tabName = "Frequencyofpatientswithtargetedtherapy", icon = icon("angle-right")),
               menuSubItem("UMAP clustering based on mutations", tabName = "Clusterandhistologyrelationship", icon = icon("angle-right")),
               menuSubItem("Heterogeneity within histologic types", tabName = "Heterogeneitywithinhistologictypes", icon = icon("angle-right"))
      ),
      menuItem("Survival after CGP", icon = icon("th"),
               hr(),
               p("Survival analysis"),
               menuSubItem("Survival and clinical information", tabName = "SurvivalafterCGP", icon = icon("angle-right")),
               menuSubItem("Custom survival analysis", tabName = "SurvivalandtreatmentafterCGP", icon = icon("angle-right")),
               menuSubItem("Survival and mutations, forest plot", tabName = "SurvivalafterCGPandmutationsforestplot", icon = icon("angle-right")),
               menuSubItem("Hazard ratio", tabName = "HazardratioforsurvivalafterCGP_1", icon = icon("angle-right")),
               menuSubItem("Survival period and treatment reach rate", tabName = "Survival_treatment_reach", icon = icon("angle-right")),
               hr()
      ),
      menuItem("CGP benefit prediction", icon = icon("th"),
               hr(),
               p("Factors leading to treatment"),
               menuSubItem("Nomogram", tabName = "FactorsleadingtoTreatmentpre-CGPNomogram", icon = icon("angle-right")),
               menuSubItem("Odds ratio", tabName = "FactorsleadingtoTreatmentpre-CGPOddsratio", icon = icon("angle-right")),
               menuSubItem("Decision curve", tabName = "FactorsleadingtoTreatmentdecisioncurve", icon = icon("angle-right")),
               menuSubItem("ROC curve of nomogram", tabName = "ROCnomogram", icon = icon("angle-right")),
               hr(),
               p("Prediction for your data"),
               menuSubItem("Input data", tabName = "Input_data", icon = icon("th")),
               hr()
      ),
      menuItem("Overall survival with risk-set adjustment", icon = icon("th"),
               hr(),
               p("Survival analysis"),
               menuSubItem("Survival and clinical information", tabName = "SurvivalandtreatmentafterCTx", icon = icon("angle-right")),
               menuSubItem("Custom survival analysis", tabName = "SurvivalandtreatmentafterCTx2", icon = icon("angle-right")),
               menuSubItem("Frequent variants and survival", tabName = "Geneticvariantsandsurvivalforestplot2", icon = icon("angle-right")),
               menuSubItem("Diagnosis and survival", tabName = "Diagnosisandsurvivalforestplot2", icon = icon("angle-right")),
               menuSubItem("Mutational cluster and survival", tabName = "DiagnosisandsurvivalKM-curve2", icon = icon("angle-right")),
               menuSubItem("Hazard ratio for survival after CTx", tabName = "HazardratioforsurvivalafterCTx_1", icon = icon("angle-right")),
               hr()
      ),
      menuItem("Survival after CTx with Bayesian inference", icon = icon("th"),
               hr(),
               menuSubItem("Survival corrected for left-truncation bias", tabName = "Survivalcorrectedforleft-truncationbias", icon = icon("angle-right")),
               menuSubItem("Custom survival analysis", tabName = "BayesCustom", icon = icon("angle-right")),
               menuSubItem("Genetic variants and survival", tabName = "GeneticvariantsandsurvivalKM-curve", icon = icon("angle-right")),
               menuSubItem("Diagnosis and survival", tabName = "DiagnosisandsurvivalKM-curve", icon = icon("angle-right")),
               hr()
      ),
      menuItem("Bias correction simulation", tabName = "SurvivalSimurationKMCurve", icon = icon("th")),
      menuItem("Drug response", icon = icon("th"),
               hr(),
               p("Settings for following analyses"),
               menuSubItem("Settings", tabName = "Drugusebylineoftreatment", icon = icon("angle-right")),
               hr(),
               p("Summary Tables"),
               menuSubItem("Selected drugs", tabName = "UseofdesignatedlinesanddrugswithToTinformationbymutationpattern", icon = icon("angle-right")),
               menuSubItem("Drug usage data", tabName = "Drugperpatient", icon = icon("angle-right")),
               hr(),
               p("Treatment time"),
               menuSubItem("Treatment time and clinical information", tabName = "Timeontreatmentandpre-treatmentforthespecifiedtreatmentscatterplot", icon = icon("angle-right")),
               menuSubItem("Treatment time comparison", tabName = "ToT_interactive", icon = icon("angle-right")),
               menuSubItem("Treatment time by gene mutation cluster", tabName = "TimeontreatmentbygenemutationclusterKM-curve", icon = icon("angle-right")),
               menuSubItem("Treatment time by mutated genes", tabName = "Timeontreatmentbymutatedgenesforestplot", icon = icon("angle-right")),
               menuSubItem("Treatment time on treatment and mutations of interest", tabName = "TimeontreatmentandmutationsofinterestKM-curve", icon = icon("angle-right")),
               menuSubItem("Hazard ratio on time on treatment - genes", tabName = "Hazardratioontimeontreatment_1", icon = icon("angle-right")),
               menuSubItem("Hazard ratio on time on treatment - clusters", tabName = "Hazardratioontimeontreatment_2", icon = icon("angle-right")),
               menuSubItem("Volcano plot for treatment time, hazard ratio", tabName = "VolcanoPlot_ToT_1", icon = icon("angle-right")),
               hr(),
               p("Response rate"),
               menuSubItem("Volcano plot for objective response rate", tabName = "VolcanoPlotORR_1", icon = icon("angle-right")),
               menuSubItem("Odds ratio on objective response rate - genes", tabName = "Oddsratioonobjectiveresponserate_1", icon = icon("angle-right")),
               menuSubItem("Odds ratio on objective response rate - clusters", tabName = "Oddsratioonobjectiveresponserate_2", icon = icon("angle-right")),
               menuSubItem("Odds ratio on disease control rate", tabName = "Oddsratioondiseasecontrolrate", icon = icon("angle-right")),
               menuSubItem("Factors for response rate", tabName = "MutatedgenesandRECIST", icon = icon("angle-right")),
               hr(),
               p("Adverse effect"),
               menuSubItem("Volcano plot for adverse effect", tabName = "VolcanoPlotORR_1_AE", icon = icon("angle-right")),
               menuSubItem("Odds ratio on adverse effect - genes", tabName = "Oddsratioonobjectiveresponserate_1_AE", icon = icon("angle-right")),
               menuSubItem("Odds ratio on adverse effect - clusters", tabName = "Oddsratioonobjectiveresponserate_2_AE", icon = icon("angle-right")),
               menuSubItem("Cumulative incidence of adverse effect", tabName = "Cumulative_incidence_AE", icon = icon("angle-right")),
               hr(),
               p("Survival after drug initiation date"),
               menuSubItem("Survival and drug", tabName = "Survival_drug", icon = icon("angle-right"))
      ),
      hr(),
      h3("Instruction"),
      menuItem("About FELIS", tabName = "Instruction", icon = icon("th")),
      menuItem("Tips", tabName = "Tips", icon = icon("th")),
      hr()
    )
  ),
  ## Body content
  dashboardBody(
    # hr(),
    useShinyjs(),
    tags$head(
      tags$style(HTML("
        /* 既存のCSS */
        .main-header {
          height: 30px;
          line-height: 30px;
          z-index: 1000;
        }
        .main-header .logo {
          height: 30px;
          line-height: 30px;
        }
        .navbar {
          min-height: 30px;
        }

        #logout_button {
          transition: background-color 0.3s;
          margin: 0;
          padding: 2px 8px;
          height: 20px;
          font-size: 12px;
          margin-bottom: 0 !important;
        }
        #memory-widget {
          margin: 1;
          padding: 2px 8px;
          height: 22px;
          font-size: 12px;
          margin-bottom: 0 !important;
          pointer-events: none !important;
          user-select: none !important;
        }

        /* レスポンシブ対応 */
        @media (max-width: 768px) {
          #memory-widget {
            font-size: 10px !important;
            padding: 6px 8px !important;
          }
        }

        .sidebar-toggle {
          height: 30px;
          padding-top: 1px !important;
          margin-bottom: 0 !important;
          padding-bottom: 0 !important;
          margin-top: 0 !important;
        }
        .sidebar-menu li a {
          white-space: normal;
          word-wrap: break-word;
          overflow-wrap: break-word;
          line-height: 1.5;
        }
        .sidebar-menu {
          word-wrap: break-word;
        }
        .main-sidebar {
          position: fixed;
          margin-top: 0px;
          left: 0;
          height: calc(100vh - 20px);
          width: 230px;
          overflow-y: auto;
          overflow-x: hidden;
          scrollbar-width: none;
          scroll-behavior: smooth;
        }
        .sidebar {
          height: 100%;
          overflow-y: auto;
          overflow-x: hidden;
        }
        .sidebar-menu {
          margin-top: 0px;
        }
        .content-wrapper, .right-side, .content {
          position: fixed;
          top: 0px;
          left: 230px;
          margin-top: 30px;
          padding-top: 0px;
          overflow-y: auto;
          height: calc(100vh - 0px);
          width: calc(100% - 230px);
          background-color: #ebf1fc;
        }
        .sidebar-collapse .content-wrapper,
        .sidebar-collapse .right-side,
        .sidebar-collapse .content {
          left: 0px !important;
          width: calc(100% - 0px) !important;
        }

        .sidebar-collapse .main-sidebar {
          width: 0px !important;
        }

        .content-wrapper:hover::-webkit-scrollbar {
          width: 8px;
        }
        .content-wrapper:hover::-webkit-scrollbar-thumb {
          background-color: #888;
          border-radius: 4px;
          border: 2px solid white;
        }
        .main-header .navbar .nav > li {
          margin: 0 !important;
          padding: 0 !important;
        }
        .main-header .sidebar-toggle:before {
          content: '\\e068';
        }
        .main-header .navbar-nav {
          margin: 0 !important;
        }
        .main-header .navbar-nav > li {
          margin: 0 !important;
          padding: 0 !important;
        }
      "))
    ),
    # JavaScript
    tags$script(src = "felis-src/jszip.min.js"),
    tags$script(HTML("
    Shiny.addCustomMessageHandler('updateMemoryWidget', function(data) {
      $('#memory-text').text(data.text);
    });
    function confirmDownload(elementId, fileLabel) {
      var message = 'Please manage according to the personal information handling rules\\n=================================\\nPlease comply with the contents of the Agreement on Utilization of C-CAT Data and the service specification conformity disclosure of/for C-CAT Research-Use Portal site, and handle it properly.\\n=================================\\nBe careful not to leave the downloaded file on your computer';

      if (confirm(message)) {

        // [追加] 通知を /@@/api/writwlog に送る
        fetch(window.location.origin + '/@@/api/writwlog', {
          method: 'POST',
          headers: { 'Content-Type': 'application/json' },
          body: JSON.stringify({
            // メッセージ内容が 'csv file' のため 'csv' を指定
            type: 'csv',
            time: new Date().toISOString(),
            file: fileLabel
          })
        }).catch(err => console.error('Download log fetch error:', err));

        // 本来のクリック（ダウンロード）動作
        document.getElementById(elementId).click();
      }
    }
    window.addEventListener('unload', function (e) {
     // リロード中は confirm が表示されるので、キャンセルされた場合は GETしない
     if (confirm('Are you closing the window? Click Cancel if you are reloading.')) {
       Shiny.setInputValue('browser_closed', true, {priority: 'event'});
     }
    });
    ")),
    tabItems(
      tabItem(tabName = "Setting",
              fluidRow(
                column(6,
                       fluidRow(box(
                         title = "How to use",
                         status = "primary",
                         solidHeader = TRUE,
                         collapsible = TRUE,
                         collapsed = TRUE,
                         width = 6,

                         tags$h5("解析対象の症例・遺伝子を選択します"),
                         tags$ul(
                           tags$li("Analyze with new dataset をNoにすると全症例または前回解析症例が取り込まれます"),
                           tags$li("Yesにすると、アップロードしたCSVファイルが取り込まれます"),
                           tags$li("必要があれば設定をデフォルトから変更します"),
                           tags$li("フィルタリング後にResultsのタブを選択し、解析結果を閲覧します")
                         )
                       )),
                       br(),
                       conditionalPanel(
                         condition = "input.new_analysis == 'Yes, input new csv files'",
                          div(
                           strong("Required: C-CAT data files (zipped files are acceptable)"),
                           actionButton("help_ccat", "",
                                        icon = icon("question-circle"),
                                        class = "help-icon",
                                        style = "border: none; background: transparent;"),
                           bsPopover("help_ccat",
                                     title = "C-CAT データファイルについて",
                                     content = "C-CAT（がんゲノム情報管理センター）から提供されるデータファイルをアップロードしてください。ZIPファイル形式でも対応しています。",
                                     placement = "top",
                                     trigger = "hover")
                         ),
                         fileInput(
                           inputId = "clinical_files",
                           label = div("Choose case CSV Files",
                                       actionButton("help_clinical_files", "",
                                                    icon = icon("question-circle"),
                                                    class = "help-icon",
                                                    style = "border: none; background: transparent;")
                           ),
                           multiple = TRUE
                         ),
                         bsPopover("help_clinical_files",
                                   title = "症例CSVファイル",
                                   content = "患者の臨床情報が含まれるCSVファイルを選択してください。複数ファイルの選択が可能です。ファイル形式：CSV、文字エンコーディング：UTF-8推奨",
                                   placement = "top",
                                   trigger = "hover"),
                         fileInput(
                           inputId = "report_files",
                           label = div("Choose report CSV Files",
                                       actionButton("help_report_files", "",
                                                    icon = icon("question-circle"),
                                                    class = "help-icon",
                                                    style = "border: none; background: transparent;")
                           ),
                           multiple = TRUE
                         ),
                         bsPopover("help_report_files",
                                   title = "レポートCSVファイル",
                                   content = "解析レポートデータが含まれるCSVファイルを選択してください。症例データと対応するレポートファイルをアップロードしてください。",
                                   placement = "top",
                                   trigger = "hover"),
                         br(),
                         br(),
                         br(),
                         hr()
                       )
                )
              ),
              fluidRow(
                column(3,strong("Filter on histology"),
                       hr(),
                       htmlOutputWithPopover(
                         "select_organ",
                         "解析する腫瘍の部位",
                         "OncoTree 1st levelで選択してください"
                       ),
                       htmlOutputWithPopover(
                         "select_histology",
                         "解析する腫瘍の詳細部位",
                         "詳細なOncoTreeの組織型を選択してください"
                       ),
                       htmlOutputWithPopover(
                         "select_minimum_pts",
                         "解析する詳細な組織型の最小患者数",
                         "この数値未満の組織型はOncoTree 1st levelに変換されます"
                       ),
                       htmlOutputWithPopover(
                         "select_histology_detail",
                         "詳細な組織型が不要な場合",
                         "Yesの場合OncoTree level 1で解析"
                       ),
                       br(),
                       hr(),
                       htmlOutputWithPopover(
                         "select_new_analysis",
                         "既存のデータか新規データかを選択して解析",
                         "Yesの場合、ファイルをアップロードして解析します"
                       ),
                       br()
                ),
                column(3,strong("Filters for clinical information"),
                       hr(),
                       htmlOutput("select_sex"),
                       br(),
                       htmlOutput("select_panel"),
                       htmlOutput("select_age"),
                       htmlOutputWithPopover(
                         "select_mid_age",
                         "年齢の閾値",
                         "指定年齢以下とそれより高齢で2群に分けます"
                       ),
                       htmlOutput(""),
                       htmlOutput("select_PS"),
                       br(),
                       htmlOutput("select_smoking"),
                       br(),
                       htmlOutput("select_year"),
                       htmlOutputWithPopover(
                         "select_year",
                         "CGP検査の実施年",
                         "2019年と直近は欠損値が多いです"
                       ),
                       br()
                ),
                column(3,strong("Filters on genes"),
                       hr(),
                       htmlOutputWithPopover(
                         "select_TMB_threshold",
                         "TMB highとする閾値",
                         "TMBとICIの効果の関連の解析等に使用"
                       ),
                       htmlOutputWithPopover(
                         "select_gene",
                         "注目する遺伝子セット",
                         "優先的に自動解析対象となる遺伝子"
                       ),
                       htmlOutputWithPopover(
                         "select_gene_group_1",
                         "上で選択した遺伝子をさらに群分け",
                         "変異パターンごとの解析に使用"
                       ),
                       htmlOutput("select_gene_group_2"),
                       htmlOutput("select_gene_group_analysis"),
                       br(),
                       br(),
                       strong("For detailed study of mutations of a gene"),
                       htmlOutputWithPopover(
                         "select_special_gene",
                         "変異体を詳細に解析したい遺伝子を選択",
                         "Hotspotやドメインなどに注目した解析が可能"
                       ),
                       htmlOutputWithPopover(
                         "select_special_gene_mutation_1",
                         "まとめて解析したい変異を選択",
                         "下で名称を指定"
                       ),
                       htmlOutput("select_special_gene_mutation_1_name"),
                       htmlOutput("select_special_gene_mutation_2"),
                       htmlOutput("select_special_gene_mutation_2_name"),
                       htmlOutput("select_special_gene_annotation"),
                       htmlOutput("select_special_gene_independent"),
                       br()
                ),
                column(3,strong("Other Settings"),
                       hr(),
                       htmlOutputWithPopover(
                         "select_gene_no",
                         "自動解析する頻度の高い遺伝子数を指定",
                         "生存期間解析等で使用"
                       ),
                       htmlOutput(""),
                       htmlOutputWithPopover(
                         "select_patho",
                         "解析対象とする病的意義を指定",
                         "CKDB evidence level Fをoncogenicと定義"
                       ),
                       htmlOutput("select_T_N"),
                       htmlOutput("select_fusion"),
                       htmlOutput("select_amplification"),
                       htmlOutputWithPopover(
                         "select_clustering",
                         "遺伝子変異に基づくクラスタリング解析の設定",
                         "C-CAT全症例での処理済みデータを使用可能。選択した症例で再計算は可能ですが時間がかかります"
                       ),
                       htmlOutputWithPopover(
                         "select_evidence_level",
                         "CGP benefit analysisにおける治療到達に関係すると考えるエビデンスレベル",
                         "A&BあるいはA&B&Cが通常と思います"
                       )

                ),
              ),
              hr(),
              h5("Figures in results are downloadable as png files."),
              h6("FELIS; Functions Especially for LIquid and Solid tumor clinical sequencing."),
              a(href="https://github.com/MANO-B/FELIS",h6("https://github.com/MANO-B/FELIS")),
              hr(),
              br(),
              br(),
              h4(""),
              fluidRow(box(
                title = "The following settings are for advanced analysis only",
                status = "primary",
                solidHeader = TRUE,
                collapsible = TRUE,
                collapsed = TRUE,
                width = 6,

                tags$h5(""),
                tags$ul(
                  tags$li("動作はしますが、とくに変更は不要です")
                )
              )),
              hr(),
              br(),
              fluidRow(
                column(3,
                       br(),
                       htmlOutputWithPopover(
                         "select_histology_group_1",
                         "組織型をまとめて一つの腫瘍としたい場合",
                         "まとめたい組織型を選択し、名称をつけます"
                       ),
                       htmlOutput("select_histology_group_1_name"),
                       htmlOutput("select_histology_group_2"),
                       htmlOutput("select_histology_group_2_name"),
                       htmlOutput("select_histology_group_3"),
                       htmlOutput("select_histology_group_3_name"),
                       htmlOutput("select_histology_group_4"),
                       htmlOutput("select_histology_group_4_name"),
                       br(),
                       downloadButton('download_test_clinical_data', 'Download sample clinical data'),
                       br(),
                       br(),
                       downloadButton('download_test_report_data', 'Download sample report data'),
                       br()
                ),
                column(3,
                       strong("Option files"),
                       hr(),
                       fileInput(inputId = "ID_histology",
                                 label = "Correspondence table between ID and histology (CSV)",
                                 multiple = TRUE,
                                 accept = c(
                                   "text/csv",
                                   "text/comma-separated-values,text/plain",
                                   ".csv")
                       ),
                       "Specify the ID of the patient whose diagnosis you want to correct and the modified histology.",
                       br(),
                       br(),
                       downloadButton('download_ID_histology', 'Download CSV file template'),
                       br(),
                       br(),
                       hr(),
                       fileInput(inputId = "ID_drug",
                                 label = "Correspondence table between ID and drug information (CSV)",
                                 multiple = TRUE,
                                 accept = c(
                                   "text/csv",
                                   "text/comma-separated-values,text/plain",
                                   ".csv")
                       ),
                       "To analyze only the drug after curation",
                       br(),
                       br(),
                       downloadButton('download_ID_drug_pre_CGP', 'Download CSV template (before CGP)'),
                       br(),
                       br(),
                       downloadButton('download_ID_drug_post_CGP', 'Download CSV template (after CGP)'),
                       br(),
                       br(),
                       hr(),
                       fileInput(inputId = "drug_rename",
                                 label = "Correspondence table for drug renaming (CSV)",
                                 multiple = TRUE,
                                 accept = c(
                                   "text/csv",
                                   "text/comma-separated-values,text/plain",
                                   ".csv")
                       ),
                       "To analyze similar drugs together by renaming the drugs",
                       br(),
                       "Reclassified as molecular targeted therapies, immune checkpoint inhibitors, etc.",
                       br(),
                       "Drugs are listed in ABC order, separated by commas",
                       br(),
                       br(),
                       downloadButton('download_drug_rename', 'Download CSV template'),
                       br(),
                       br(),
                       fileInput(inputId = "drug_combination_rename",
                                 label = "Correspondence table for drug combination renaming (CSV)",
                                 multiple = TRUE,
                                 accept = c(
                                   "text/csv",
                                   "text/comma-separated-values,text/plain",
                                   ".csv")
                       ),
                       "To rename drug combinations to groups",
                       br(),
                       "Reclassified as molecular targeted therapies, immune checkpoint inhibitors, etc.",
                       br(),
                       "Drugs are listed in ABC order, separated by commas",
                       br(),
                       br(),
                       downloadButton('download_drug_combination_rename', 'Download CSV template'),
                       br(),
                       br()
                ),
                column(3,
                       strong("Option files"),
                       hr(),
                       fileInput(inputId = "mutation_rename",
                                 label = "Correspondence table for mutation renaming (CSV)",
                                 multiple = TRUE,
                                 accept = c(
                                   "text/csv",
                                   "text/comma-separated-values,text/plain",
                                   ".csv")
                       ),
                       "To analyze similar mutations together by renaming mutations",
                       br(),
                       "Reclassified as Exon 19 mutation, Exon 20 mutation, gene amplification, etc.",
                       br(),
                       "'Other' if there is an unspecified mutation in the designated gene",
                       br(),
                       br(),
                       downloadButton('download_mutation_rename', 'Download CSV template'),
                       br(),
                       br(),
                       hr(),
                       fileInput(inputId = "reaanotation",
                                 label = "Correspondence table for mutation reannotation (CSV)",
                                 multiple = TRUE,
                                 accept = c(
                                   "text/csv",
                                   "text/comma-separated-values,text/plain",
                                   ".csv")
                       ),
                       "To reannotate variants",
                       br(),
                       "F: pathogenic variants, G: neutral variants",
                       br(),
                       br(),
                       downloadButton('download_reannotation', 'Download CSV template'),
                       br(),
                       br(),
                       hr(),
                       fileInput(inputId = "histology_rename",
                                 label = "Correspondence table of histological type renaming (CSV)",
                                 multiple = TRUE,
                                 accept = c(
                                   "text/csv",
                                   "text/comma-separated-values,text/plain",
                                   ".csv")
                       ),
                       "To analyze similar tissue types together by renaming them",
                       br(),
                       "Reclassified as differentiated gastric cancer, undifferentiated gastric cancer, etc.",
                       br(),
                       br(),
                       downloadButton('download_histology_rename', 'Download CSV template'),
                       br(),
                       br(),
                       hr(),
                       fileInput(inputId = "regimen_rename",
                                 label = "Correspondence table for drug renaming based on regimen (CSV)",
                                 multiple = TRUE,
                                 accept = c(
                                   "text/csv",
                                   "text/comma-separated-values,text/plain",
                                   ".csv")
                       ),
                       "Provide the drug name if the drug name is unknown and the regimen is known.",
                       br(),
                       "Conversion from MAP therapy to 'Cisplatin,Doxorubicin,Methotrexate'",
                       br(),
                       br(),
                       downloadButton('download_regimen_rename', 'Download CSV template'),
                       br(),
                       br()
                ),
                column(3,
                       strong("Analysis setting"),
                       hr(),
                       htmlOutputWithPopover(
                         "select_eps",
                         "UMAPクラスタリングの設定",
                         "EPSが大きい場合、再計算時のクラスタ数が減少します"
                       ),
                       htmlOutputWithPopover(
                         "select_RMST_CGP",
                         "CGP検査後の生存期間解析の設定",
                         "Restricted mean survival timeの測定日を設定"
                       ),
                       htmlOutputWithPopover(
                         "select_RMST_CTx",
                         "Risk-set adjustmentでの全生存期間解析の設定",
                         "Restricted mean survival timeの測定日を設定"
                       ),
                       htmlOutputWithPopover(
                         "select_RMST_drug",
                         "投薬開始後の治療継続期間解析の設定",
                         "Restricted mean survival timeの測定日を設定"
                       ),
                       htmlOutput("select_MSI"),
                       htmlOutput("select_MMR"),
                       htmlOutput("select_HER2"),
                       htmlOutput("select_PD_L1"),
                       htmlOutput("select_AE_analysis"),
                       htmlOutput("select_machine_learning"),
                       htmlOutput("select_importance"),
                       htmlOutput("select_drug_volcano"),
                       htmlOutput("select_volcano_pathology"),
                       # htmlOutputWithPopover(
                       #   "select_URL",
                       #   "同意撤回されたIDを記載したリストのURL",
                       #   "1列のcsvファイルとしてダウンロード"
                       # ),
                       h6("If you select No, csv files may not be necessary."),
                       htmlOutput("select_intermediate_file"),
                       h6("If you select Yes, faster when performing the same analysis repeatedly.")
                )
              )
      ),
      tabItem(tabName = "Summarizedbymutationpattern",
              fluidRow(
                column(4,
                       htmlOutput("select_table_summary"),
                       htmlOutput("select_table_summary_1_gene_group_1"),
                       htmlOutput("select_table_summary_1_gene_group_2")),
                column(4,
                       htmlOutput("select_table_summary_no"),
                       htmlOutput("select_table_summary_layout")
                ),
                column(4,
                       strong("Click button after setting modification"),
                       br(),
                       actionButton("summary_base_reload", "Reload")
                )
              ),
              hr(),
              gt_output("table_summary_1"),
              hr(),
              h6(paste("Table. Characteristics for selected patients.",
                       "The present, retrospective cohort study was performed with clinicogenomic,",
                       "real-world data on the patients who were registered in the C-CAT database",
                       "from June 1, 2019. The patients were registered by hospitals throughout Japan",
                       "and provided written informed consent to the secondary use of their clinicogenomic",
                       "data for research."))
      ),
      tabItem(tabName = "Oncoprint",
              h6("Figure. Recurrent oncogenic mutations in selected cases. The 30 genes with the highest frequency of oncogenic mutations are shown. Mutational landscapes were created using ComplexHeatmap package for R."),
              fluidRow(
                column(4,
                       htmlOutput("select_oncoprint_option")),
                column(4,
                       htmlOutputWithPopover(
                         "select_oncoprint_gene_no",
                         "表示遺伝子数",
                         "下部で遺伝子の並び順を指定してください"
                       ),
                       htmlOutputWithPopover(
                         "select_oncoprint_max_samples",
                         "表示対象症例数",
                         "症例数が多い場合、表示の高速化のためサンプリングして描画します"
                       ),
                       htmlOutput("select_oncoprint_gene")
                ),
                column(4,
                       strong("Click button after setting modification"),
                       br(),
                       actionButton("oncoprint_reload", "Reload")
                )
              ),
              hr(),
              plotOutput('figure_oncoprint',
                         height = "1000px",
                         width = "1000px"),
              hr(),
              wellPanel(
                uiOutput("dynamic_order")
              )
      ),
      tabItem(tabName = "Lolliplotfortheselectedgene",
              htmlOutput("select_lolliplot_gene"),
              htmlOutput("select_lolliplot_no"),
              hr(),
              girafeOutput('figure_lolliplot2',
                           height = "500px",
                           width = "1500px"),
              hr(),
              h6("Figure. Frequency of oncogenic mutations in the selected gene. The most frequent oncogenic mutations are shown with amino acid change."),
              h6("Mutplot by Zhang W, PMID:31091262. If error occurs, correct 'source/UniPlot.txt'."),
              h6("Protein structure source: Uniprot"),
              HTML("<p>Github for Mutplot. <a href='https://github.com/VivianBailey/Mutplot'>Link for the website</a></p>")
      ),
      tabItem(tabName = "Tableofclinicalandmutationinformationperpatient",
              DT::dataTableOutput("table_patient"),
              hr(),
              actionButton('show_download_confirm', 'Download all mutation data',
                           icon = icon("download"),
                           onclick = "confirmDownload('download_mutation_data', 'FELIS downloaded raw data in Table of clinical and mutation information per patient tab')",
                           class = "btn-primary"),

              # 非表示のダウンロードボタン
              div(style = "visibility:hidden; position:absolute;",
                  downloadButton('download_mutation_data', 'Hidden Download')
              )
      ),
      tabItem(tabName = "Mutualexclusivity",
              strong("Figure for probability"),
              girafeOutput('figure_mutually_exclusive_1',
                           height = "1500px",
                           width = "1500px"),
              hr(),
              strong("Figure for odds ratio"),
              girafeOutput('figure_mutually_exclusive_2',
                           height = "1500px",
                           width = "1500px"),
              hr(),
              DT::dataTableOutput("table_mutually_exclusive"),
              h6(paste0("Figure. Alterations among mutually exclusive or co-occurring pairs. The 30 genes with the highest frequency of oncogenic mutations were selected to determine ",
                        "whether oncogenic mutations are likely to occur simultaneously between the two genes. Blue boxes indicates mutually exclusivity and red boxes indicates co-occurrence. ",
                        "An asterisk shows a significant correlation (p < 0.001). Analysis was performed with Rediscover package in R language. Odds ratios were estimated by Fisher exact test. ",
                        "An odds ratio less than 1 does not necessarily correspond to mutual exclusivity as evaluated by the negative binomial distribution. ",
                        "This statement highlights that a low odds ratio (OR < 1) indicates a negative association between two events but does not inherently imply mutual exclusivity. ",
                        "In statistical modeling, particularly with count data, the negative binomial distribution is often employed to account for overdispersion. However, the interpretation of ",
                        "mutual exclusivity requires careful consideration beyond the OR value alone. ",
                        "For instance, in the context of count data, the negative binomial regression model can be used to estimate the odds of an event occurring. ",
                        "However, the OR derived from such models may not fully capture the complexity of mutual exclusivity between events. Factors such as overdispersion and the underlying data ",
                        "distribution can influence the interpretation of the OR. ",
                        "Therefore, while an OR less than 1 suggests a negative association, it should not be solely relied upon to infer mutual exclusivity, especially when using models like the",
                        "negative binomial distribution. A comprehensive analysis considering the specific context and model assumptions is essential for accurate interpretation."))
      ),
      tabItem(tabName = "Variationbyhistology",
              h6("Figure. Recurrent oncogenic mutations across subtypes. The 30 genes with the highest frequency of oncogenic mutations were displayed."),
              hr(),
              girafeOutput('figure_mut_subtype_plot')
      ),
      tabItem("Basicdata",
              plotOutput('figure_base',
                         height = "1500px",
                         width = "1000px"),
              hr(),
              h6("Figure. Distribution of age, sex, detected oncogenic mutations, tumor mutation burden (TMB), metastasis pattern, patients with treatment option recommended by the expert panel, patients received recommended chemotherapy, mediantime from the initiation date of the first palliative chemotherapy to CGP, and median time from CGP to final observation. In the boxplots of age and TMB, the box borders indicate the 25th and 75th percentiles, the inner line the median, and the whiskers 1.5× the interquartile range."),
              DT::dataTableOutput("table_basic_data")
      ),
      tabItem("Clusterandhistologyrelationship",
              htmlOutput("select_color_var_cluster"),
              girafeOutput('figure_cluster_subtype',
                           height = "800px",
                           width = "1600px"),
              hr(),
              h4("Summary, cluster and mutated gene"),
              DT::dataTableOutput("table_mutation"),
              hr(),
              h4("Summary, cluster and histology"),
              DT::dataTableOutput("table_disease"),
              h4("Raw data"),
              DT::dataTableOutput("table_basic_data_raw"),
              h6("Figure. Unsupervised clustering of the patients based on the detected oncogenic mutations. Two-dimensional mutational pattern mapping was generated using Uniform Manifold Approximation and Projection (UMAP). The three variants and histotypes with the highest odds ratios that were more common than the other clusters at p<0.05. Clustering analysis was performed as follows. All pathogenic mutations detected by the cancer-related genes were assembled into a binary matrix format per patient. The dimension of this input matrix was reduced using Uniform Manifold Approximation and Projection (UMAP) via the umap package for R (with default hyperparameters). Clustering analysis was performed using the Hierarchical Density-Based Spatial Clustering of Applications with Noise (HDBSCAN) via the dbscan package for R (EPS: 1.0; minimum points: 3)."),
              HTML("<p>Mochizuki T, et al., Factors predictive of second-line chemotherapy in soft tissue sarcoma: An analysis of the National Genomic Profiling Database. Cancer Science, 2023. <a href='https://doi.org/10.1111/cas.16050'>Link for the paper</a></p>")
      ),
      tabItem("Heterogeneitywithinhistologictypes",
              girafeOutput('figure_entropy'),
              hr(),
              h6("Figure. Unsupervised clustering based on oncogenic mutations. For each histology, the percentage of cases belonging to one of the clusters is shown as a bar graph. Characteristic oncogenic mutations were found in each cluster. There was a tendency for each histologic type to cluster in specific clusters. Heterogeneity of genetic variation within histologic types was assessed by Shannon's entropy. Low Shannon entropy values indicate that the genetic mutation pattern of the tumor is uniform, while high values indicate diversity.")
      ),
      tabItem("Frequencyofpatientswithtargetedtherapy",
              htmlOutput("select_figure_drug_evidence_var"),
              girafeOutput('figure_drug_evidence',
                           height = "1000px",
                           width = "1000px"),
              hr(),
              DT::dataTableOutput("table_drug_evidence_y"),
              DT::dataTableOutput("table_drug_evidence_z"),
              DT::dataTableOutput("table_evidence"),
              h6("Figure. Level of evidence for targeted therapy for detected gene mutations. The highest level of evidence was extracted for each patient. Evidence levels of C-CAT are defined as A for biomarkers that predict a response to Japanese Pharmaceuticals and Medical Devices Agency (PMDA)– or FDA-approved therapies or are described in professional guidelines, B for biomarkers that predict a response based on well-powered studies with consensus of experts in the field, C for biomarkers that predict a response to therapies approved by the PMDA or FDA in another type of tumor or that predict a response based on clinical studies, D for biomarkers that predict a response based on case reports, E for biomarkers that show plausible therapeutic significance based on preclinical studies.")
      ),
      tabItem("SurvivalafterCGP",
              h5("Survival difference will be evaluated with restricted mean survival time in this section. Analysis with hazard ratio is also provided in 'Overall survival with risk-set adjustment' section (Survival analysis start date = CGP test date)."),
              h6("Figure. Survival analysis after CGP test using the conventional Kaplan–Meier estimator, log–rank test were undertaken with survival package for R. EP: expert panel. RMST, restricted mean survival time."),
              htmlOutput("select_color_var_surv_CGP"),
              hr(),
              plotOutput('figure_surv_CGP',
                         height = "800px",
                         width = "800px"),
              hr(),
              h6("Raw data"),
              DT::dataTableOutput("figure_survival_CGP_3_rawdata")
      ),
      tabItem("SurvivalandtreatmentafterCGP",
              h6("Figure. Survival analysis after CGP test using the conventional Kaplan–Meier estimator, log–rank test were undertaken with survival package for R. EP: expert panel. RMST, restricted mean survival time."),
              plotOutput('figure_survival_CGP_1',
                         height = "800px",
                         width = "800px"),
              hr(),
              fluidRow(
                column(3,
                       h6("Group 1"),
                       htmlOutput("select_gene_survival_CGP_1_P_1"),
                       htmlOutput("select_gene_survival_CGP_1_P_2"),
                       htmlOutput("select_gene_survival_CGP_1_W"),
                       htmlOutput("select_gene_survival_CGP_1_D"),
                       htmlOutput("select_gene_survival_CGP_1_ND"),
                       htmlOutput("select_gene_survival_CGP_1_L"),
                       htmlOutput("select_gene_survival_CGP_1_R"),
                       htmlOutput("select_gene_survival_CGP_1_EP_option"),
                       htmlOutput("select_gene_survival_CGP_1_EP_treat")
                ),
                column(3,
                       h6("   "),
                       htmlOutput("select_gene_survival_CGP_1_A"),
                       htmlOutput("select_gene_survival_CGP_1_S"),
                       htmlOutput("select_gene_survival_CGP_1_H"),
                       htmlOutput("select_gene_survival_CGP_1_C"),
                       htmlOutput("select_gene_survival_CGP_1_M"),
                       htmlOutput("select_gene_survival_CGP_1_P")
                ),
                column(3,
                       h6("Group 2"),
                       htmlOutput("select_gene_survival_CGP_2_P_1"),
                       htmlOutput("select_gene_survival_CGP_2_P_2"),
                       htmlOutput("select_gene_survival_CGP_2_W"),
                       htmlOutput("select_gene_survival_CGP_2_D"),
                       htmlOutput("select_gene_survival_CGP_2_ND"),
                       htmlOutput("select_gene_survival_CGP_2_L"),
                       htmlOutput("select_gene_survival_CGP_2_R"),
                       htmlOutput("select_gene_survival_CGP_2_EP_option"),
                       htmlOutput("select_gene_survival_CGP_2_EP_treat")
                ),
                column(3,
                       h6("   "),
                       htmlOutput("select_gene_survival_CGP_2_A"),
                       htmlOutput("select_gene_survival_CGP_2_S"),
                       htmlOutput("select_gene_survival_CGP_2_H"),
                       htmlOutput("select_gene_survival_CGP_2_C"),
                       htmlOutput("select_gene_survival_CGP_2_M"),
                       htmlOutput("select_gene_survival_CGP_2_P")
                )
              ),
              br(),
              htmlOutputWithPopover(
                "select_propensity_survival_CGP",
                "指定した因子で傾向スコアマッチングを行います",
                "ロジットPSでのマッチ、caliper=0.2 * sd_logit,ペア単位の2000回のブートストラップでRMSTの差の信頼区間を推定 "
              ),
              br(),
              downloadButton("dl_love_plot_PSM", "Download love plot of PS-matching"),
              hr(),
              h6("To reduce confounding between the two groups, we performed propensity score matching. The propensity score was estimated using a logistic regression model including prespecified clinically relevant covariates (CGP platform, sex, age, PS, histology, treatment lines before CGP, and the best treatment effect before CGP). Patients were matched 1:1 using nearest‐neighbor matching without replacement on the logit of the propensity score (MatchIt package, method = “nearest”, distance = “logit”). A caliper width of 0.2 on the logit scale was applied to restrict matches to comparable individuals. Matched sets were identified using the MatchIt subclass variable, and each subclass was treated as a matched pair for following analyses."),
              h6("Covariate balance before and after matching was evaluated using standardized mean differences (SMDs) with the cobalt package. Adequate balance was defined a priori as an absolute SMD < 0.1 for all covariates. Balance diagnostics were visualized using Love plots. The maximum absolute SMD after matching was additionally reported to provide a single summary measure of balance."),
              h6("Because propensity score matching induces dependence within matched pairs, confidence intervals for the RMST difference were obtained using nonparametric bootstrap resampling at the matched‐pair level. Specifically, matched pairs were resampled with 2000-time replacement, RMST differences were recalculated for each bootstrap replicate, and the 2.5th and 97.5th percentiles of the bootstrap distribution were used to derive a two‐sided 95% confidence interval."),
              hr(),
              br(),
              br(),
              br(),
              br(),
              br(),
              hr()
      ),
      tabItem("SurvivalafterCGPandmutationsforestplot",
              h5("Survival difference is evaluated with restricted mean survival time in this section. Analysis with hazard ratio will be also provided in 'Overall survival with risk-set adjustment' section (Survival analysis start date = CGP test date)."),
              h6("Figure. Suvival periods after CGP and gene mutations estimated with conventional Kaplan-Meier estimator. Restricted mean survival time in two years (days) were estimated with survRM2 package in R."),
              plotOutput('figure_survival_CGP_3',
                         height = "1000px",
                         width = "1200px"),
              hr(),
              plotOutput('figure_survival_CGP_4',
                         height = "3600px",
                         width = "2000px"),
              hr(),
              DT::dataTableOutput("figure_survival_CGP_3_data")
      ),
      tabItem("HazardratioforsurvivalafterCGP_1",
              fluidRow(
                column(4,
                       htmlOutput("select_figure_survival_CGP_5_1_var")),
                column(4,
                       htmlOutput("select_figure_survival_CGP_5_1_max_samples")
                ),
                column(4,
                       actionButton("figure_survival_CGP_5_1_reload", "Reload")
                )
              ),
              htmlOutput(""),
              hr(),
              gt_output('figure_survival_CGP_5_1'),
              hr(),
              h6("Table. Univariable and multiple variable regression analysis were performed with gtsummary package for R. Factors in the multivariable analysis were selected by variable decreasing method based on the Akaike information criterion with stat package for R. Variables with variance inflation factor > 10 were excluded to eliminate multicollinearity.")
      ),
      tabItem("Survival_treatment_reach",
              h6(paste0("Figure. Treatment reach rate and survival period after CGP. We estimated the cumulative incidence function (CIF) using Gray’s method for competing risk analysis.",
                        "Treatment initiation was defined as the primary outcome, death was considered a competing event, and censoring was treated as informative.",
                        "The cumulative incidence at each time point t was calculated as the probability of experiencing the specified event by t, reflecting the conditional probability in the presence of competing events.",
                        "If the recommended treatment was confirmed to have been administered but the treatment initiation date was unknown, the treatment was assumed to have started at the median of the entire observation period.",
                        "Statistical analyses were performed using the cmprsk package in R, and CIFs were estimated with 95% confidence intervals.",
                        "We used competing risk analysis instead of the conventional Kaplan–Meier method for the following reasons:",
                        "Presence of competing events: Patients who die permanently lose the opportunity to receive treatment, making death a competing event.",
                        "Bias avoidance: Treating deaths as simple censored observations in Kaplan–Meier analysis could overestimate the treatment initiation rate.",
                        "Clinical interpretability: CIFs provide probabilities that more accurately reflect event occurrence as observed in real-world clinical settings.")),
              fluidRow(
                column(6,
                       girafeOutput('figure_survival_treatment_reach_1',
                                    height = "800px",
                                    width = "800px")
                ),
                column(6,
                       girafeOutput('figure_survival_treatment_reach_2',
                                    height = "800px",
                                    width = "800px")
                )
              ),
              hr(),
              fluidRow(
                column(3,
                       htmlOutput("select_gene_survival_reach_1_P_1"),
                       htmlOutput("select_gene_survival_reach_1_P_2"),
                       htmlOutput("select_gene_survival_reach_1_W"),
                       htmlOutput("select_gene_survival_reach_1_D"),
                       htmlOutput("select_gene_survival_reach_1_ND"),
                       htmlOutput("select_gene_survival_reach_1_L"),
                       htmlOutput("select_gene_survival_reach_1_R"),
                       htmlOutput("select_gene_survival_reach_1_EP_option")
                ),
                column(3,
                       htmlOutput("select_gene_survival_reach_1_A"),
                       htmlOutput("select_gene_survival_reach_1_S"),
                       htmlOutput("select_gene_survival_reach_1_H"),
                       htmlOutput("select_gene_survival_reach_1_C"),
                       htmlOutput("select_gene_survival_reach_1_M"),
                       htmlOutput("select_gene_survival_reach_1_P"),
                       htmlOutput("select_gene_survival_reach_1_Panel")
                )
              ),
              br(),
              br(),
              br(),
              br(),
              hr()
      ),
      tabItem("FactorsleadingtoTreatmentpre-CGPNomogram",
              h6("Figure. Based on clinical information, a nomogram was developed to predict treatment reach. The nomogram was created with the lrm function of the rms package for R with a setting of penalty=1. Nagelkerke R2 was calculated with blorr package for R. Possible sampling bias was corrected with 500-time bootstrap sampling and then concordance index was estimated. Best_effect: the best treatment effect of CTx before CGP."),
              fluidRow(
                column(4,
                       htmlOutput("select_nomogram_type")),
                column(4,
                       htmlOutput("select_nomogram_max_samples")
                ),
                column(4,
                       actionButton("CGP_benefit_analysis_reload", "Reload"))
              ),
              hr(),
              plotOutput('figure_nomogram_treat_pre_CGP',
                         height = "2000px",
                         width = "1000px"),
              hr(),
              verbatimTextOutput("histology_nomogram_pre_CGP")
      ),
      tabItem("FactorsleadingtoTreatmentpre-CGPOddsratio",
              h6("Table. Univariable and multiple variable regression analysis were performed with gtsummary package for R. Factors selected in the multivariable analysis were selected by variable decreasing method based on the Akaike information criterion with stat package for R. Variables with variance inflation factor > 10 were excluded to eliminate multicollinearity."),
              htmlOutput("select_nomogram_type_2"),
              hr(),
              gt_output('figure_survival_CGP_7'),
              hr(),
      ),
      tabItem("FactorsleadingtoTreatmentdecisioncurve",
              h6("Figure. Decision curve analysis was performed to verify the usefulness of the nomogram with dcurves package for R. Ten-fold cross-validation was performed to prevent overfitting. If the blue line is located above the other lines, then the nomogram-based decision to perform CGP testing may be worthwhile."),
              hr(),
              fluidRow(
                column(4,htmlOutputWithPopover(
                "select_nomogram_type_3",
                "治療到達、推奨治療提案、Evidence level A/B/Cの変異検出についての予測モデルの有用性を評価",
                HTML('Decision curve analysisの詳細な説明は<a href="https://github.com/MANO-B/FELIS/blob/main/decision_curve_analysis.md" target="_blank" rel="noopener noreferrer">https://github.com/MANO-B/FELIS/blob/main/decision_curve_analysis.md</a>をご確認ください。',),
                trigger = "click",
                placement = "bottom"
              ))),
              hr(),
              plotOutput('figure_DCA_pre_CGP_1',
                         height = "1000px",
                         width = "1300px"),
              hr(),
              gt_output('table_DCA_pre_CGP_2')
      ),
      tabItem("ROCnomogram",
              h6("Figure. Predicted treatment reach rate and Receiver Operatorating Characteristic curve of the nomogram using pre-CGP information by pROC package for R. The nomogram was based on logistic regression analysis, random forest model, and LightGBM model, all of which calculated sensitivity and specificity using prediction results from 5-fold cross-validation and plotted ROC curves. The random forest model and LightGBM model used 5-fold cross-validation on a single training set and performed a grid search with n=8 for each parameter to determine the optimal parameters."),
              hr(),
              htmlOutput("select_nomogram_type_4"),
              hr(),
              plotOutput('figure_DCA_pre_CGP_ROC',
                         height = "3000px",
                         width = "1000px"),
              hr(),
              h6("Download raw data"),
              DT::dataTableOutput("table_prediction"),
              hr(),
              p("Not shown when machine learning is not performed"),
              plotOutput('figure_DCA_pre_CGP_Machine_learning',
                         height = "2000px",
                         width = "1500px"),
              hr()
      ),
      tabItem(tabName = "Input_data",
              htmlOutput("select_prediction_var"),
              hr(),
              fluidRow(
                column(4,strong("Clinical information 1"),
                       hr(),
                       htmlOutput("input_age"),
                       htmlOutput("input_sex"),
                       htmlOutput("input_Smoking_history"),
                       htmlOutput("input_Alcoholic_history"),
                       # htmlOutput("input_Enroll_date"),
                       # htmlOutput("input_Enroll_date_raw"),
                       htmlOutput("input_time_diagnosis_enroll"),
                       htmlOutput("input_PS"),
                       #htmlOutput("input_PS_raw"),
                       br(),
                       hr()
                ),
                column(4,strong("Clinical information 2"),
                       hr(),
                       htmlOutput("input_Lines"),
                       #htmlOutput("input_Lines_raw"),
                       htmlOutput("input_Best_effect"),
                       htmlOutput("input_organ"),
                       #htmlOutput("input_organ_raw"),
                       htmlOutput("input_Panel"),
                       #htmlOutput("input_Panel_raw"),
                       br(),
                       hr()
                ),
                column(4,strong("Metastasis information"),
                       hr(),
                       htmlOutput("input_Lymph_met"),
                       htmlOutput("input_Brain_met"),
                       htmlOutput("input_Lung_met"),
                       htmlOutput("input_Bone_met"),
                       htmlOutput("input_Liver_met"),
                       br(),
                       hr()
                )
              ),
              br(),
              hr(),
              br(),
              actionButton("start_prediction", "Start prediction"),
              hr(),
              htmlOutput("prediction_nomogram"),
              br(),
              br(),
              hr()
      ),
      tabItem("SurvivalandtreatmentafterCTx",
              h6("Figure. Overall survival after the first survival-prolonging chemotherapy after adjusting for left-truncation bias. To evaluate the association between prognostic factors and survival, a risk-set (number at risk) adjustment model was applied to adjust for left-truncation bias with survival package. Note that the analysis assumes quasi-independent left-truncation (conditional Kendall tau = 0)."),
              h6("  Take care of left-truncation bias."),
              fluidRow(
                column(4,
                       htmlOutputWithPopover(
                         "select_color_var_surv_CTx_3",
                         "左側切断バイアスの補正の有無：変更後はリロード",
                         "R survival packageによるrisk-set adjustment。Kendall tauが0のindependent truncationの場合に妥当性があります。"
                       ),
                       htmlOutputWithPopover(
                         "select_color_var_surv_CTx_1",
                         "生存期間の測定開始日：変更後はリロード",
                         "CGP test dateにした場合は左側切断のないCGP検査後の生存期間解析です。"
                       ),
                       actionButton("color_var_surv_CTx_reload", "Reload")),
                column(4,
                       htmlOutput("select_color_var_surv_CTx_2")
                ),
                column(4,
                       )
              ),
              hr(),
              plotOutput('figure_survival_CTx_1_2',
                         height = "800px",
                         width = "800px"),
              hr(),
              h6("Raw data"),
              DT::dataTableOutput("figure_survival_CTx_2_data_2_raw")
      ),
      tabItem("SurvivalandtreatmentafterCTx2",
              h6("Figure. Overall survival after the first survival-prolonging chemotherapy after adjusting for left-truncation bias. To evaluate the association between prognostic factors and survival, a risk-set (number at risk) adjustment model was applied to adjust for left-truncation bias with survival package. Note that the analysis assumes quasi-independent left-truncation (conditional Kendall tau = 0)."),
              h6("  Take care of left-truncation bias."),
              plotOutput('figure_survival_CTx_interactive_1',
                         height = "800px",
                         width = "800px"),
              hr(),
              fluidRow(
                column(3,
                       h6("Group 1"),
                       htmlOutput("select_gene_survival_interactive_1_P_1"),
                       htmlOutput("select_gene_survival_interactive_1_P_2"),
                       htmlOutput("select_gene_survival_interactive_1_W"),
                       htmlOutput("select_gene_survival_interactive_1_D"),
                       htmlOutput("select_gene_survival_interactive_1_L"),
                       htmlOutput("select_gene_survival_interactive_1_R")
                ),
                column(3,
                       h6("   "),
                       htmlOutput("select_gene_survival_interactive_1_A"),
                       htmlOutput("select_gene_survival_interactive_1_S"),
                       htmlOutput("select_gene_survival_interactive_1_H"),
                       htmlOutput("select_gene_survival_interactive_1_C"),
                       htmlOutput("select_gene_survival_interactive_1_M"),
                       htmlOutput("select_gene_survival_interactive_1_P")
                ),
                column(3,
                       h6("Group 2"),
                       htmlOutput("select_gene_survival_interactive_2_P_1"),
                       htmlOutput("select_gene_survival_interactive_2_P_2"),
                       htmlOutput("select_gene_survival_interactive_2_W"),
                       htmlOutput("select_gene_survival_interactive_2_D"),
                       htmlOutput("select_gene_survival_interactive_2_L"),
                       htmlOutput("select_gene_survival_interactive_2_R")
                ),
                column(3,
                       h6("   "),
                       htmlOutput("select_gene_survival_interactive_2_A"),
                       htmlOutput("select_gene_survival_interactive_2_S"),
                       htmlOutput("select_gene_survival_interactive_2_H"),
                       htmlOutput("select_gene_survival_interactive_2_C"),
                       htmlOutput("select_gene_survival_interactive_2_M"),
                       htmlOutput("select_gene_survival_interactive_2_P")
                )
              ),
              br(),
              br(),
              br(),
              br(),
              hr()
      ),
      tabItem("Geneticvariantsandsurvivalforestplot2",
              h6("Figure. Overall survival after the first survival-prolonging chemotherapy after adjusting for left-truncation bias. To evaluate the association between oncogenic mutations and survival, a risk-set adjustment model was performed to adjust for left-truncation bias with survival package."),
              h6("  Take care of left-truncation bias."),
              plotOutput('figure_survival_CTx_2_2',
                         height = "1000px",
                         width = "1250px"),
              hr(),
              plotOutput('figure_survival_CTx_3_2',
                         height = "3600px",
                         width = "2000px"),
              hr(),
              DT::dataTableOutput("figure_survival_CTx_2_data_2")
      ),
      tabItem("Diagnosisandsurvivalforestplot2",
              plotOutput('figure_survival_CTx_4_2',
                         height = "1200px",
                         width = "1000px"),
              hr()
      ),
      tabItem("DiagnosisandsurvivalKM-curve2",
              h6("  Take care of left-truncation bias."),
              plotOutput('figure_survival_CTx_5_2',
                         height = "1200px",
                         width = "1000px"),
              hr(),
              h6("Figure. Overall survival after the first survival-prolonging chemotherapy after adjusting for left-truncation bias. To evaluate the association between oncogenic mutations and survival, a risk-set adjustment model was performed to adjust for left-truncation bias with survival package.")
      ),
      tabItem("HazardratioforsurvivalafterCTx_1",
              h6("  It takes minutes."),
              h6("  Take care of left-truncation bias."),
              plotOutput('figure_survival_CTx_5_1_2',
                         height = "1200px",
                         width = "1000px"),
              hr(),
              h6("Figure. Hazard ratio estimated by cox model with survival package.")
      ),
      tabItem("Survivalcorrectedforleft-truncationbias",
              h6("  It takes minutes."),
              fluidRow(
                column(4,
                       htmlOutput("select_figure_Bayes_max_samples"),
                       h6("This setting also applies to Bayesian estimation in other tabs."),
                       htmlOutput("select_Bayes_basic_var"),
                       actionButton("figure_Bayes_reload", "Reload")
                ),
                column(4,
                       actionButton("start_Bayes_basic_prediction", "Start overall cohort prediction")
                ),
                column(4
                )
              ),
              hr(),
              plotOutput('figure_survival_CTx_1',
                         height = "2500px",
                         width = "1000px"),
              hr(),
      ),
      tabItem("BayesCustom",
              plotOutput('figure_survival_CTx_compare',
                         height = "800px",
                         width = "800px"),
              hr(),
              actionButton("start_Bayes_compare_prediction", "Start 2-group prediction"),
              htmlOutput("select_Bayes_prediction_var"),
              hr(),
              fluidRow(
                column(3,
                       h6("Group 1"),
                       htmlOutput("select_gene_survival_interactive_1_P_1_Bayes"),
                       htmlOutput("select_gene_survival_interactive_1_P_2_Bayes"),
                       htmlOutput("select_gene_survival_interactive_1_W_Bayes"),
                       htmlOutput("select_gene_survival_interactive_1_D_Bayes"),
                       htmlOutput("select_gene_survival_interactive_1_L_Bayes"),
                       htmlOutput("select_gene_survival_interactive_1_R_Bayes")
                ),
                column(3,
                       h6("   "),
                       htmlOutput("select_gene_survival_interactive_1_A_Bayes"),
                       htmlOutput("select_gene_survival_interactive_1_S_Bayes"),
                       htmlOutput("select_gene_survival_interactive_1_H_Bayes"),
                       htmlOutput("select_gene_survival_interactive_1_C_Bayes"),
                       htmlOutput("select_gene_survival_interactive_1_M_Bayes"),
                       htmlOutput("select_gene_survival_interactive_1_P_Bayes")
                ),
                column(3,
                       h6("Group 2"),
                       htmlOutput("select_gene_survival_interactive_2_P_1_Bayes"),
                       htmlOutput("select_gene_survival_interactive_2_P_2_Bayes"),
                       htmlOutput("select_gene_survival_interactive_2_W_Bayes"),
                       htmlOutput("select_gene_survival_interactive_2_D_Bayes"),
                       htmlOutput("select_gene_survival_interactive_2_L_Bayes"),
                       htmlOutput("select_gene_survival_interactive_2_R_Bayes")
                ),
                column(3,
                       h6("   "),
                       htmlOutput("select_gene_survival_interactive_2_A_Bayes"),
                       htmlOutput("select_gene_survival_interactive_2_S_Bayes"),
                       htmlOutput("select_gene_survival_interactive_2_H_Bayes"),
                       htmlOutput("select_gene_survival_interactive_2_C_Bayes"),
                       htmlOutput("select_gene_survival_interactive_2_M_Bayes"),
                       htmlOutput("select_gene_survival_interactive_2_P_Bayes")
                )
              ),
              br(),
              br(),
              br(),
              br(),
              hr()
      ),
      tabItem("GeneticvariantsandsurvivalKM-curve",
              h6("  It takes minutes."),
              fluidRow(
                column(4,
                       htmlOutput("select_Bayes_mutation_var"),
                       htmlOutput("select_Bayes_mutation_no"),
                ),
                column(4,
                       actionButton("start_Bayes_mutation_prediction", "Start prediction"),
                ),
                column(4
                )
              ),
              hr(),
              plotOutput('figure_survival_CTx_3',
                         height = "3600px",
                         width = "1500px"),
              hr(),
              plotOutput('figure_survival_CTx_2',
                         height = "1000px",
                         width = "1250px"),
              hr(),
              DT::dataTableOutput("figure_survival_CTx_2_data"),
              h6(paste0("Figure. Overall survival after the first survival-prolonging chemotherapy after adjusting for left-truncation bias. To evaluate the association between oncogenic mutations and survival, a Bayesian survival simulation based on a semi-independent, two-hit model was performed to adjust for left-truncation bias. ",
                        "Two survival curves from the date of commencement of chemotherapy for prolonging survival to the date of CGP and from the CGP testing date to the date of death were fitted with Weibull distribution and log-logistic distribution, respectively. The overall survival curve from the first chemotherapy induction was approximated by merging these survival curves. ",
                        "Survival curves were obtained from each of the 8000 iterations of inference, and the median survival and 95% equal-tailed CIs were calculated. Bayesian inference was performed with the rstan package for R. P value of conditional Kendall tau statistics was calculated, and the survival curves were adjusted for length bias, using a structural transformation method with tranSurv package for R.")),
              HTML("<p>Tamura T, et al., Selection bias due to delayed comprehensive genomic profiling in Japan. Cancer Science, 2022. <a href='https://doi.org/10.1111/cas.15651'>Link for the paper</a></p>")
      ),
      tabItem("DiagnosisandsurvivalKM-curve",
              h6("  It takes minutes."),
              fluidRow(
                column(4,
                       htmlOutput("select_Bayes_diagnosis_var"),
                       htmlOutput("select_Bayes_diagnosis_no"),
                ),
                column(4,
                       actionButton("start_Bayes_diagnosis_prediction", "Start prediction"),
                ),
                column(4
                )
              ),
              hr(),
              plotOutput('figure_survival_CTx_5',
                         height = "3200px",
                         width = "3000px"),
              hr(),
              plotOutput('figure_survival_CTx_4',
                         height = "600px",
                         width = "1250px"),
              hr(),
              DT::dataTableOutput("figure_survival_CTx_5_table"),
              DT::dataTableOutput("figure_survival_CTx_4_data"),
              h6("Figure. Overall survival after the first survival-prolonging chemotherapy after adjusting for left-truncation bias. To evaluate the association between oncogenic mutations and survival, a Bayesian survival simulation based on a semi-independent, two-hit model was performed to adjust for left-truncation bias. Two survival curves from the date of commencement of chemotherapy for prolonging survival to the date of CGP and from the CGP testing date to the date of death were fitted with Weibull distribution and log-logistic distribution, respectively. The overall survival curve from the first chemotherapy induction was approximated by merging these survival curves. Survival curves were obtained from each of the 8000 iterations of inference, and the median survival and 95% equal-tailed CIs were calculated. Bayesian inference was performed with the rstan package for R."),
              HTML("<p>Tamura T, et al., Selection bias due to delayed comprehensive genomic profiling in Japan. Cancer Science, 2022. <a href='https://doi.org/10.1111/cas.15651'>Link for the paper</a></p>")
      ),
      tabItem("SurvivalSimurationKMCurve",
              fluidRow(
                column(3,strong("Simulation settings"),
                       hr(),
                       htmlOutput("input_shape_llogis"),
                       htmlOutput("input_scale_llogis"),
                       br(),
                       hr()
                ),
                column(3,strong(" "),
                       hr(),
                       htmlOutput("input_shape_CGP"),
                       htmlOutput("input_scale_CGP"),
                       br(),
                       hr()
                ),
                column(3,strong(" "),
                       hr(),
                       htmlOutput("input_shape_weibull"),
                       htmlOutput("input_scale_weibull"),
                       br(),
                       hr()
                ),
                column(3,strong(" "),
                       hr(),
                       htmlOutput("input_censor_rate"),
                       htmlOutput("input_test_patients"),
                       br(),
                       hr()
                )
              ),
              br(),
              hr(),
              actionButton("survival_simuration", "Left-truncation bias adjustment simulation"),
              br(),
              hr(),
              plotOutput('figure_bias_correction_simulation',
                         height = "1000px",
                         width = "1000px"),
              hr(),
              h6("Figure. Overall survival after the first survival-prolonging chemotherapy. "),
              HTML("<p>Tamura T, et al., Selection bias due to delayed comprehensive genomic profiling in Japan. Cancer Science, 2022. <a href='https://doi.org/10.1111/cas.15651'>Link for the paper</a></p>"),
              br(),
              verbatimTextOutput("text_bias_correction_simulation"),
              br(),
              br(),
              hr()
      ),

      tabItem("Drugperpatient",
              DT::dataTableOutput("table_drug_all_summary"),
              hr(),
              actionButton('show_download_confirm_drug', 'Download all drug data',
                           icon = icon("download"),
                           onclick = "confirmDownload('download_drug_data', 'FELIS downloaded raw data in Drug usage data tab')",
                           class = "btn-primary"),

              # 非表示のダウンロードボタン
              div(style = "visibility:hidden; position:absolute;",
                  downloadButton('download_drug_data', 'Hidden Download')
              )
      ),
      tabItem("Drugusebylineoftreatment",
              fluidRow(
                column(3,
                       strong("Drug response analysis"),
                       hr(),
                       htmlOutput("select_target_line"),
                       hr(),
                       htmlOutput("select_ToT_TTF"),
                       br(),
                       htmlOutput("select_multiple_used_drug"),
                       br(),
                       htmlOutput("select_drug_multi_gene"),
                       br(),
                       htmlOutput("select_drug"),
                       br(),
                       htmlOutput("select_drug_group_1"),
                       br(),
                       htmlOutput("select_drug_group_1_name"),
                       br(),
                       htmlOutput("select_drug_group_2"),
                       br(),
                       htmlOutput("select_drug_group_2_name"),
                       br(),
                       htmlOutput("select_drug_group_3"),
                       br(),
                       htmlOutput("select_drug_group_3_name"),
                       br(),
                       htmlOutput("select_drug_group_4"),
                       br(),
                       htmlOutput("select_drug_group_4_name"),
                       br()
                ),
                column(9,
                       h4("Overall drug usage"),
                       fluidRow(
                         column(3,
                                htmlOutput("select_table_var_drug_1")
                         ),
                         column(3,
                                htmlOutput("select_drug_table_layout")
                         ),
                         column(3,
                                htmlOutput("select_minimum_courses"),
                                actionButton("drug_table_reload", "Reload dataset with a new threshold")
                         )
                       ),
                       hr(),
                       gt_output("table_drug_all_1")
                )
              ),
              hr()
      ),
      tabItem("UseofdesignatedlinesanddrugswithToTinformationbymutationpattern",
              fluidRow(
                column(3,
                       htmlOutput("select_table_var_drug_2")
                ),
                column(3,
                       h6("Patients without treatment time excluded in treatment time dataset"),
                       h6("Patients with RECIST-NE excluded in objective response dataset"),
                       h6("Patients without treatment time  or with RECIST-NE excluded in adverse effect dataset"),
                       htmlOutput("select_table_var_drug_dataset")
                )
              ),
              hr(),
              htmlOutput("select_gene_survival_interactive_1_P_1"),

              hr(),
              gt_output("table_drug_ToT_1")
      ),
      tabItem("Timeontreatmentandpre-treatmentforthespecifiedtreatmentscatterplot",
              htmlOutput("select_ToT_var_1"),
              hr(),
              plotOutput('figure_drug_1',
                         height = "800px",
                         width = "800px"),
              hr(),
              plotOutput('figure_drug_3_forest',
                         height = "600px",
                         width = "1000px"),
              hr(),
              DT::dataTableOutput("figure_drug_3_forest_table"),
              h6(paste0("Figure. Time on treatment analysis for the survival-prolonging chemotherapy. ",
                        "Time on treatment represents the period from the start date of chemotherapy to the end date; ",
                        "if the patient was on medication at the time of CGP testing, the patient was censored; ",
                        "otherwise, the patient was terminated, and a survival curve was generated using the Kaplan-Meier method.\n"
              ))
      ),
      tabItem("ToT_interactive",
              h6("Figure. Treatment time."),
              hr(),
              plotOutput('figure_survival_drug_interactive_1',
                         height = "800px",
                         width = "800px"
              ),
              hr(),
              fluidRow(
                column(3,
                       h6("Group 1"),
                       htmlOutput("select_gene_survival_drug_interactive_1_P_1"),
                       htmlOutput("select_gene_survival_drug_interactive_1_P_2"),
                       htmlOutput("select_gene_survival_drug_interactive_1_W"),
                       htmlOutput("select_gene_survival_drug_interactive_1_D"),
                       htmlOutput("select_gene_survival_drug_interactive_1_L"),
                       htmlOutput("select_gene_survival_drug_interactive_1_AE"),
                       htmlOutput("select_gene_survival_drug_interactive_1_AED")
                ),
                column(3,
                       h6("   "),
                       htmlOutput("select_gene_survival_drug_interactive_1_A"),
                       htmlOutput("select_gene_survival_drug_interactive_1_S"),
                       htmlOutput("select_gene_survival_drug_interactive_1_H"),
                       htmlOutput("select_gene_survival_drug_interactive_1_C"),
                       htmlOutput("select_gene_survival_drug_interactive_1_M"),
                       htmlOutput("select_gene_survival_drug_interactive_1_R"),
                       htmlOutput("select_gene_survival_drug_interactive_1_AEG"),
                       htmlOutput("select_gene_survival_drug_interactive_1_PD_L1"),
                       htmlOutput("select_gene_survival_drug_interactive_1_EAE")
                ),
                column(3,
                       h6("Group 2"),
                       htmlOutput("select_gene_survival_drug_interactive_2_P_1"),
                       htmlOutput("select_gene_survival_drug_interactive_2_P_2"),
                       htmlOutput("select_gene_survival_drug_interactive_2_W"),
                       htmlOutput("select_gene_survival_drug_interactive_2_D"),
                       htmlOutput("select_gene_survival_drug_interactive_2_L"),
                       htmlOutput("select_gene_survival_drug_interactive_2_AE"),
                       htmlOutput("select_gene_survival_drug_interactive_2_AED")
                ),
                column(3,
                       h6("   "),
                       htmlOutput("select_gene_survival_drug_interactive_2_A"),
                       htmlOutput("select_gene_survival_drug_interactive_2_S"),
                       htmlOutput("select_gene_survival_drug_interactive_2_H"),
                       htmlOutput("select_gene_survival_drug_interactive_2_C"),
                       htmlOutput("select_gene_survival_drug_interactive_2_M"),
                       htmlOutput("select_gene_survival_drug_interactive_2_R"),
                       htmlOutput("select_gene_survival_drug_interactive_2_AEG"),
                       htmlOutput("select_gene_survival_drug_interactive_2_PD_L1"),
                       htmlOutput("select_gene_survival_drug_interactive_2_EAE")
                )
              ),
              br(),
              br(),
              br(),
              br(),
              hr()
      ),
      tabItem("TimeontreatmentbygenemutationclusterKM-curve",
              plotOutput('figure_drug_cluster_2',
                         height = "700px",
                         width = "1100px"),
              hr(),
              plotOutput('figure_drug_cluster',
                         height = "3500px",
                         width = "2000px"),
              DT::dataTableOutput("figure_drug_cluster_table")
      ),
      tabItem("Timeontreatmentbymutatedgenesforestplot",
              plotOutput('figure_drug_4',
                         height = "700px",
                         width = "1100px"),
              hr(),
              plotOutput('figure_drug_5',
                         height = "3500px",
                         width = "2000px"),
              DT::dataTableOutput("figure_drug_4_table")
      ),
      tabItem("TimeontreatmentandmutationsofinterestKM-curve",
              plotOutput('figure_drug_pattern',
                         height = "3500px",
                         width = "2000px")
      ),
      tabItem("Hazardratioontimeontreatment_1",
              gt_output('figure_drug_6_1'),
              hr(),
              h6("Table. Univariable and multiple variable regression analysis were performed with gtsummary package for R. Factors selected in the multivariable analysis were selected by variable decreasing method based on the Akaike information criterion with stat package for R. Variables with variance inflation factor > 10 were excluded to eliminate multicollinearity.")
      ),
      tabItem("Hazardratioontimeontreatment_2",
              gt_output('figure_drug_6_2'),
              hr(),
              h6("Table. Univariable and multiple variable regression analysis were performed with gtsummary package for R. Factors selected in the multivariable analysis were selected by variable decreasing method based on the Akaike information criterion with stat package for R. Variables with variance inflation factor > 10 were excluded to eliminate multicollinearity.")
      ),
      tabItem("VolcanoPlot_ToT_1",
              h6("Volcano plots for frequent regimens"),
              hr(),
              plotOutput('figure_volcano_ToT_1',
                         height = "600px",
                         width = "800px"),
              hr(),
              htmlOutput("select_table_var_volcano_1"),
              hr(),
              h6("Raw data for volcano plots of all regimens"),
              h6("Hazard ratio and p-value were calculated by cox regression model concerning pathology."),
              h6("All regimens, and genes in which more than or equal to 8 of the treated patients had mutations were included in the analysis."),
              hr(),
              DT::dataTableOutput("table_volcano_ToT")
      ),
      tabItem("VolcanoPlotORR_1",
              h6("Volcano plots for frequent regimens"),
              hr(),
              plotOutput('figure_volcano_1',
                         height = "600px",
                         width = "800px"),
              hr(),
              htmlOutput("select_table_var_volcano_2"),
              hr(),
              DT::dataTableOutput("table_volcano"),
              hr(),
              h6("Figure. Patients with objective response data treated with the specified drugs were divided into two groups: those who obtained an Objective response and those who did not, and their odds ratios and p-values were calculated and a volcano plot was plotted."),
              h6("Odds ratio and p-value were calculated by logistic regression model concerning pathology."),
              h6("All regimens, and genes in which more than or equal to 8 of the treated patients had mutations were included in the analysis."),
              h6("Objective response: CR or PR.")
      ),
      tabItem("Oddsratioonobjectiveresponserate_1",
              gt_output('figure_drug_7_1'),
              hr(),
              h6("Table. Univariable and multiple variable regression analysis were performed with gtsummary package for R. Factors selected in the multivariable analysis were selected by variable decreasing method based on the Akaike information criterion with stat package for R. Variables with variance inflation factor > 10 were excluded to eliminate multicollinearity."),
              h6("Candidate genes: arbitrary selected genes, the five most frequently mutated genes, genes with significance in odds ratio of objective response rate or time on treatment."),
              h6("Objective response: CR or PR.")
      ),
      tabItem("Oddsratioonobjectiveresponserate_2",
              gt_output('figure_drug_7_2'),
              hr(),
              h6("Table. Univariable and multiple variable regression analysis were performed with gtsummary package for R. Factors selected in the multivariable analysis were selected by variable decreasing method based on the Akaike information criterion with stat package for R. Variables with variance inflation factor > 10 were excluded to eliminate multicollinearity."),
              h6("Candidate genes: arbitrary selected genes, the five most frequently mutated genes, genes with significance in odds ratio of objective response rate or time on treatment."),
              h6("Objective response: CR or PR.")
      ),
      tabItem("Oddsratioondiseasecontrolrate",
              gt_output('figure_drug_8'),
              hr(),
              h6("Table. Univariable and multiple variable regression analysis were performed with gtsummary package for R. Factors selected in the multivariable analysis were selected by variable decreasing method based on the Akaike information criterion with stat package for R. Variables with variance inflation factor > 10 were excluded to eliminate multicollinearity."),
              h6("Candidate genes: arbitrary selected genes, the five most frequently mutated genes, genes with significance in odds ratio of objective response rate or time on treatment."),
              h6("Disease control: CR, PR, or SD.")
      ),
      tabItem("MutatedgenesandRECIST",
              htmlOutput("select_drug_ORR_table_var"),
              hr(),
              DT::dataTableOutput("table_outcome_3"),
              hr(),
              h6("OR: Objective response (CR or PR), ORR: Objective response rate, DC: Disease control (CR, PR, or SD), DCR: Disease control rate"),
              h6("95% confidence intervals were calculated using the Clopper-Pearson method.")
      ),
      tabItem("VolcanoPlotORR_1_AE",
              h6("Volcano plots for frequent regimens"),
              hr(),
              plotOutput('figure_volcano_1_AE',
                         height = "600px",
                         width = "800px"),
              hr(),
              htmlOutput("select_table_var_volcano_2_AE"),
              hr(),
              DT::dataTableOutput("table_volcano_AE"),
              hr(),
              h6("Figure. Patients with objective response data treated with the specified drugs were divided into two groups: those who obtained an Objective response and those who did not, and their odds ratios and p-values were calculated and a volcano plot was plotted."),
              h6("Odds ratio and p-value were calculated by logistic regression model concerning pathology."),
              h6("All regimens, and genes in which more than or equal to 8 of the treated patients had mutations were included in the analysis.")
      ),
      tabItem("Oddsratioonobjectiveresponserate_1_AE",
              gt_output('figure_drug_9_1'),
              hr(),
              h6("Table. Univariable and multiple variable regression analysis were performed with gtsummary package for R. Factors selected in the multivariable analysis were selected by variable decreasing method based on the Akaike information criterion with stat package for R. Variables with variance inflation factor > 10 were excluded to eliminate multicollinearity."),
              h6("Candidate genes: arbitrary selected genes, the five most frequently mutated genes, genes with significance in odds ratio of objective response rate or time on treatment."),
              h6("Objective response: CR or PR.")
      ),
      tabItem("Oddsratioonobjectiveresponserate_2_AE",
              gt_output('figure_drug_9_2'),
              hr(),
              h6("Table. Univariable and multiple variable regression analysis were performed with gtsummary package for R. Factors selected in the multivariable analysis were selected by variable decreasing method based on the Akaike information criterion with stat package for R. Variables with variance inflation factor > 10 were excluded to eliminate multicollinearity."),
              h6("Candidate genes: arbitrary selected genes, the five most frequently mutated genes, genes with significance in odds ratio of objective response rate or time on treatment."),
              h6("Objective response: CR or PR.")
      ),
      tabItem("Cumulative_incidence_AE",
              h6("Cumulative Incidence of Adverse Events"),
              hr(),
              htmlOutput("select_figure_CI_AE_var"),
              htmlOutput("select_AE_max_samples"),
              hr(),
              plotOutput('figure_CI_AE',
                         height = "800px",
                         width = "800px"),
              hr(),
              h6("Cumulative incidence of adverse events stratified by drug effectiveness, accounting for the competing risk of treatment completion. The observation period was defined as 90 days from treatment completion to capture both acute and delayed adverse events. Blue line: patients with objective response; Red line: patients without response. Shaded areas represent 95% confidence intervals. Statistical comparison was performed using Gray's test. The subdistribution hazard ratio with 95% confidence interval and p-value from the Fine-Gray competing risk regression model adjusted for age, sex, smoking status, diagnosis, and treatment line are displayed."),
              h6("No Response: PD/SD, Response: PR/CR")
      ),
      tabItem("Survival_drug",
              fluidRow(
                column(3,
                       htmlOutput("select_figure_survival_drug_var"),
                       htmlOutput("select_figure_survival_drug_var_2")
                ),
                column(3,
                       htmlOutput("select_drug_survival_gene"),
                       htmlOutput("select_figure_Drug_Bayes_max_samples")
                ),
                column(3,
                       actionButton("figure_survival_drug_reload", "Analyze")
                )
              ),
              plotOutput('figure_survival_drug',
                         height = "800px",
                         width = "800px")
      ),
      tabItem("Instruction",
              includeMarkdown("www/README.md")
      ),
      tabItem("Tips",
              includeMarkdown("www/Tips.md")
      )
    ),
    hr(),
    br()
  )
)
