# グローバル環境に最新セッション開始時刻を保存（全セッションで共有）
if (!exists(".last_session_start_time", envir = .GlobalEnv)) {
  .GlobalEnv$.last_session_start_time <- NULL
}

# セッション初期化
session$userData$logged_out <- FALSE
session$userData$session_start_time <- Sys.time()

# グローバル変数を更新
.GlobalEnv$.last_session_start_time <- session$userData$session_start_time

message('[R] Session initialized at: ', session$userData$session_start_time)

# JavaScript
shinyjs::runjs("
  console.log('[JS] Script initialized');

  $(document).ready(function() {
    console.log('[JS] Document ready');

    // navigation typeを送信（デバッグ用）
    var navType = 'unknown';
    if (performance.navigation) {
      switch(performance.navigation.type) {
        case 0:
          navType = 'navigate';
          break;
        case 1:
          navType = 'reload';
          break;
        case 2:
          navType = 'back_forward';
          break;
      }
    }
    console.log('[JS] Navigation type:', navType);
  });
")

# ログアウトボタン
observeEvent(input$logout_button, {
  message('[R] Logout button clicked')

  shinyjs::runjs('
    console.log("[JS] Logout button event");
    if (confirm("Want to log out?")) {
      console.log("[JS] Logout confirmed");
      Shiny.setInputValue("confirm_logout", true, {priority: "event"});
    }
  ')
})

observeEvent(input$confirm_logout, {
  message('[R] confirm_logout fired - calling do_logout()')
  do_logout(session, reason = 'explicit_button')
})

# ログアウト処理
do_logout <- function(session, reason = 'unknown') {
  message('[R] ========== do_logout() called ==========')
  message('[R] Reason: ', reason)

  if (is.null(session) || isTRUE(session$userData$logged_out)) {
    message('[R] Already logged out - returning')
    return(invisible(NULL))
  }

  session$userData$logged_out <- TRUE
  message('[R] logged_out set to TRUE')

  # セッションが有効な場合のみJavaScriptを実行
  if (reason != 'session_ended') {
    tryCatch({
      message('[R] Showing JavaScript alert')
      shinyjs::runjs('
      const logoutUrl = window.location.origin + "/@@/logout";
      console.log("[JS] Sending logout request to:", logoutUrl);

      fetch(logoutUrl, {
        method: "GET",
        headers: { "Content-Type": "application/json" },
        credentials: "same-origin",
        body: JSON.stringify({ reason: "logout" })
      })
        .then(resp => {
          console.log("[JS] Logout request sent. Status:", resp.status);
          alert("Session finished. Please close this tab. Logging out via " + logoutUrl);
          setTimeout(() => window.close(), 1500);
        })
        .catch(err => {
          console.error("[JS] Logout request failed:", err);
          alert("Logout failed: " + err);
        });
    ')
    }, error = function(e) {
      message("[R] Could not run JavaScript (session may be closed): ", e$message)
    })
  }

  later::later(function() {
    message('[R] Closing session - Reason: ', reason)
    try(session$close(), silent = TRUE)
    if (interactive()) stopApp()
  }, delay = 1)

  message('[R] ========== do_logout() finished ==========')
}

# セッション終了時の処理
session$onSessionEnded(function() {
  message('[R] ========== session$onSessionEnded() called ==========')
  message('[R] This session started at: ', session$userData$session_start_time)

  if (is.null(session)) {
    return(NULL)
  }

  message('[R] logged_out flag: ', session$userData$logged_out)

  # 0.8秒待って、新しいセッションが開始されたかチェック
  later::later(function() {
    message('[R] --- Checking for new session (after 0.8s delay) ---')

    is_reload <- FALSE

    if (!is.null(.GlobalEnv$.last_session_start_time)) {
      message('[R] Global last session start: ', .GlobalEnv$.last_session_start_time)
      message('[R] This session started: ', session$userData$session_start_time)

      # グローバルの最新セッション開始時刻が、このセッションより後なら、リロード
      if (.GlobalEnv$.last_session_start_time > session$userData$session_start_time) {
        time_diff <- as.numeric(difftime(.GlobalEnv$.last_session_start_time,
                                         session$userData$session_start_time, units = "secs"))
        message('[R] New session started ', round(time_diff, 2), ' seconds after this session')
        is_reload <- TRUE
      } else {
        message('[R] No new session detected - this is a window close')
      }
    }

    # リロードの場合はログアウトをスキップ
    if (is_reload) {
      message('[R] Reload detected - skipping logout')

      # クリーンアップのみ実行
      message('[R] Cleanup started')

      large_objects <- c(
        'Data_case_raw', 'Data_case', 'Data_drug_raw',
        'Data_drug_raw_rename', 'Data_report_raw', 'Data_report',
        'OUTPUT_DATA', 'Data_cluster_ID', 'tmp_post', 'analysis_env'
      )

      for (obj_name in large_objects) {
        if (exists(obj_name, envir = .GlobalEnv)) {
          rm(list = obj_name, envir = .GlobalEnv)
          message('[R] Removed: ', obj_name)
        }
      }

      all_objects <- ls(envir = .GlobalEnv)
      for (obj_name in all_objects) {
        tryCatch({
          obj <- get(obj_name, envir = .GlobalEnv)
          size_mb <- as.numeric(object.size(obj)) / (1024^2)
          if (size_mb > 10) {
            message('[R] Removing large object: ', obj_name, ' (', round(size_mb, 1), 'MB)')
            rm(list = obj_name, envir = .GlobalEnv)
          }
        }, error = function(e) {})
      }

      session$reactlog(FALSE)
      gc()

      message('[R] Cleanup completed')
      message('[R] ========== session$onSessionEnded() finished (reload skipped) ==========')
      return(NULL)
    }

    # 通常のクローズまたはログアウト
    if (!isTRUE(session$userData$logged_out)) {
      message('[R] Window/tab closed detected - calling do_logout()')
      do_logout(session, reason = 'window_closed')
    }

    message('[R] Cleanup started')

    large_objects <- c(
      'Data_case_raw', 'Data_case', 'Data_drug_raw',
      'Data_drug_raw_rename', 'Data_report_raw', 'Data_report',
      'OUTPUT_DATA', 'Data_cluster_ID', 'tmp_post', 'analysis_env'
    )

    for (obj_name in large_objects) {
      if (exists(obj_name, envir = .GlobalEnv)) {
        rm(list = obj_name, envir = .GlobalEnv)
        message('[R] Removed: ', obj_name)
      }
    }

    all_objects <- ls(envir = .GlobalEnv)
    for (obj_name in all_objects) {
      tryCatch({
        obj <- get(obj_name, envir = .GlobalEnv)
        size_mb <- as.numeric(object.size(obj)) / (1024^2)
        if (size_mb > 10) {
          message('[R] Removing large object: ', obj_name, ' (', round(size_mb, 1), 'MB)')
          rm(list = obj_name, envir = .GlobalEnv)
        }
      }, error = function(e) {})
    }

    session$reactlog(FALSE)
    gc()

    message('[R] Cleanup completed')
    message('[R] ========== session$onSessionEnded() finished ==========')
  }, delay = 0.8)
})
