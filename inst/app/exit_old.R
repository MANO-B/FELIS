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

# =========================================================
# 1. ログアウトボタン（確認のみ行い、R側へ通知）
# =========================================================
observeEvent(input$logout_button, {
  message('[R] Logout button clicked')

  # 確認ダイアログのみJSで表示
  shinyjs::runjs('
    if (confirm("ログアウトしますか？")) {
      console.log("[JS] Logout confirmed");
      // R側の input$confirm_logout を true に設定してイベント発火
      Shiny.setInputValue("confirm_logout", true, {priority: "event"});
    }
  ')
})

# =========================================================
# 2. 確認後の処理（do_logoutを呼び出し）
# =========================================================
observeEvent(input$confirm_logout, {
  message('[R] confirm_logout fired - calling do_logout()')
  do_logout(session, reason = 'explicit_button')
})

# =========================================================
# 3. ログアウト実行関数（クリーンアップ ＋ 画面遷移）
# =========================================================
do_logout <- function(session, reason = 'unknown') {
  message('[R] ========== do_logout() called ==========')
  message('[R] Reason: ', reason)

  if (is.null(session) || isTRUE(session$userData$logged_out)) {
    return(invisible(NULL))
  }

  # フラグを立てる
  session$userData$logged_out <- TRUE

  # -------------------------------------------------------
  # A. クリーンアップ処理 (メモリ解放)
  # -------------------------------------------------------
  # ※ onSessionEnded と重複しても問題ないよう tryCatch で囲むか、
  #    共通関数化するとより良いですが、ここでは直接実行します。
  message('[R] Performing explicit cleanup before redirect...')

  large_objects <- c(
    'Data_case_raw', 'Data_case', 'Data_drug_raw',
    'Data_drug_raw_rename', 'Data_report_raw', 'Data_report',
    'OUTPUT_DATA', 'Data_cluster_ID', 'tmp_post', 'analysis_env'
  )

  # グローバル環境の変数を削除
  for (obj_name in large_objects) {
    if (exists(obj_name, envir = .GlobalEnv)) {
      rm(list = obj_name, envir = .GlobalEnv)
    }
  }

  # ガベージコレクション
  gc()

  # -------------------------------------------------------
  # B. 画面遷移 (リダイレクト)
  # -------------------------------------------------------
  # 条件判定（サーバー環境かつCCATフラグがTRUEの場合のみ遷移など）
  do_redirect <- isTRUE(ENV_ == "server" && isTRUE(CCAT_FLAG))

  if (do_redirect) {
    message('[R] Redirecting to logout URL...')
    shinyjs::runjs('
      const logoutUrl = window.location.origin + "/@@/logout";
      console.log("[JS] Redirecting to:", logoutUrl);
      // 同じタブで移動
      window.location.href = logoutUrl;
    ')
  } else {
    # 条件に合致しない場合（ローカル環境など）は、単にアプリを閉じる等の処理
    message('[R] Not redirecting (Condition not met). Closing session.')
    shinyjs::runjs('alert("ログアウトしました");')
    session$close() # 遷移しない場合のみ明示的に閉じる
  }

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
