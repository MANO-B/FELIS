# =========================================================
# Session init (logout/cleanup architecture)
# - Logout: explicit button only
# - Browser close / reload: NO logout (cleanup only)
# - Hang / OOM: avoid forced logout in onSessionEnded
# - Large objects: store in session$userData$cache_env (NOT .GlobalEnv)
# =========================================================

# -----------------------------
# 0) Session state initialize
# -----------------------------
session$userData$logged_out <- FALSE
session$userData$session_start_time <- Sys.time()

# App専用キャッシュ環境（ここに巨大オブジェクトを入れる）
session$userData$cache_env <- new.env(parent = emptyenv())

message("[R] Session initialized at: ", session$userData$session_start_time)

# -----------------------------
# 1) Utility: safe cleanup
# -----------------------------
cleanup_session <- function(session) {
  message("[R] Cleanup started")

  tryCatch({
    e <- session$userData$cache_env
    if (!is.null(e) && is.environment(e)) {
      objs <- ls(envir = e, all.names = TRUE)
      if (length(objs) > 0) {
        rm(list = objs, envir = e)
        message("[R] Removed from cache_env: ", paste(objs, collapse = ", "))
      } else {
        message("[R] cache_env is already empty")
      }
    } else {
      message("[R] cache_env not found (or not environment)")
    }
  }, error = function(e) {
    message("[R] Cleanup error: ", conditionMessage(e))
  })

  # 最低限のGC
  invisible(gc())

  message("[R] Cleanup completed")
  invisible(NULL)
}

# -----------------------------
# 2) Logout executor (explicit only)
# -----------------------------
do_logout <- function(session, reason = "explicit_button") {
  message("[R] ========== do_logout() called ==========")
  message("[R] Reason: ", reason)

  if (is.null(session) || isTRUE(session$userData$logged_out)) {
    message("[R] do_logout skipped (session null or already logged out)")
    return(invisible(NULL))
  }

  # logoutフラグ
  session$userData$logged_out <- TRUE

  # cleanup (safe, idempotent)
  cleanup_session(session)

  # redirect condition (your existing logic)
  do_redirect <- isTRUE(ENV_ == "server" && isTRUE(CCAT_FLAG))

  if (do_redirect) {
    message("[R] Redirecting to logout URL...")
    shinyjs::runjs('
      const logoutUrl = window.location.origin + "/@@/logout";
      console.log("[JS] Redirecting to:", logoutUrl);
      window.location.href = logoutUrl;
    ')
  } else {
    message("[R] Not redirecting (Condition not met). Closing session.")
    shinyjs::runjs('alert("ログアウトしました");')
    session$close()
  }

  message("[R] ========== do_logout() finished ==========")
  invisible(NULL)
}

# -----------------------------
# 3) Logout button flow
#    - confirm only in JS
#    - notify R via input$confirm_logout
# -----------------------------
observeEvent(input$logout_button, {
  message("[R] Logout button clicked")

  shinyjs::runjs('
    if (confirm("ログアウトしますか？")) {
      console.log("[JS] Logout confirmed");
      // 値は true だと二度目以降発火しないケースがあるので、タイムスタンプ推奨
      Shiny.setInputValue("confirm_logout", Date.now(), {priority: "event"});
    } else {
      console.log("[JS] Logout cancelled");
    }
  ')
})

observeEvent(input$confirm_logout, {
  message("[R] confirm_logout fired - calling do_logout()")
  do_logout(session, reason = "explicit_button")
})

# -----------------------------
# 4) onSessionEnded
#    - IMPORTANT: DO NOT logout here
#    - only cleanup (safe)
# -----------------------------
session$onSessionEnded(function() {
  message("[R] ========== session$onSessionEnded() called ==========")
  message("[R] This session started at: ", session$userData$session_start_time)
  message("[R] logged_out flag: ", session$userData$logged_out)

  # ブラウザclose / reload / hang などは区別しない
  # いずれにせよ logout はしない（ここ重要）
  cleanup_session(session)

  message("[R] ========== session$onSessionEnded() finished ==========")
})
