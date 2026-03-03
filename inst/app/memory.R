get_system_memory <- function() {
  tryCatch({
    if(.Platform$OS.type == "unix") {
      # Linux: /proc/meminfo
      if(file.exists("/proc/meminfo")) {
        meminfo <- readLines("/proc/meminfo")
        total_line <- grep("^MemTotal:", meminfo, value = TRUE)
        available_line <- grep("^MemAvailable:", meminfo, value = TRUE)
        
        if(length(total_line) > 0 && length(available_line) > 0) {
          total_kb <- as.numeric(gsub("\\D", "", total_line))
          available_kb <- as.numeric(gsub("\\D", "", available_line))
          used_kb <- total_kb - available_kb
          
          return(list(
            total_mb = total_kb / 1024,
            used_mb = used_kb / 1024,
            free_mb = available_kb / 1024,
            success = TRUE
          ))
        }
      }
      
      # macOS: 実際の物理メモリを取得
      hw_memsize <- tryCatch({
        as.numeric(system("sysctl -n hw.memsize", intern = TRUE)) / (1024^2)
      }, error = function(e) NULL)
      
      if(!is.null(hw_memsize)) {
        # vm_statで使用量を取得
        vm_stat <- tryCatch({
          system("vm_stat", intern = TRUE)
        }, error = function(e) NULL)
        
        if(!is.null(vm_stat)) {
          page_size <- 4096
          free_pages <- as.numeric(gsub("\\D", "", grep("Pages free:", vm_stat, value = TRUE)))
          active_pages <- as.numeric(gsub("\\D", "", grep("Pages active:", vm_stat, value = TRUE)))
          inactive_pages <- as.numeric(gsub("\\D", "", grep("Pages inactive:", vm_stat, value = TRUE)))
          wired_pages <- as.numeric(gsub("\\D", "", grep("Pages wired down:", vm_stat, value = TRUE)))
          
          if(length(free_pages) > 0) {
            used_pages <- active_pages + inactive_pages + wired_pages
            used_mb <- (used_pages * page_size) / (1024^2)
            free_mb <- (free_pages * page_size) / (1024^2)
            
            return(list(
              total_mb = hw_memsize,
              used_mb = used_mb,
              free_mb = free_mb,
              success = TRUE
            ))
          }
        }
      }
    } else {
      # Windows: より正確な方法
      total_cmd <- 'wmic computersystem get TotalPhysicalMemory /value'
      available_cmd <- 'wmic OS get FreePhysicalMemory /value'
      
      total_output <- tryCatch({
        system(total_cmd, intern = TRUE)
      }, error = function(e) NULL)
      
      available_output <- tryCatch({
        system(available_cmd, intern = TRUE)
      }, error = function(e) NULL)
      
      if(!is.null(total_output) && !is.null(available_output)) {
        total_line <- grep("TotalPhysicalMemory=", total_output, value = TRUE)
        available_line <- grep("FreePhysicalMemory=", available_output, value = TRUE)
        
        if(length(total_line) > 0 && length(available_line) > 0) {
          total_bytes <- as.numeric(gsub(".*=", "", total_line))
          available_kb <- as.numeric(gsub(".*=", "", available_line))
          
          total_mb <- total_bytes / (1024^2)
          available_mb <- available_kb / 1024
          used_mb <- total_mb - available_mb
          
          return(list(
            total_mb = total_mb,
            used_mb = used_mb,
            free_mb = available_mb,
            success = TRUE
          ))
        }
      }
    }
    
    # フォールバック: 失敗
    return(list(
      total_mb = NA,
      used_mb = NA,
      free_mb = NA,
      success = FALSE
    ))
    
  }, error = function(e) {
    return(list(
      total_mb = NA,
      used_mb = NA,
      free_mb = NA,
      success = FALSE
    ))
  })
}

observe({
  invalidateLater(10000)  # 10秒毎
  
  # システムメモリ取得
  sys_mem <- get_system_memory()
  
  if(sys_mem$success) {
    # システムメモリを使用
    total_mb <- sys_mem$total_mb
    used_mb <- sys_mem$used_mb

    # 表示テキスト（システムメモリ）
    text <- paste0(round(used_mb/1024, 1), "GB/", round(total_mb/1024, 1), "GB")
    
  } else {
    # フォールバック: Rメモリ
    gc_info <- gc(verbose = FALSE)
    r_used_mb <- gc_info[1, 2] + gc_info[2, 2]
    r_limit_mb <- gc_info[2, 5]
    
    if(!is.na(r_limit_mb)) {
      text <- paste0(round(r_used_mb, 0), "MB")
    } else {
      text <- "N/A"
    }
  }
  
  # UI更新
  session$sendCustomMessage("updateMemoryWidget", list(
    text = text
  ))
})
