#' Get status icon
#'
#' Internal helper to standardize log icons.
#'
#' @param type Character. One of "step", "running", "completed", "info", "warning", "error".
#' @return A string with the icon.
#' @keywords internal
#' @noRd

get_icon <- function(type) {
  switch(type,
    "step"      = "▶️ ", # 补齐 1 格
    "running"   = "⏳",
    "completed" = "✅",
    "info"      = "❕",
    "warning"   = "⚠️ ", # 补齐 1 格
    "error"     = "❌"
  )
}
