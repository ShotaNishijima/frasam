#' 条件付き負の対数尤度を成分別に集計する
#'
#' SAMの推定結果から、資源尾数過程、漁獲死亡係数過程、および観測系列ごとの
#' 条件付き負の対数尤度を取得します。
#'
#' @param samres [sam()] が返す `sam` オブジェクト。
#' @param index_name 資源量指標の表示名。`NULL` の場合は
#'   `"Index_1"`, `"Index_2"`, ... を使用します。漁獲量系列を除く
#'   観測系列数と同じ長さの文字ベクトルを指定してください。
#'
#' @details
#' この関数が返す値は、推定されたランダム効果に条件付けた負の対数尤度です。
#' TMBのLaplace近似によってランダム効果を積分した周辺負の対数尤度では
#' ありません。また、初期状態の密度、体重・成熟モデル、ランダム効果
#' `logB` の密度、およびペナルティ項は戻り値に含まれません。
#'
#' C++モデルが返す観測ごとの `ans_obs` をfleetごとに合計します。
#' fleet 1を年齢別漁獲尾数、fleet 2以降を資源量指標として扱います。
#'
#' @return 次の列を持つ[tibble][tibble::tibble]。
#' \describe{
#'   \item{type}{尤度成分または観測系列を表すfactor。}
#'   \item{nll}{条件付き負の対数尤度。}
#' }
#'
#' @examples
#' \dontrun{
#' data("samres_example", package = "frasam")
#' get_cond_nll(samres)
#' }
#'
#' @export
get_cond_nll <- function(samres, index_name = NULL) {
  if (!inherits(samres, "sam")) {
    stop("'samres' must be an object of class 'sam'.", call. = FALSE)
  }
  if (is.null(samres$obj) || !is.function(samres$obj$report)) {
    stop("'samres$obj$report' is not available.", call. = FALSE)
  }
  if (is.null(samres$data$obs)) {
    stop("'samres$data$obs' is not available.", call. = FALSE)
  }

  report <- samres$obj$report()
  required <- c("ans_n", "ans_f", "ans_obs")
  missing <- required[vapply(required, function(x) is.null(report[[x]]), logical(1))]
  if (length(missing) > 0L) {
    stop(
      "The TMB report is missing: ", paste(missing, collapse = ", "), ". ",
      "Recompile the model from a compatible sam2.cpp.",
      call. = FALSE
    )
  }

  obs <- as.data.frame(samres$data$obs)
  if (!"fleet" %in% names(obs)) {
    stop("'samres$data$obs' must contain a 'fleet' column.", call. = FALSE)
  }

  ans_obs <- as.numeric(report$ans_obs)
  if (length(ans_obs) != nrow(obs)) {
    stop(
      "The length of 'ans_obs' does not match nrow(samres$data$obs).",
      call. = FALSE
    )
  }
  obs$nll <- ans_obs

  obs_nll <- dplyr::summarise(
    dplyr::group_by(obs, fleet),
    nll = sum(nll),
    .groups = "drop"
  )
  obs_nll <- dplyr::arrange(obs_nll, fleet)

  if (!1 %in% obs_nll$fleet) {
    stop("Fleet 1 (catch at age) is missing from the observations.", call. = FALSE)
  }

  index_fleets <- obs_nll$fleet[obs_nll$fleet != 1]
  nindex <- length(index_fleets)
  if (is.null(index_name)) {
    index_name <- paste0("Index_", seq_len(nindex))
  } else if (!is.character(index_name) || length(index_name) != nindex || anyNA(index_name)) {
    stop(
      "'index_name' must be a character vector with one name per index fleet.",
      call. = FALSE
    )
  }

  obs_nll$type <- NA_character_
  obs_nll$type[obs_nll$fleet == 1] <- "Catch_at_age"
  obs_nll$type[match(index_fleets, obs_nll$fleet)] <- index_name

  process_nll <- tibble::tibble(
    type = c("Process_N", "Process_F"),
    nll = c(as.numeric(report$ans_n), as.numeric(report$ans_f))
  )

  result <- dplyr::bind_rows(
    process_nll,
    dplyr::select(obs_nll, type, nll)
  )
  result$type <- factor(result$type, levels = unique(result$type))
  result
}
