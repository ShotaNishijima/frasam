
#' do.callの引数が与えられていない場合にディフォルト値を与えて計算する関数
#'
#' 余計な引数があった場合、その関数を削除して実行する
#'
#' @param func 関数名
#' @param args_list 引数名
#' @encoding UTF-8
#'
#' @examples
#' example_func <- function(a, b = 2, c = 3) {
#'   a + b + c
#' }
#'
#' args <- list(a = 1, d = 4)  # `d` is an extra argument
#' safe_do_call(example_func, args)  # Returns 6 (1 + 2 + 3)
#'
#' @export
#'
safe_do_call <- function(func, args_list, default_args = list()) {
  # 関数の引数を取得
  func_args <- formals(func)

  # 使用可能な引数のみをフィルタリング
  valid_args <- intersect(names(args_list), names(func_args))
  filtered_args <- args_list[valid_args]

  # リストにデフォルト引数を追加 (引数がない場合)
  for (arg_name in names(func_args)) {
    if (!arg_name %in% names(filtered_args)) {
      if (arg_name %in% names(default_args)) {
        filtered_args[[arg_name]] <- default_args[[arg_name]]
      } else if (!is.null(func_args[[arg_name]])) {
        filtered_args[[arg_name]] <- eval(func_args[[arg_name]])
      }
    }
  }

  # 関数を実行
  do.call(func, filtered_args)
}

