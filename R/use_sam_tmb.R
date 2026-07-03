
#' SAMでTMBで実行するためにcppファイルのコンパイル等をする関数
#'
#' @importFrom frasyr use_rvpa_tmb
#'
#' @param TmbFile Cppファイルの名前
#' @param CppDir Cppファイルが格納されているディレクトリ
#'
#' @encoding UTF-8
#'
#' @examples
#' \dontrun{
#' use_sam_tmb()
#' }
#'
#' @export

use_sam_tmb <- function(TmbFile = "sam2",
                        CppDir = system.file("executable", package = "frasam"),
                        RunDir = getwd(),
                        overwrite = FALSE,
                        compile = c("auto", "always", "never"),
                        auto_update = NULL,
                        ...) {
  compile <- match.arg(compile)

  if (!is.null(auto_update)) {
    warning("'auto_update' is deprecated. Please use 'compile = \"auto\"' instead.",
            call. = FALSE)
  }

  test <- try(
    frasyr::use_rvpa_tmb(
      TmbFile = TmbFile,
      CppDir = CppDir,
      RunDir = RunDir,
      overwrite = overwrite,
      compile = compile,
      ...
    )
  )

  if (inherits(test, "DLLInfo")) TRUE else FALSE
}
