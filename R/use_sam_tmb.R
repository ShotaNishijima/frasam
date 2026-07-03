
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
                        overwrite = FALSE,
                        auto_update = TRUE,
                        ...) {
  src <- file.path(CppDir, paste0(TmbFile, ".cpp"))
  dst <- paste0(TmbFile, ".cpp")

  if (isTRUE(auto_update) && file.exists(src) && file.exists(dst)) {
    overwrite <- file.info(src)$mtime > file.info(dst)$mtime
  }

  test <- try(
    frasyr::use_rvpa_tmb(
      TmbFile = TmbFile,
      CppDir = CppDir,
      overwrite = overwrite,
      ...
    )
  )

  if (inherits(test, "DLLInfo")) TRUE else FALSE
}
