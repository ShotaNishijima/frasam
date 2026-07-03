
#' SAMでTMBで実行するためにcppファイルのコンパイル等をする関数
#'
#' @param TmbFile Cppファイルの名前
#' @param CppDir Cppファイルが格納されているディレクトリ
#' @param RunDir CppファイルとDLLを配置するディレクトリ
#' @param overwrite RunDirのCppファイルを上書きするかどうか
#' @param compile DLLをコンパイルするかどうか。"auto"ではCppファイルがDLLより新しい場合にコンパイルする
#' @param auto_update 互換性のために残している引数。compile = "auto"を使用してください
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
                        auto_update = NULL) {
  if (!requireNamespace("TMB", quietly = TRUE)) {
    stop("Please install TMB package!", call. = FALSE)
  }

  compile <- match.arg(compile)

  if (!is.null(auto_update)) {
    warning("'auto_update' is deprecated. Please use 'compile = \"auto\"' instead.",
            call. = FALSE)
  }

  cpp_src <- file.path(CppDir, paste0(TmbFile, ".cpp"))
  cpp_dst <- file.path(RunDir, paste0(TmbFile, ".cpp"))
  dll_dst <- file.path(RunDir, paste0(TmbFile, .Platform$dynlib.ext))

  if (!file.exists(cpp_src)) {
    stop("Cpp file not found: ", cpp_src, call. = FALSE)
  }

  file.copy(
    from = cpp_src,
    to = cpp_dst,
    overwrite = overwrite
  )

  if (!file.exists(cpp_dst)) {
    stop("Failed to copy Cpp file to: ", cpp_dst, call. = FALSE)
  }

  needs_compile <- switch(
    compile,
    always = TRUE,
    never = FALSE,
    auto = !file.exists(dll_dst) ||
      file.info(cpp_dst)$mtime > file.info(dll_dst)$mtime
  )

  owd <- setwd(RunDir)
  on.exit(setwd(owd), add = TRUE)

  if (needs_compile) {
    if (file.exists(dll_dst)) {
      try(dyn.unload(dll_dst), silent = TRUE)
    }

    TMB::compile(paste0(TmbFile, ".cpp"))
  }

  if (file.exists(dll_dst)) {
    try(dyn.unload(dll_dst), silent = TRUE)
  }

  dyn.load(TMB::dynlib(TmbFile))
  TRUE
}
