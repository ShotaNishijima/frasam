
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
                         CppDir = system.file("executable",package="frasam"),...) {
  test <- try(frasyr::use_rvpa_tmb(TmbFile=TmbFile,CppDir=CppDir,...))
  if (class(test) == "DLLInfo") {
    return(TRUE)
  } else {
      return(FALSE)
    }
}
