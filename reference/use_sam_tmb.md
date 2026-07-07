# SAMでTMBで実行するためにcppファイルのコンパイル等をする関数

SAMでTMBで実行するためにcppファイルのコンパイル等をする関数

## Usage

``` r
use_sam_tmb(
  TmbFile = "sam2",
  CppDir = system.file("executable", package = "frasam"),
  RunDir = getwd(),
  overwrite = FALSE,
  compile = c("auto", "always", "never"),
  auto_update = NULL
)
```

## Arguments

- TmbFile:

  Cppファイルの名前

- CppDir:

  Cppファイルが格納されているディレクトリ

- RunDir:

  CppファイルとDLLを配置するディレクトリ

- overwrite:

  RunDirのCppファイルを上書きするかどうか

- compile:

  DLLをコンパイルするかどうか。"auto"ではCppファイルがDLLより新しい場合にコンパイルする

- auto_update:

  互換性のために残している引数。compile = "auto"を使用してください

## Examples

``` r
if (FALSE) { # \dontrun{
use_sam_tmb()
} # }
```
