# do.callの引数が与えられていない場合にディフォルト値を与えて計算する関数

余計な引数があった場合、その関数を削除して実行する

## Usage

``` r
safe_do_call(func, args_list, default_args = list())
```

## Arguments

- func:

  関数名

- args_list:

  引数名

## Examples

``` r
example_func <- function(a, b = 2, c = 3) {
  a + b + c
}

args <- list(a = 1, d = 4)  # `d` is an extra argument
safe_do_call(example_func, args)  # Returns 6 (1 + 2 + 3)
#> [1] 6
```
