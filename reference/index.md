# Package index

## All functions

- [`Calcu_Bmsy()`](Calcu_Bmsy.md) : Calculate the stock biomass at the
  maximum sustainable yield (Bmsy).
- [`Calcu_F0.1()`](Calcu_F0.1.md) : Calculate F0.1 value.
- [`Calcu_Fmax()`](Calcu_Fmax.md) : Calculate the fishing mortality that
  maximizes the yield per recruit.
- [`Calcu_Fmsy()`](Calcu_Fmsy.md) : Calculate the fishing mortality that
  maximizes the maximum sustainable yield (Fmsy) with fixed steepness
  parameters (h) in the Bevertoh-Holt stock-recruitment relationship.
- [`Calcu_N()`](Calcu_N.md) : Calculate the relative number at age at
  the equilibrium.
- [`Calcu_SBR()`](Calcu_SBR.md) : Calculate spawning biomass per recruit
  (SBR) with a given fishing mortality value.
- [`Calcu_SBmsy()`](Calcu_SBmsy.md) : Calculate the spawning stock
  biomass at the maximum sustainable yield (SSBmsy).
- [`Calcu_SPR_X()`](Calcu_SPR_X.md) : Calculate F%SPR: the fishing
  mortality that reduces the spawning biomass per recruit (SBR) to X %
  of those without fishing (SPR0).
- [`Calcu_YPR()`](Calcu_YPR.md) : Calculate yield per recruit with a
  given fishing mortality value.
- [`Est_SR()`](Est_SR.md) : Estimate the Beverton-Holt and hockey-stick
  stock recruitment relationships from data of spawning stock biomass
  (SBy) and recruits (Ry) with a give steepness (h) value.
- [`caa_plot()`](caa_plot.md) : Catch at
  ageの当てはまりについてプロットする関数
- [`calc_mase()`](calc_mase.md) :
  レトロの結果から各Indexに対する予測値を抽出して、Mean Absolute Scaled
  Errorを計算する関数
- [`calc_metrics()`](calc_metrics.md) : Calculate regression/forecast
  accuracy and bias metrics
- [`convert_sam_tibble()`](convert_sam_tibble.md) :
  SAMの結果オブジェクトをtibble形式に変換する関数
- [`divide_sigma()`](divide_sigma.md) :
  観測誤差とプロセス誤差をどこかの年齢間で分けて推定し直す関数
- [`do_jitter()`](do_jitter.md) : Do jitter analysis
- [`do_loo_index()`](do_loo_index.md) : Do leave-one-out index analysis
- [`do_osa_resid()`](do_osa_resid.md) : SAMのOSA residualを計算する関数
- [`est_fixed()`](est_fixed.md) : SAMの固定効果だけを推定する関数
- [`est_mixed()`](est_mixed.md) :
  TMB::MakeADFunに必要な引数から固定効果とランダム効果を推定する関数
- [`fit2PSdata()`](fit2PSdata.md) :
  生成された疑似データにVPA/SAMを推定し、Self-test/Cross-testを実行
- [`get_pm()`](get_pm.md) : utilities for extracting statistics
- [`get_predSR()`](get_predSR.md) : SAMで推定された再生産関係の予測値
- [`index_plot()`](index_plot.md) :
  Indexの当てはまりについてプロットする関数
- [`index_plot2()`](index_plot2.md) :
  Indexの当てはまりについてプロットする関数
- [`make_assess_result()`](make_assess_result.md) :
  SAMの結果を読み込んで主要なパラメータの予測値とSD、信頼区間などを出力する
- [`make_named_list()`](make_named_list.md) : Making a combined list
  with original list name(s)
- [`out_par()`](out_par.md) : SAMのFixed effect parametersの表を出力する
- [`out_sam()`](out_sam.md) : output sam object results
- [`plot_SR_simple()`](plot_SR_simple.md) :
  再生産関係についてプロットする関数
- [`plot_boosam()`](plot_boosam.md) :
  ブートストラップについてプロットする関数
- [`plot_hindcastCV()`](plot_hindcastCV.md) :
  レトロの結果を使って資源量指標値に対するhindcast cross
  validationをプロットする関数
- [`plot_osa_resid()`](plot_osa_resid.md) : OSA
  residualをプロットする関数
- [`plot_popsim()`](plot_popsim.md) : Popsimの結果ををプロットする関数
- [`plot_samvpa()`](plot_samvpa.md) : SAM or VPAの結果を描くグラフ
- [`popsim_vpasam()`](popsim_vpasam.md) :
  Popsimと同じ手法でSAMとVPAの疑似データを生成する
- [`retro_plot()`](retro_plot.md) :
  レトロスペクティブ解析の結果をプロットする
- [`retro_sam()`](retro_sam.md) : レトロスペクティブ解析を実施する
- [`rmvnorm_prec()`](rmvnorm_prec.md) : Simulate multivariate normal
  variables given a mean vector and precision matrix
- [`safe_do_call()`](safe_do_call.md) :
  do.callの引数が与えられていない場合にディフォルト値を与えて計算する関数
- [`sam()`](sam.md) : SAMによる資源計算を実施する
- [`samprofile()`](samprofile.md) :
  SAMで推定されたパラメータに対してプロファイル尤度を計算する
- [`select_sigma()`](select_sigma.md) :
  観測誤差やプロセス誤差のステップ形式のモデル選択（一つの変数について）
- [`select_sigma_grid()`](select_sigma_grid.md) :
  観測誤差やプロセス誤差のステップ形式のモデル選択（複数の変数について）
- [`sumup_popsim()`](sumup_popsim.md) : Self-test,
  Cross-testの結果をまとめるための関数
- [`update_sam()`](update_sam.md) : Update SAM result by new parameters
  of fixed and random effects
- [`use_sam_tmb()`](use_sam_tmb.md) :
  SAMでTMBで実行するためにcppファイルのコンパイル等をする関数
