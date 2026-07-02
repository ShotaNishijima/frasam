#' SAMによる資源計算を実施する
#'
#' @param dat samに使用するdataでrvpaと同じフォーマットで利用可能
#' @param rec.age 加入年齢 (default: 0)
#' @param alpha 最高年齢-1歳へのFに対する最高年齢のFの比
#' @param upper 推定パラメータの上限値。固定効果の数のLengthを持つ必要あり。NULL（デフォルト）の場合Inf
#' @param lower 推定パラメータの下限値。固定効果の数のLengthを持つ必要あり。NULL（デフォルト）の場合-Inf
#' @param abund Indexの種類。用いるIndexの長さのベクトル。"B": 資源量、"SSB"：親魚資源量、"N"：尾数、"Bs"：資源量×fleetごとの選択率、"Bf"：資源量×選択率．合計する年齢の幅については\code{min.age}と\code{max.age}
#' @param catch_prop abund="Bs"のとき、対象とするfleetのcatch at ageの全体に対する比率.
#' @param min.age Indexの最低年齢 (\code{frasyr::vpa()}と同じで最小の年齢を0とする) 用いるIndexの長さのベクトル
#' @param max.age Indexの最高年齢 (\code{frasyr::vpa()}と同じで最小の年齢を0とする) 用いるIndexの長さのベクトル
#' @param b.est Indexと資源量の間の非線形関係を考慮しない（FALSE, デフォルト）、考慮する(TRUE)
#' @param b.fix b.est=TRUEの場合、非線形パラメータbを推定するか（NA）、固定するか（固定する値） ??
#' @param index.b.key Indexのbの制約 ??  どうやって使う？？
#' @param sel.def 選択率の定義。"max"（デフォルト）の場合、最大年齢を１とする。
#' @param use.index この仕様は設定ミスを引き起こしやすいので廃止しました。\code{sam()}を使用するまえに、Indexのデータを使用するもののみにsubsetするようにしてください
#' @param varC CAAの観測誤差
#' @param varN ??
#' @param varF ??
#' @param varNfix ??
#' @param rho.mode Fの多変量Random walkの非対角成分の相関係数のタイプ。１ならすべての年齢間で相関係数1, 0なら完全にランダム, 2は任意の年齢間で共通のrhoを推定、3は年齢i,jの相関を\eqn{rho^|i-j|}で推定
#' @param no_est 推定しない(TRUE)、パラメータ推定する（FALSE, デフォルト）
#' @param getJointPrecision JointPrecision matrixを計算しない（FALSE, デフォルト）、計算する（TRUE)
#' @param loopnum 最適化を繰り返す回数．デフォルトは2
#' @param index.key Index間の観測誤差sigmaの制約。NULLの場合はIndexごとに別々、rep(0, length(abund))の場合は全Indexで共通。
#' @param est.method Deprecated. Use \code{index.key} instead. If \code{"ls"} is supplied and \code{index.key = NULL}, all index observation-error sigmas are shared.
#' @param plus.group プラスグループを考慮する（TRUE, デフォルト）、考慮しない（FALSE）
#' @param SR 再生産関係："RW"(Random walk), "BH", "RI", "HS", "Mesnil", or "Const"
#' @param AR 再生産関係の残差の自己相関パラメータを推定する（１）、推定しない（0, デフォルト） ?? これでいい？
#' @param b_random 再生産パラメータbをランダム効果として推定するかどうか??
#' @param b_range 再生産パラメータbの値の範囲
#' @param lambda ??
#' @param lambda_Mesnil SRを"Mesnil"にした場合のラムダの値??
#' @param q.init パラメータqの初期値。NULLの場合にはexp(-5)が用いられる。
#' @param sdFsta.init Fのランダムウォークのσの初期値。NULLの場合（デフォルト）にはlog_sigma=-0.693147が用いられる。初期値を与える場合には、ノーマルスケールでの値を与える。
#' @param sdLogN.init Nのプロセス誤差のσの初期値。NULLの場合（デフォルト）にはlog_sigma=0.35が用いられる（←コード上はこうなっているがこれで良い？）。初期値を与える場合には、ノーマルスケールでの値を与える。
#' @param sdLogObs.init Indexの観察誤差の初期値。NULLの場合（デフォルト）にはlog_sigma=-0.356675が用いられる。初期値を与える場合には、ノーマルスケールでの値を与える。
#' @param rho.init rhoの初期値。NULLの場合（デフォルト）には0が用いられる。
#' @param a.init 再生産関係パラメータaの初期値。NULLの場合（デフォルト）、log_a=8 (SR="Const")またはlog_a=-4（SR="Const"以外）
#' @param b.init 再生産関係パラメータbの初期値。NULLの場合（デフォルト）、log_b=7 (SR="HS", "Mesnil", "BHS") またはlog_b=-8（それ以外）
#' @param ref.year 管理基準値を計算するときの参照年。最終年から何年分さかのぼるか。デフォルトは1:5（最新年からさかのぼって５年分）。
#' @param bias.correct 固定効果パラメータ？？のバイアスを補正する（TRUE：デフォルト）、補正しない（FALSE）
#' @param bias.correct.sd ランダム効果のSDパラメータ？？のバイアスを補正する（TRUE）、補正しない（FALSE：デフォルト）
#' @param get.random.vcov ランダム効果の分散共分散行列を推定する（TRUE、時間かかります）、推定しない（FALSE：デフォルト）
#' @param silent MakeADfunのときの標準出力あり（TRUE: デフォルト）、なし（FALSE）
#' @param model_wm weightとmaturityの成長をモデリングするかどうか
#' @importFrom glmmTMB glmmTMB
#' @param last.catch.zero 最終年の漁獲量がない場合。デフォルトはFALSE（最終年の漁獲量が利用できて用いる）
#' @param cpp.file.name 推定に用いるcppファイル。デフォルトは最新版の"sam2"
#' @param scale 資源量のスケーリングファクター。資源量はscaleで割った値となる
#' @param scale_number 尾数のスケーリングファクター。尾数はscaleで割った値となる
#' @param remove.Fprocess.year Years to remove from the F process model.
#' @param RW.Forder Order of the random-walk process for fishing mortality.
#' @param map.add Additional TMB map settings.
#' @param p0.list Initial parameter list.
#' @param gamma Gamma parameter setting.
#' @param FreeADFun \code{TMB::FreeADFun}を使う場合、TRUEにする。See \code{?TMB::FreeADFun}.
#' @param add_random Additional random effects.
#' @param tmbdata TMB data list.
#' @param map TMB parameter map.
#' @param w0_factor Factor setting for initial weight.
#' @param weight_factor Factor setting for weight-at-age.
#' @param family_w Error distribution for weight observations.
#' @param maturity_factor Factor setting for maturity-at-age.
#' @param weight_weight Weight assigned to weight-at-age likelihood components.
#' @param maturity_weight Weight assigned to maturity-at-age likelihood components.
#' @param g_fix Fixed value setting for growth parameter g.
#' @param CV_w_fix Fixed value setting for weight coefficient of variation.
#' @param w_link Link function for weight model.
#' @param sep_omicron Whether to separate omicron parameters.
#' @param growth_regime Growth regime setting.
#'
#' @export
#'

sam <- function(dat,
                last.catch.zero = FALSE,
                cpp.file.name = "sam2",
                tmb.run = FALSE,
                abund = c("B"),
                rec.age = 0,
                min.age = rep(0,length(abund)),
                max.age = rep(0,length(abund)),
                plus.group = TRUE,
                alpha = 1.0,
                b.est = FALSE,
                b.fix = rep(NA,length(abund)),
                varC = 0,
                varN = 0,
                varF = 0,
                varN.fix = NULL,
                est.method = NULL,
                SR = "BH", #"RW": random walk, "RI": Ricker
                AR = 0,
                rho.mode = 2,
                # MA = 0,
                q.init = NULL,
                sdFsta.init = NULL,
                sdLogN.init = NULL,
                sdLogObs.init = NULL,
                rho.init = NULL,
                a.init = NULL,
                b.init = NULL,
                ref.year = 1:5,
                bias.correct = TRUE,
                bias.correct.sd = FALSE,
                get.random.vcov = FALSE,
                silent = FALSE,
                remove.Fprocess.year = NULL,
                RW.Forder= 0,
                map.add = NULL,
                p0.list = NULL,
                scale=1000,
                gamma=10,
                sel.def="max",
                use.index = NULL,
                upper = NULL,
                lower = NULL,
                index.key = NULL,
                index.b.key = NULL,
                b_random = FALSE,
                b_range = NULL,
                lambda = 0,
                FreeADFun = FALSE,
                add_random = NULL,
                lambda_Mesnil = 0,
                tmbdata = NULL,
                map = NULL,
                model_wm = c(FALSE,FALSE),
                w0_factor =c("none","SSB")[1],
                weight_factor = c("none","N_total")[1],
                family_w = c("lognormal","gamma")[2],
                maturity_factor = c("none","cohort_plus")[1],
                scale_number=1000,
                weight_weight = NULL,
                maturity_weight = NULL,
                g_fix = NULL,
                CV_w_fix = NULL,
                w_link = "log",
                sep_omicron=TRUE,
                growth_regime = NULL,
                catch_prop = NULL,
                no_est=FALSE,
                getJointPrecision = FALSE,
                loopnum = 2,
                obj_overwrite=NULL,
                ignore.parm.uncertainty = FALSE
){

  argname <- ls()
  arglist <- lapply(argname,function(xx) eval(parse(text=xx)))
  names(arglist) <- argname

  if (!is.null(est.method)) {
    warning(
      "'est.method' is deprecated. Please use 'index.key' instead. ",
      "For the old est.method = 'ls' behavior, use index.key = rep(0, length(abund)).",
      call. = FALSE
    )

    if (identical(est.method, "ls") && is.null(index.key)) {
      index.key <- rep(0, length(abund))
    }
  }
  arglist$index.key <- index.key

  index.age <- min.age

  if (is.null(tmbdata)) {
    if(isTRUE(b.est)){
      if (!is.null(index.b.key)) {
        if (length(unique(index.b.key)) != length(b.fix)){
          stop("'length(unique(index.b.key)) == length(b.fix)' must be satisfied")
        }
      }
    }

    # check and convert data (why? see: test-vignette-data.R)
    if(!is.data.frame(dat$caa)) dat$caa <- as.data.frame(dat$caa)
    if(!is.data.frame(dat$maa)) dat$maa <- as.data.frame(dat$maa)
    if(!is.data.frame(dat$M  )) dat$M   <- as.data.frame(dat$M)
    if(!is.data.frame(dat$waa)) dat$waa <- as.data.frame(dat$waa)
    if(!is.data.frame(dat$index)) dat$index <- as.data.frame(dat$index)
    if (!is.null(dat$waa.catch) && !is.data.frame(dat$waa.catch)) {
      dat$waa.catch <- as.data.frame(dat$waa.catch)
    }
    arglist$dat <- dat

    caa <- dat$caa
    if (last.catch.zero) caa[ncol(caa)] <- NULL

    maa <- dat$maa
    waa <- dat$waa
    M <- dat$M
    index <- dat$index

    if (!is.null(use.index)) {
      warning(
        "'use.index' is deprecated. Please subset dat$index, abund, min.age, max.age, index.key, and related arguments
  before calling sam().",
        call. = FALSE
      )
    }

    for(i in 1:length(abund)) {
      if (abund[i]=="SSB") {
        index.age[i] <- rec.age
        max.age[i] <- rec.age + nrow(waa)-1
      }
    }

    obs <- suppressWarnings(
      cbind(expand.grid(as.numeric(rownames(caa)),as.numeric(colnames(caa))),1,unlist(caa))[,c(2,3,1,4)]
    )
    obs <- cbind(obs,obs[,3])

    nindex <- nrow(index)
    for (i in 1:nindex){
      index1 <- index[i,]
      index2 <- cbind(i+1,as.numeric(names(index1)[!is.na(index1)]),index.age[i],c(index1[!is.na(index1)]),max.age[i])[,c(2,1,3,4,5)]
      obs <- rbind(as.matrix(obs), index2)
    }

    colnames(obs) <- c("year","fleet","age","obs","maxage")
    rownames(obs) <- NULL

    obs <- as.matrix(obs)
    obs[,"age"] <- obs[,"age"]+rec.age
    obs[,"maxage"] <- obs[,"maxage"]+rec.age

    if (tmb.run) {
      # library(TMB)
      compile(paste(cpp.file.name, ".cpp", sep = ""))
      dyn.load(TMB::dynlib(cpp.file.name))
    }

    data <- list()

    data$obs <- obs
    data$noFleets <- max(obs[,2])
    data$fleetTypes <- data$sampleTimes <- numeric(nindex+1)
    assertthat::assert_that(all(abund %in% c("B","SSB","N","Bs","Bf")))
    for (i in 1:nindex) {
      if (is.null(abund[i])) data$fleetTypes[i+1] <- data$freetTypes[i]
      if (abund[i] == "B") data$fleetTypes[i+1] <- 2
      if (abund[i] == "SSB") data$fleetTypes[i+1] <- 3
      if (abund[i] == "N") data$fleetTypes[i+1] <- 4
      if (abund[i] == "Bs") data$fleetTypes[i+1] <- 6
      if (abund[i] == "Bf") data$fleetTypes[i+1] <- 7
      # if (!(abund[i] %in% c("B","SSB","N","Bs","Bf"))) {
      #   stop("abund code not recognized")
      # }
    }

    data$noYears <- ncol(waa)
    data$years <- as.numeric(colnames(waa))
    # if (retro.years>0) {
    #   for(i in 1:retro.years) data$years <- c(data$years,max(data$years)+1)
    # }
    data$iy <- data$obs[,1]-min(data$years) #year since the first year
    data$nobs <- nrow(data$obs)
    # colnames(maa) <- colnames(waa) <- colnames(M) <- range(data$obs[,1])[1]:((range(data$obs[,1])[2])+retro.years)
    data$propMat <- data$propMat2 <- as.matrix(t(maa))
    data$stockMeanWeight <- data$catchMeanWeight <- as.matrix(t(waa))
    data$natMor <- as.matrix(t(M))
    data$minAge <- as.matrix(rec.age)
    data$maxAge <- as.matrix(rec.age+nrow(caa)-1) #cppファイルには使わないデータ
    data$maxAgePlusGroup <- ifelse(isTRUE(plus.group), 1, 0)
    data$rhoMode <- rho.mode
    ncol1 <- as.numeric(data$maxAge-data$minAge+1)
    data$landFrac <- data$disMeanWeight <- data$landMeanWeight <- data$propF <- data$propM <- matrix(0, nrow=data$noYears, ncol=ncol1)

    basemat <- matrix(-1, ncol = ncol1, nrow = max(data$obs[,2]))
    data$keyLogFsta <- data$keyLogQ <- data$keyLogB <- data$keyVarObs <- data$keyVarF <- data$keyVarLogN <- basemat

    data$keyLogFsta[1,] <- c(0:(data$maxAge-1-rec.age),(data$maxAge-1-rec.age))
    for (i in 1:nindex) data$keyLogQ[i+1,index.age[i]+1] <- i-1
    if (isTRUE(b.est)) {
      if(is.null(index.b.key)) {
        for (i in 1:nindex) {
          data$keyLogB[i+1,index.age[i]+1] <- i-1
        }
      }else{
        for (i in 1:nindex) {
          data$keyLogB[i+1,index.age[i]+1] <- index.b.key[i]-min(index.b.key)
        }
      }
    }

    data$keyVarObs[1,] <- varC
    if (is.null(index.key)) {
      for (i in 1:nindex) data$keyVarObs[i+1,index.age[i]+1] <- max(varC)+i
    } else {
      if (length(index.key) != nindex) {
        stop("'index.key' must have the same length as the number of indices")
      }
      for (i in 1:nindex) {
        data$keyVarObs[i+1,index.age[i]+1] <- max(varC)+index.key[i]-min(index.key)+1
      }
    }

    data$keyVarF[1,] <- varF
    if(rev(data$keyVarF[1,])[1]!=rev(data$keyVarF[1,])[2]) {
      stop("'varF' in the maximum age class (A) should be the same as that in A-1")
    }

    data$keyVarLogN[1,] <- varN

    if (SR == "BH") SR.mode <- 2 #Beverton-Holt
    if (SR == "RI") SR.mode <- 1 #Ricker
    if (SR == "RW") SR.mode <- 0 #Random walk
    if (SR == "HS") SR.mode <- 3 #Hockey-stick
    if (SR == "Mesnil") SR.mode <- 4 #Hockey-stick
    if (SR == "Const") SR.mode <- 5 # Constant R0(=a)
    if (SR == "Prop") SR.mode <- 6
    if (SR == "BHS") SR.mode <- 7
    if (SR == "MR") SR.mode <- 8
    data$stockRecruitmentModelCode <- matrix(SR.mode)

    data$scale <- scale
    data$gamma <- gamma
    if(sel.def=="max"){
      data$sel_def <- 0
    }
    if(sel.def=="mean"){
      data$sel_def <- 1
    }
    if(sel.def=="maxage"){
      data$sel_def <- 2
    }
    if(isTRUE(b_random)){
      data$b_random <- 1
    }else{
      data$b_random <- 0
    }

    data$nlogF = max(data$keyLogFsta)+1
    data$nlogN = as.numeric(data$maxAge-data$minAge)+1
    # data$AR <- AR
    data$alpha <- alpha
    # data$MA <- MA

    data$Fprocess_weight = rep(1,ncol(waa))
    if (!is.null(remove.Fprocess.year)) {
      data$Fprocess_weight[which(as.numeric(colnames(dat$waa)) %in% remove.Fprocess.year)] <- 0
    }
    data$F_RW_order <- RW.Forder
    data$lambda <- lambda
    data$lambda_Mesnil <- lambda_Mesnil

    data$model_weight_maturity <- 2*(as.numeric(model_wm)-0.5)
    data$scale_number <- scale_number
    if (family_w == "lognormal") {
      dist_wobs = 0
    } else {
      if (family_w == "gamma") {
        dist_wobs = 1
      } else {
        stop("'family_w' not recognized")
      }
    }
    data$dist_wobs = dist_wobs
    if(is.null(weight_weight)) {
      weight_weight = matrix(1,ncol=ncol(dat$waa),nrow=nrow(dat$waa))
      if (isTRUE(last.catch.zero)) weight_weight[,ncol(dat$waa)] <- 0
    }
    if(!isTRUE(all( dim(weight_weight) == dim(dat$waa)))) {
      stop("'weight_weight' has a wrong dimension!")
    }
    data$weight_weight <- weight_weight
    if(is.null(maturity_weight)) {
      maturity_weight = matrix(1,ncol=ncol(dat$maa),nrow=nrow(dat$maa))
      if (isTRUE(last.catch.zero)) maturity_weight[,ncol(dat$maa)] <- 0
    }
    if(!isTRUE(all( dim(maturity_weight) == dim(dat$maa)))) {
      stop("'maturity_weight' has a wrong dimension!")
    }
    data$maturity_weight <- maturity_weight
    if(is.null(g_fix)) {
      g <- 1
      g_fix <- c()
      for(i in 1:nrow(dat$maa)) {
        if(sd(as.numeric(dat$maa[i,]),na.rm=TRUE)==0) {
          g_fix <- c(g_fix,mean(as.numeric(dat$maa[i,]),na.rm=FALSE))
        } else {
          g_fix <- c(g_fix,-g)
          g <- g + 1
        }
      }
    }
    data$g_fix <- g_fix
    if(w_link=="id"){
      alpha_w_link = 0
    } else {
      if (w_link=="log") {
        alpha_w_link = 1
      } else {
        stop("'alpha_w_link' not recognized")
      }
    }
    data$alpha_w_link <- alpha_w_link
    if (is.null(growth_regime)) {
      growth_regime = rep(0,ncol(dat$waa))
    } else {
      if (length(growth_regime) != ncol(dat$waa)) {
        stop("The length of 'growth regime' must be the number of years")
        growth_regime <- growth_regime - min(growth_regime)
      }
    }
    data$gr <- growth_regime
    if(!is.null(catch_prop)){
      if(!isTRUE(all(dim(catch_prop)==c(dim(dat$waa),max(data$obs[,2]))))) {
        stop("Dimension of 'catch_prop' is incorrect!")
      } else {
        data$catch_prop4index <- catch_prop
      }
    } else {
      catch_prop <- array(1,dim=c(dim(dat$waa),max(data$obs[,2])))
      data$catch_prop4index <- catch_prop
    }
    assertthat::assert_that(all(dim(data$catch_prop4index) == c(dim(dat$waa),max(data$obs[,2]))))
    data$logobs = log(obs[,4])
  } else {
    message("'dat' and related arguments are ignored when using 'tmbdata'")
    data = tmbdata
    ncol1 <- as.numeric(data$maxAge-data$minAge+1)
    nindex = max(data$obs[,"fleet"]-1)
    SR.mode = as.numeric(data$stockRecruitmentModelCode)
  }

  if(isTRUE(model_wm[1]) && is.null(p0.list)) {
    waa_dat = expand.grid(Age=as.numeric(rownames(dat$waa)),
                          Year=as.numeric(colnames(dat$waa))) %>%
      dplyr::mutate(Weight=as.numeric(unlist(dat$waa)),
                    Maturity=as.numeric(unlist(dat$maa)))
    A = max(waa_dat$Age)
    waa_prev = sapply(1:nrow(waa_dat), function(i) {
      if (waa_dat$Age[i]==0 | waa_dat$Year[i]==min(waa_dat$Year)) {
        value <- NA
      } else {
        if (waa_dat$Age[i]<A) {
          value <- as.numeric(
            dat$waa[waa_dat$Age[i],as.character(waa_dat$Year[i]-1)])
        } else {
          value1 <- as.numeric(
            dat$waa[waa_dat$Age[i],as.character(waa_dat$Year[i]-1)])
          value2 <- as.numeric(
            dat$waa[waa_dat$Age[i]+1,as.character(waa_dat$Year[i]-1)]) # plus group
          n1 <- as.numeric(
            dat$caa[waa_dat$Age[i],as.character(waa_dat$Year[i]-1)])
          n2 <- as.numeric(
            dat$caa[waa_dat$Age[i]+1,as.character(waa_dat$Year[i]-1)]) # plus group
          value <- (value1*n1+value2*n2)/(n1+n2)
        }
      }
      value
    })
    maa_prev = sapply(1:nrow(waa_dat), function(i) {
      if (waa_dat$Age[i]==0 | waa_dat$Year[i]==min(waa_dat$Year)) {
        value <- NA
      } else {
        value <- as.numeric(
          dat$maa[waa_dat$Age[i],as.character(waa_dat$Year[i]-1)])
      }
      value
    })
    delta_mat = function(Maturity,Maturity_prev) (Maturity-Maturity_prev)/(1-Maturity_prev)

    waa_dat = waa_dat %>%
      mutate(Weight_prev=waa_prev,
             Maturity_prev=maa_prev) %>%
      mutate(YearClass=Year-Age)
    waa_dat2 = waa_dat %>% filter(Age>min(Age) & Year > min(Year))

    waa_dat %>% filter(Age %in% c(0,1))
    maa_dat = waa_dat %>% filter(Age %in% (which(g_fix<0)-1)) %>%
      mutate(y=delta_mat(Maturity,Maturity_prev)) %>%
      na.omit()

    modw = glm(Weight~1+Weight_prev,data=waa_dat2,family=Gamma("identity"))
    nr = max(growth_regime)+1 #number of growth regime
    alpha_w = matrix(rep(as.numeric(coef(modw)[1]),nr),nrow=1)
    if (w_link=="log") alpha_w = log(alpha_w)
    rho_w = matrix(rep(as.numeric(coef(modw)[2]),nr),nrow=1)
    logCV_w <- log(sqrt(summary(modw)$dispersion))
    beta_w0 = waa_dat %>% dplyr::filter(Age==min(Age)) %>% dplyr::pull(.,Weight) %>% mean %>% log
    beta_w0 = matrix(rep(beta_w0,nr),nrow=1)
    if (weight_factor!="none") {
      alpha_w = rbind(alpha_w,0)
    }
    if (w0_factor!="none") {
    beta_w0 = rbind(beta_w0,0)
    }
    # tmp = waa_dat %>% filter(Age==min(Age)) %>% pull(.,Weight) %>% log %>% sd
    # logCV_w <- c(log(tmp),logCV_w)
    if(!is.null(CV_w_fix)) logCV_w = rep_len(CV_w_fix,length(logCV_w))
  } else {
    nr=1
    alpha_w = matrix(0)
    rho_w = matrix(0)
    logCV_w <- 0
    beta_w0 = matrix(0)
  }

  if (isTRUE(model_wm[2]) && is.null(p0.list)) {
    mod_mat0 = glmmTMB(y~factor(Age) + 0,data=maa_dat,family=ordbeta)
    alpha_g = mod_mat0$fit$par[names(mod_mat0$fit$par)=="beta"] %>% as.numeric()
    alpha_g = matrix(rep(alpha_g,nr),ncol=nr)
    psi = as.numeric(mod_mat0$fit$par[names(mod_mat0$fit$par)=="psi"])
    logdisp = as.numeric(mod_mat0$fit$par[names(mod_mat0$fit$par)=="betad"])[1]
  } else {
    alpha_g <- matrix(0)
    psi = rep(0,2)
    logdisp = 0
  }
  logwaa = as.matrix(log(dat$waa))
  beta_g = matrix(0,nrow=nrow(alpha_g),ncol=ncol(alpha_g))

  U_init = rbind(
    matrix(5,nrow=ncol1,ncol=data$noYears),
    matrix(-1,nrow=max(data$keyLogFsta)+1,ncol=data$noYears)
    )

  check_init_length <- function(x, name, expected_length, positive = FALSE, open_unit = FALSE) {
    if (is.null(x)) return(invisible(TRUE))
    if (!is.numeric(x)) {
      stop(sprintf("'%s' must be numeric.", name), call. = FALSE)
    }
    if (length(x) != expected_length) {
      stop(
        sprintf("'%s' must have length %d, but length %d was supplied.", name, expected_length, length(x)),
        call. = FALSE
      )
    }
    if (any(!is.finite(x))) {
      stop(sprintf("'%s' must contain only finite values.", name), call. = FALSE)
    }
    if (isTRUE(positive) && any(x <= 0)) {
      stop(sprintf("'%s' must contain only positive values.", name), call. = FALSE)
    }
    if (isTRUE(open_unit) && any(x <= 0 | x >= 1)) {
      stop(sprintf("'%s' must be between 0 and 1.", name), call. = FALSE)
    }
    invisible(TRUE)
  }

  check_init_length(q.init, "q.init", nindex, positive = TRUE)
  check_init_length(sdFsta.init, "sdFsta.init", max(data$keyVarF)+1, positive = TRUE)
  check_init_length(sdLogN.init, "sdLogN.init", max(data$keyVarLogN)+1, positive = TRUE)
  check_init_length(sdLogObs.init, "sdLogObs.init", max(data$keyVarObs)+1, positive = TRUE)
  check_init_length(a.init, "a.init", 1, positive = TRUE)
  check_init_length(b.init, "b.init", 1, positive = TRUE)
  check_init_length(rho.init, "rho.init", 1, open_unit = TRUE)

  validate_p0_list <- function(p0.list, expected, skip_structure = character()) {
    if (!is.list(p0.list)) {
      stop("'p0.list' must be a list.", call. = FALSE)
    }

    missing_names <- setdiff(names(expected), names(p0.list))
    extra_names <- setdiff(names(p0.list), names(expected))
    if (length(missing_names) > 0 || length(extra_names) > 0) {
      msg <- c()
      if (length(missing_names) > 0) msg <- c(msg, paste0("missing: ", paste(missing_names, collapse = ", ")))
      if (length(extra_names) > 0) msg <- c(msg, paste0("unexpected: ", paste(extra_names, collapse = ", ")))
      stop(
        sprintf("'p0.list' does not match the current model parameter names (%s).", paste(msg, collapse = "; ")),
        call. = FALSE
      )
    }

    format_dim <- function(x) {
      if (is.null(dim(x))) "NULL" else paste(dim(x), collapse = " x ")
    }

    problems <- character()
    for (nm in names(expected)) {
      x <- p0.list[[nm]]
      ref <- expected[[nm]]
      if (!is.numeric(x)) {
        problems <- c(problems, sprintf("%s must be numeric", nm))
        next
      }
      if (any(!is.finite(x))) {
        problems <- c(problems, sprintf("%s contains non-finite values", nm))
        next
      }
      if (nm %in% skip_structure) next
      if (!identical(dim(x), dim(ref))) {
        problems <- c(problems, sprintf("%s has dim %s; expected %s", nm, format_dim(x), format_dim(ref)))
        next
      }
      if (length(x) != length(ref)) {
        problems <- c(problems, sprintf("%s has length %d; expected %d", nm, length(x), length(ref)))
        next
      }
    }

    if (length(problems) > 0) {
      stop(
        sprintf("'p0.list' does not match the current model parameter structure:\n- %s", paste(problems, collapse = "\n- ")),
        call. = FALSE
      )
    }
    invisible(TRUE)
  }

  logSdLogN_init <- if (is.null(sdLogN.init)) rep(0.356675,max(data$keyVarLogN)+1) else log(sdLogN.init)

  if(SR == "Const") logSdLogN_init[1] <- log(2)
  if (!is.null(varN.fix)) {
    if (length(varN.fix) != max(data$keyVarLogN) + 1) {
      stop("'varN.fix' must have length max(data$keyVarLogN) + 1.", call. = FALSE)
    }
    if (any(!is.na(varN.fix) & varN.fix <= 0)) {
      stop("'varN.fix' must contain positive values or NA.", call. = FALSE)
    }
    map_logSdLogN = 0:max(data$keyVarLogN)
    for(i in 1:(max(data$keyVarLogN)+1)) {
      if(!is.na(varN.fix[i])){
        logSdLogN_init[i] <- 0.5*log(varN.fix[i])
        map_logSdLogN[i] <- NA
      }
    }
  }



  ## p0.listとvarN.fixが両方使われてる場合に、varN.fixが効くように修正する (2026/07/02)

  if (!is.null(p0.list) && !is.null(varN.fix) && any(!is.na(varN.fix))) {
    if (length(p0.list$logSdLogN) != length(logSdLogN_init)) {
      stop("'p0.list$logSdLogN' and 'varN.fix' have incompatible lengths.", call. = FALSE)
    }

    fixed_varN <- !is.na(varN.fix)
    p0.list$logSdLogN[fixed_varN] <- logSdLogN_init[fixed_varN]

    message("Initial values for fixed 'logSdLogN' were overwritten by 'varN.fix'.")
  }

  # b.fixを使うときにp0.listも使っていると、b.fixが効かなくなるので追加&修正　2026/7/2
  if (is.null(index.b.key)) {
    logB_init <- sapply(
      seq_len(nindex),
      function(i) ifelse(is.na(b.fix[i]), 0, log(b.fix[i]))
    )
  } else {
    logB_init <- sapply(
      seq_along(b.fix),
      function(i) ifelse(is.na(b.fix[i]), 0, log(b.fix[i]))
    )
  }

  if (!is.null(p0.list) && isTRUE(b.est) && !is.null(b.fix) && any(!is.na(b.fix))) {
    if (length(p0.list$logB) != length(logB_init)) {
      stop("'p0.list$logB' and 'b.fix' have incompatible lengths.", call. = FALSE)
    }

    fixed_b <- !is.na(b.fix)
    p0.list$logB[fixed_b] <- logB_init[fixed_b]

    message("Initial values for fixed 'logB' were overwritten by 'b.fix'.")
  }

  default_params <- list(
      logQ      = if (is.null(q.init)) rep(-5,nindex) else log(q.init),
      logB      = logB_init,
      logSdLogFsta = if (is.null(sdFsta.init)) rep(-0.693147,max(data$keyVarF)+1) else log(sdFsta.init),
      logSdLogN    = logSdLogN_init,
      logSdLogObs  = if (is.null(sdLogObs.init)) rep(-0.356675,max(data$keyVarObs)+1) else log(sdLogObs.init),
      rec_loga     = if (is.null(a.init)) { if (SR == "Const") {
        8
      } else{
        -4
      }
        } else {
          log(a.init)
          },
      rec_logb     = if (is.null(b.init)) {if (SR %in% c("HS","Mesnil","BHS","MR")) 7 else -8} else {log(b.init)},
      logit_rho  = if (is.null(rho.init)) 0 else log(rho.init/(1-rho.init)),
      # logScale     = numeric(data$noScaledYears),
      # logScaleSSB  = if(any(data$fleetTypes %in% c(3,4))) {numeric(0)} else {numeric(0)},
      # logPowSSB    = if(any(data$fleetTypes == 4))        {numeric(0)} else {numeric(0)},
      # logSdSSB     = if(any(data$fleetTypes %in% c(3,4))) {numeric(0)} else {numeric(0)},
      U = U_init,
      trans_phi1 = 0,
      logSD_b = log(0.3),
      logwaa = logwaa,
      beta_w0 = beta_w0,
      alpha_w = alpha_w,
      rho_w = rho_w,
      omicron = if(isTRUE(sep_omicron)) c(log(0.1),log(0.1)) else log(0.1),
      logCV_w = logCV_w,
      alpha_g = alpha_g,
      psi = psi,
      logdisp = logdisp,
      beta_g = beta_g,
      rec_logk = log(1)
    )

  if (is.null(p0.list)){
    params <- default_params
  } else {
    skip_structure <- character()
    if (isTRUE(model_wm[1])) skip_structure <- c(skip_structure, "beta_w0", "alpha_w", "rho_w", "omicron", "logCV_w")
    if (isTRUE(model_wm[2])) skip_structure <- c(skip_structure, "alpha_g", "psi", "logdisp", "beta_g")
    validate_p0_list(p0.list, default_params, skip_structure = skip_structure)
    params <- p0.list[names(default_params)]
  }

  if(is.null(map)) {
    map <- list()
    if (!isTRUE(b_random)) {
      map$logSD_b <- factor(NA)
      if (!isTRUE(b.est)) map$logB <- factor(rep(NA,nindex)) else {
        if (!is.null(b.fix)) {
          map$logB <- 0:max(data$keyLogB)

          if (is.null(index.b.key)) {
            for (i in seq_len(nindex)) {
              if (!is.na(b.fix[i])) map$logB[i] <- NA
            }
          } else {
            for (i in seq_along(b.fix)) {
              if (!is.na(b.fix[i])) map$logB[i] <- NA
            }
          }

          map$logB <- factor(map$logB)
        }
      }
    }

    if (SR=="RW") {
      map$rec_loga <- map$rec_logb <- factor(NA)
      map$trans_phi1 <- factor(NA)
    } else {
      if (AR==0) map$trans_phi1 <- factor(NA)
      if (SR=="Const" | SR=="Prop") map$rec_logb <- factor(NA)
      # if (AR==1) map$phi2 <- factor(NA)
    }

    if (rho.mode==0) map$logit_rho <- factor(NA)
    if (rho.mode==1) map$logit_rho <- factor(NA)

    if(!is.null(varN.fix)) map$logSdLogN <- factor(map_logSdLogN)

    if(!isTRUE(model_wm[1])) {
      map$logwaa <- factor(matrix(NA,ncol=ncol(dat$waa),nrow=nrow(dat$waa)))
      map$alpha_w <- rep(factor(NA),prod(dim(params$alpha_w)))
      map$beta_w0 <- rep(factor(NA),prod(dim(params$beta_w0)))
      map$rho_w <- rep(factor(NA),prod(dim(params$rho_w)))
      map$omicron <- rep(factor(NA),length(params$omicron))
      map$logCV_w <- rep(factor(NA),length(params$logCV_w))
    }
    if(!isTRUE(model_wm[2])) {
      map$alpha_g = rep(factor(NA),prod(dim(params$alpha_g)))
      map$psi = c(factor(NA),factor(NA))
      map$logdisp = factor(NA)
      map$beta_g = rep(factor(NA),prod(dim(params$beta_g)))
    } else {
      if (maturity_factor=="none") {
        map$beta_g = rep(factor(NA),prod(dim(params$beta_g)))
      }
    }
    if (!(SR %in% c("BHS","MR"))) map$rec_logk <- factor(NA)

    if(!is.null(CV_w_fix)) map$logCV_w = rep(factor(NA),length(logCV_w))

    if (!is.null(map.add)) {
      tmp = make_named_list(map.add)
      for(i in 1:length(tmp$map.add)) {
        map[[names(tmp$map.add)[i]]] <- map.add[[i]]
      }
    }
  }

  # stop("Tentative Stop!!")
  random = c("U","logwaa")

  if(isTRUE(b_random)) {
    random = c(random,"logB")
    # obj <- TMB::MakeADFun(data, params, map = map, random=c("U","logB"), DLL=cpp.file.name,silent=silent)
  }
  if(!is.null(add_random)) {
    random = c(random,add_random)
  }

    obj <- TMB::MakeADFun(data, params, map = map, random=random, DLL=cpp.file.name,silent=silent)
    obj$fn(obj$par)

    if(!is.null(obj_overwrite)) obj <- obj_overwrite

  if(isTRUE(FreeADFun)) {
    TMB::FreeADFun(obj)
  }

  if (is.null(lower)) lower <- obj$par*0-Inf
  if (is.null(upper)) upper <- obj$par*0+Inf

  if (!is.null(b_range) & "rec_logb" %in% names(obj$par)) {
    lower["rec_logb"] <- log(b_range[1])
    upper["rec_logb"] <- log(b_range[2])
    if (params$rec_logb < log(b_range[1]) | params$rec_logb > log(b_range[2])) {
      stop("Set the initial value of 'rec_logb' within 'b_range'")
    }
  }

  if(!isTRUE(no_est)) {
    nlminb.control = list(eval.max = 1e4,
                          iter.max = 1e4,
                          trace = 0)

    if(!is.null(obj_overwrite)) {
      opt <- list()
      opt$par <- obj$par
      opt$convergence <- 0
      opt$objective <- NA
    } else {
      # inital optimization
      opt <- nlminb(obj$par, obj$fn, obj$gr, lower=lower, upper=upper,control=nlminb.control)

      # Re-run to further decrease final gradient (https://github.com/kaskr/TMB_contrib_R/blob/master/TMBhelper/R/fit_tmb.R)
      for( i in seq(2,loopnum,length=max(0,loopnum-1)) ){
        # Temp = parameter_estimates[c('iterations','evaluations')]
        opt2 = nlminb( start=obj$par, objective=obj$fn, gradient=obj$gr, control=nlminb.control, lower=lower, upper=upper )
        if(opt2$objective<=opt$objective) {
          opt <- opt2
        } else {
          break
        }
      }
    }

    if (opt$convergence!=0) warning("May not converge")
    rep <- TMB::sdreport(obj,bias.correct = bias.correct,bias.correct.control = list(sd=bias.correct.sd), getReportCovariance=get.random.vcov,getJointPrecision=getJointPrecision
                         ,ignore.parm.uncertainty=ignore.parm.uncertainty
                         )
    if (max(rep$gradient.fixed)>1e-2) warning("Large maximum gradient component")

    # stop("Finish Optimization")

    SR.name <- "RW"
    is.SR <- ifelse(SR.mode==0, FALSE, TRUE)
    # par_list = obj$env$parList(opt$par)
    par_list = obj$env$parList() #最適化しないように修正

    BioRefPt <- function(data, rep, opt, ref.year, is.SR, obj){
      if (is.null(rep$unbiased)) {
        NF <- matrix(rep$par.random[names(rep$par.random)=="U"],nrow=max(data$keyLogFsta)+1 + data$maxAge-data$minAge+1,ncol=data$noYears)
        logN <- NF[1:data$nlogN,]
        logF <- NF[(data$nlogN+1):(data$nlogN+data$nlogF),]
        logF <- logF[data$keyLogFsta[1,]+1,]
        logF[nrow(logF),] <- log(data$alpha)+logF[nrow(logF),]
        tmp = rep$value[names(rep$value)=="stockMeanWeight_true"]
        waa_est = t(matrix(tmp,nrow=ncol(logN)))
        # waa_obs = dat$waa
        colnames(logN) <- colnames(logF) <- colnames(waa_est) <- data$years
        rownames(logN) <- rownames(logF) <- rownames(waa_est) <- data$minAge:data$maxAge
        naa <- exp(logN)
        faa <- exp(logF)
      } else {
        # logN <- matrix(rep$unbiased$value[names(rep$unbiased$value)=="logN"],ncol=data$noYears)
        # logF <- matrix(rep$unbiased$value[names(rep$unbiased$value)=="logF"],ncol=data$noYears)
        naa <- matrix(rep$unbiased$value[names(rep$unbiased$value)=="exp_logN"],ncol=data$noYears)
        faa <- matrix(rep$unbiased$value[names(rep$unbiased$value)=="exp_logF"],ncol=data$noYears)
        faa <- faa[data$keyLogFsta[1,]+1,]
        faa[nrow(faa),] <- data$alpha*faa[nrow(faa),]
        tmp = rep$unbiased$value[names(rep$value)=="stockMeanWeight_true"]
        waa_est = t(matrix(tmp,nrow=ncol(naa)))
        colnames(naa) <- colnames(faa) <- colnames(waa_est) <- data$years
        rownames(naa) <- rownames(faa) <- rownames(waa_est) <- data$minAge:data$maxAge
        logN <- log(naa)
        logF <- log(faa)
      }
      waa_obs = dat$waa

      ref.year1 <- data$noYears-ref.year+1

      # naa <- exp(logN)
      # faa <- exp(logF)
      baa <- naa*waa_est
      ssb <- baa*t(data$propMat)
      if(sel.def=="max"){
        saa <- sweep(faa,2,apply(faa,2,max),FUN="/")
      }
      if(sel.def=="mean"){
        saa <- sweep(faa,2,apply(faa,2,mean),FUN="/")
      }
      if(sel.def=="maxage"){
        saa <- sweep(faa,2,faa[nrow(faa),],FUN="/")
      }
      zaa <- faa+t(data$natMor)
      caa <- faa/zaa*naa*(1-exp(-zaa))

      SR.rec <- rec.par <- BRP0 <- BRPmsy <- NULL
      if (is.SR){
        # SR

        SR.name <- case_when(as.numeric(data$stockRecruitmentModelCode)==0 ~ "RW",
                             as.numeric(data$stockRecruitmentModelCode)==1 ~ "RI",
                             as.numeric(data$stockRecruitmentModelCode)==2 ~ "BH",
                             as.numeric(data$stockRecruitmentModelCode)==3 ~ "HS",
                             as.numeric(data$stockRecruitmentModelCode)==4 ~ "Mesnil",
                             as.numeric(data$stockRecruitmentModelCode)==5 ~ "Const",
                             as.numeric(data$stockRecruitmentModelCode)==7 ~ "BHS",
                             as.numeric(data$stockRecruitmentModelCode)==8 ~ "MR")
        #
        # a <- exp(rep$par.fixed[names(rep$par.fixed)=="rec_loga"])
        # b <- exp(rep$par.fixed[names(rep$par.fixed)=="rec_logb"])
        a <- exp(obj$env$parList()[["rec_loga"]])
        b <- exp(obj$env$parList()[["rec_logb"]])
        if (SR.mode==5) b <- NA
        if (SR.mode==6) b <- 0
        if (SR.mode %in% c(7,8)) k <- exp(obj$env$parList()[["rec_logk"]])
        #
        rec.par <- c(a,b)
        names(rec.par) <- c("a","b")
        if(SR.mode==7) {
          rec.par <- c(rec.par,k)
          names(rec.par)[3] <- "k"
        }
        #
      }

      loglik <- -opt$objective
      aic <- 2*opt$objective+2*length(opt$par)

      # q1 <- exp(as.numeric(rep$par.fixed[names(rep$par.fixed) == "logQ"]))
      q1 <- exp(as.numeric(obj$env$parList()[["logQ"]]))
      if (isTRUE(b_random)){
        b1 <- exp(as.numeric(rep$par.random[names(rep$par.random) == "logB"]))
      }else{
        if(!isTRUE(b.est)) b1 <- rep(1,nindex) else {
          b1 <- exp(obj$env$parList()[["logB"]])
          # if(is.null(b.fix)) b1 <- exp(as.numeric(rep$par.fixed[names(rep$par.fixed) == "logB"]))
          # if(!is.null(b.fix)) {
          #   b1 <- b.fix
          #   b1[is.na(b1)] <- exp(as.numeric(rep$par.fixed[names(rep$par.fixed) == "logB"]))
          # }
          if(!is.null(index.b.key)) {
            # b1 <- exp(obj$env$parList()[["logB"]])
            b1 <- b1[index.b.key-min(index.b.key)+1]
          }
        }
      }
      # sigma0 <- exp(as.numeric(rep$par.fixed[names(rep$par.fixed) == "logSdLogObs"]))
      sigma0 <- exp(as.numeric(par_list[["logSdLogObs"]]))
      sigma1 <- sigma0[1:(max(data$keyVarObs[1,])+1)]
      sigma1 <- sapply(1:ncol1, function(i) sigma1[data$keyVarObs[1,i]+1])
      sigma2 <- sigma0[(max(data$keyVarObs[1,])+2):length(sigma0)]
      if(!is.null(index.key)) sigma2 <- sigma2[index.key-min(index.key)+1]
      sigma3 <- exp(as.numeric(par_list[["logSdLogFsta"]]))
      sigma3 <- sapply(1:ncol1, function(i) sigma3[data$keyVarF[1,i]+1])
      # sigma4 <- exp(obj$env$parList(opt$par)$logSdLogN)
      sigma4 <- exp(obj$env$parList()[["logSdLogN"]]) #修正2025/11/12
      # sigma4 <- exp(as.numeric(rep$par.fixed[names(rep$par.fixed) == "logSdLogN"]))
      sigma4 <- sapply(1:ncol1, function(i) sigma4[data$keyVarLogN[1,i]+1])
      rho1 <- ifelse(as.numeric(data$rhoMode) > 1, 1/(1+exp(-as.numeric(rep$par.fixed[names(rep$par.fixed) == "logit_rho"]))),as.numeric(data$rhoMode))
      phi<-numeric(1)
      if (AR>0) phi[1] = rep$value[names(rep$value)=="phi1"][1]
      pred.index <- data.frame(matrix(0,nrow=nindex,ncol=nrow(data$stockMeanWeight)))
      colnames(pred.index) <- rownames(data$stockMeanWeight)
      for (i in 1:nrow(pred.index)){
        if (abund[i]=="N") {
          pred.index[i,] <- 0
          for (j in index.age[i]:max.age[i]) {
            pred.index[i,] <- pred.index[i,]+as.numeric(naa[j+1,])
          }
        }
        if (abund[i]=="B") {
          pred.index[i,] <- 0
          for (j in index.age[i]:max.age[i]) {
            pred.index[i,] <- pred.index[i,]+as.numeric(baa[j+1,])/scale
          }
        }
        if (abund[i]=="SSB") pred.index[i,] <- as.numeric(colSums(ssb))/scale
        if (abund[i]=="Bs") {
          pred.index[i,] <- 0
          for (j in index.age[i]:max.age[i]) {
            pred.index[i,] <- pred.index[i,]+as.numeric(baa[j+1,]*saa[j+1,])/scale
          }
        }
        if (abund[i]=="Bf") {
          pred.index[i,] <- 0
          for (j in index.age[i]:max.age[i]) {
            pred.index[i,] <- pred.index[i,]+as.numeric(baa[j+1,]*obj$env$report()[["saa_f"]][j+1,,i+1])/scale
          }
        }
        if (is.na(b1[i])) {
          pred.index[i,] <- q1[i]*pred.index[i,]
        } else {
          pred.index[i,] <- q1[i]*pred.index[i,]^b1[i]
        }
      }

      output <- list(data=data, obj=obj, opt=opt, rep=rep, q=q1, b=b1, sigma=sigma2, sigma.logC=sigma1, sigma.logFsta=sigma3, sigma.logN=sigma4, rho=rho1, phi=phi, loglik=loglik,aic=aic,N=exp(logN),F=exp(logF),ref.year=rev(ref.year1),SR=SR.name,rec.par=rec.par, BRP0=BRP0,BRPmsy=BRPmsy,naa=naa,faa=faa,baa=baa,ssb=ssb,saa=saa,zaa=zaa,caa=caa,pred.index=pred.index,waa_est=waa_est)
      return(output)
    }

    brp1 <- BioRefPt(data, rep, opt, ref.year, is.SR, obj)
    brp1$par_list = obj$env$parList(opt$par)
  } else {
    brp1 <- list(data=data, obj=obj)
  }
  brp1$input <- arglist
  # brp1$tmbdata <- data
  brp1$init <- params
  brp1$map <- map

  class(brp1) <- "sam"

  return(brp1)
}

