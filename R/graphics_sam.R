#'
#' @import ggplot2
#' @import magrittr
#' @import dplyr
#' @import tidyr
#' @import tibble
#' @import frasyr
#' @import stringr
#' @import assertthat
#' @import purrr
#' @import TMB
#' @import frasyr
#' @importFrom magrittr %>%
#' @importFrom dplyr filter
#' @importFrom dplyr select
#' @importFrom stats na.omit
#' @importFrom stats qnorm
#' @importFrom stats resid
#' @importFrom stats rnorm
#' @importFrom stats nlminb
#' @importFrom frasyr theme_SH
#' @importFrom frasyr convert_vpa_tibble
#' @importFrom frasyr plot_vpa
#'
NULL

#' 条件付き負の対数尤度を成分別に描画する
#'
#' [get_cond_nll()] の出力を横向きの棒グラフで表示します。
#'
#' @param x [get_cond_nll()] が返すデータフレーム。
#' @param show_value 棒の外側に負の対数尤度を小数第1位まで表示するか。
#'   既定値は `TRUE`。
#'
#' @return `ggplot` オブジェクト。
#'
#' @examples
#' \dontrun{
#' data("samres_example", package = "frasam")
#' plot_cond_nll(get_cond_nll(samres))
#' }
#'
#' @export
plot_cond_nll <- function(x, show_value = TRUE) {
  if (!is.data.frame(x) || !all(c("type", "nll") %in% names(x))) {
    stop("'x' must be the output of get_cond_nll().", call. = FALSE)
  }
  if (!is.numeric(x$nll) || anyNA(x$nll) || any(!is.finite(x$nll))) {
    stop("'x$nll' must contain only finite numeric values.", call. = FALSE)
  }
  if (!is.logical(show_value) || length(show_value) != 1L || is.na(show_value)) {
    stop("'show_value' must be TRUE or FALSE.", call. = FALSE)
  }

  plot_data <- x
  plot_data$type <- factor(plot_data$type, levels = unique(as.character(plot_data$type)))
  plot_data$label <- sprintf("%.1f", plot_data$nll)
  plot_data$label_hjust <- ifelse(plot_data$nll >= 0, -0.1, 1.1)
  plot_data$type <- forcats::fct_rev(plot_data$type)

  p <- ggplot2::ggplot(
    plot_data,
    ggplot2::aes(x = nll, y = type, fill = type)
  ) +
    ggplot2::geom_col() +
    ggplot2::geom_vline(xintercept = 0, linewidth = 0.3) +
    ggplot2::scale_x_continuous(
      expand = ggplot2::expansion(mult = c(0.15, 0.15))
    ) +
    ggplot2::labs(x = "Conditional negative log-likelihood", y = NULL) +
    ggplot2::guides(fill = "none") +
    ggplot2::theme_bw()

  if (show_value) {
    p <- p + ggplot2::geom_text(
      ggplot2::aes(label = label, hjust = label_hjust),
      size = 3.5
    )
  }

  p
}

#' Update a future projection plot legend for SAM
#'
#' `frasyr::plot_futures()` labels the historical assessment line as VPA.
#' This helper updates only the legend labels of the returned ggplot object so
#' the historical line is shown as SAM when a SAM result was supplied.
#'
#' @param plot A ggplot object returned by \code{frasyr::plot_futures()}.
#' @param sam_label Label used for the historical SAM line.
#' @param scenario_labels Optional labels for future scenarios. Supply either a
#'   named character vector whose names match the current scenario names, or an
#'   unnamed vector with the same length as the non-SAM scenarios.
#' @param legend_title Legend title.
#' @param ncol_legend Number of columns in the colour legend.
#'
#' @return A ggplot object with updated colour and fill legend labels.
#' @export
#'
plot_update2sam <- function(plot,
                            sam_label = "SAM",
                            scenario_labels = NULL,
                            legend_title = "",
                            ncol_legend = 2) {
  if (!inherits(plot, "ggplot")) {
    stop("'plot' must be a ggplot object.", call. = FALSE)
  }
  if (is.null(plot$data) || !all(c("scenario", "col") %in% names(plot$data))) {
    stop("'plot' must contain 'scenario' and 'col' columns in plot$data.", call. = FALSE)
  }

  style_def <- plot$data %>%
    dplyr::ungroup() %>%
    dplyr::select(scenario, col, dplyr::any_of("lty")) %>%
    dplyr::filter(!is.na(col)) %>%
    dplyr::distinct(col, .keep_all = TRUE)

  labels <- as.character(style_def$scenario)
  is_sam_line <- is.na(labels) | labels == "VPA" | style_def$col == "black"
  labels[is_sam_line] <- sam_label

  if (!is.null(scenario_labels)) {
    scenario_labels <- as.character(scenario_labels)
    future_idx <- which(!is_sam_line)
    if (!is.null(names(scenario_labels)) && any(nzchar(names(scenario_labels)))) {
      matched <- match(labels[future_idx], names(scenario_labels))
      replace_idx <- future_idx[!is.na(matched)]
      labels[replace_idx] <- scenario_labels[matched[!is.na(matched)]]
    } else {
      if (length(scenario_labels) != length(future_idx)) {
        stop(
          "'scenario_labels' must have the same length as the non-SAM scenarios.",
          call. = FALSE
        )
      }
      labels[future_idx] <- scenario_labels
    }
  }

  lty <- if ("lty" %in% names(style_def)) style_def$lty else rep("solid", nrow(style_def))
  lty[is.na(lty)] <- "solid"

  plot +
    ggplot2::scale_color_identity(
      guide = "legend",
      breaks = style_def$col,
      labels = labels
    ) +
    ggplot2::scale_fill_identity(
      guide = "legend",
      breaks = style_def$col,
      labels = labels
    ) +
    ggplot2::guides(
      color = ggplot2::guide_legend(
        title = legend_title,
        ncol = ncol_legend,
        override.aes = list(linetype = lty, color = style_def$col, lwd = 0.7)
      ),
      fill = ggplot2::guide_legend(title = legend_title, ncol = 1)
    )
}

#' SAMで推定された再生産関係の予測値
#'
#' @param samres samの結果オブジェクト
#' @export
#'
get_predSR <- function(samres,max.ssb.pred=1.3,length=100){
  SR = samres$SR
  a = samres$rec.par["a"]
  b = samres$rec.par["b"]
  gamma = samres$input$gamma
  if(!is.null(samres$par_list[["rec_logk"]])) k = exp(samres$par_list[["rec_logk"]])

  if (SR=="HS") SRF <- function(x,a,b) ifelse(x>b,b*a,x*a)
  if (SR=="BH") SRF <- function(x,a,b) a*x/(1+b*x)
  if (SR=="RI") SRF <- function(x,a,b) a*x*exp(-b*x)
  if (SR=="Mesnil") SRF <- function(x,a,b) 0.5*a*(x+sqrt(b^2+gamma^2/4)-sqrt((x-b)^2+gamma^2/4))
  if (SR=="BHS") SRF <- function(x,a,b) ifelse(x<b,a*b*(x/b)^(1-(x/b)^k),a*b)

  assertthat::assert_that(SR %in% c("HS","BH","RI","Mesnil","BHS"))

  data_SR = frasyr::get.SRdata(samres) %>% as.data.frame()
  data_SR$SSB <- data_SR$SSB/samres$input$scale
  data_SR$R <- data_SR$R/samres$input$scale_number

  # rec.age>0のときにRをずらす
  # yearはSSBに合わせる（つまり何年生まれかを表す）
  data_SR <- shift_SRdata_rec_age(data_SR, samres$input$rec.age)

  ssb_c = seq(from=0,to=max(data_SR$SSB, na.rm = TRUE)*max.ssb.pred,length=length)
  R_c = purrr::map_dbl(ssb_c,SRF,a=a,b=b)
  pred_data = data.frame(SSB=ssb_c,R=R_c)

  list(pred=pred_data,est=data_SR)
  # g1 = ggplot(data=NULL,aes(x=SSB,y=R))+
  #   geom_path(data=pred_data)+
  #   geom_point(data=data_SR)+
  #   frasyr::theme_SH()
  # g1

}

#' Plot SAM or VPA results
#'
#' @param vpa_sam_list A SAM or VPA result object, or a list of result objects.
#' @param CI Confidence interval width. Set 0 to omit confidence intervals.
#' @param scenario_name Scenario names used in the legend.
#' @param alpha Alpha value for confidence interval ribbons.
#' @param size Line width.
#' @param base_size Base font size.
#' @param log_scale If \code{TRUE}, use a log scale on the y-axis.
#' @param legend_name Legend title.
#' @param legend_nrow Number of rows in the legend.
#' @param legend_position Legend position.
#' @param what.plot Statistics to plot.
#' @param years Years to include in the plot.
#' @param ncol Number of columns in the facet plot.
#' @param scale_recruitment Divisor for recruitment values in the plot.
#' @param scale_biomass Divisor for biomass values in the plot.
#' @param scale_ssb Divisor for spawning stock biomass values in the plot.
#' @param scale_catch Divisor for catch values in the plot.
#' @importFrom forcats fct_inorder
#' @export
#' @encoding UTF-8

plot_samvpa <- function(vpa_sam_list,CI=0.95,scenario_name=NULL,
                     alpha=0.4,size=1,base_size=14,log_scale=FALSE,
                     legend_name="Scenario",legend_nrow=1, legend_position="top",
                     what.plot = c("biomass","SSB","Recruitment","U"),
                     years = NULL,
                     ncol=2,
                     scale_recruitment = 1000,
                     scale_biomass = 1000,
                     scale_ssb = 1000,
                     scale_catch = 1000
){

  if(class(vpa_sam_list)[1] %in% c("sam","vpa")) {
    vpa_sam_list <- list(vpa_sam_list)
  }
  plot_scales <- c(
    scale_recruitment = scale_recruitment,
    scale_biomass = scale_biomass,
    scale_ssb = scale_ssb,
    scale_catch = scale_catch
  )
  if (any(!is.finite(plot_scales) | plot_scales <= 0)) {
    stop("All scale arguments must be positive finite values.")
  }

  g0 = frasyr::plot_vpa(vpa_sam_list)

  data = g0$data %>%
    mutate(id = forcats::fct_inorder(id))

  # data$stat %>% unique()
  # data$id %>% unique()
  data2 = data %>% dplyr::filter(stat %in% what.plot) %>%
    group_by(id,year,stat) %>%
    summarise(value=mean(value)) %>%
    ungroup() %>%
    mutate(value = case_when(
      stat == "Recruitment" ~ value / scale_recruitment,
      stat == "biomass" ~ value / scale_biomass,
      stat == "SSB" ~ value / scale_ssb,
      stat == "catch" ~ value / scale_catch,
      TRUE ~ value
    ))
  # data2 %>% filter(stat == "U")
  # class(data2$stat)
  data2 = data2 %>%
    mutate(stat = as.character(stat)) %>%
    mutate(stat2 = case_when(stat=="biomass" ~ "Biomass",
                             stat=="fishing_mortality"~"F",
                             stat=="U"~"Exploitation_rate",
                             stat=="catch"~"Catch",
                             TRUE ~ stat))
  what.plot_f = case_when(what.plot=="biomass" ~ "Biomass",
                          what.plot=="fishing_mortality"~"F",
                          what.plot=="U"~"Exploitation_rate",
                          what.plot=="catch"~"Catch",
                          TRUE ~ what.plot)
  data2 = data2 %>%
    mutate(stat_f = factor(stat2,levels=what.plot_f,labels=what.plot_f))

  if (!is.null(scenario_name)) {
    data2 = data2 %>% mutate(model = scenario_name[sapply(1:nrow(data2), function(i) which(data2$id[i]==unique(data2$id)))])
  } else{
    scenario_name = unique(as.character(data2$id))
    data2 = data2 %>% mutate(model = scenario_name[sapply(1:nrow(data2), function(i) which(data2$id[i]==unique(data2$id)))])
  }

  CVdata_all = tibble()
  stat_order = case_when(what.plot=="biomass" ~ "Biomass",
                          what.plot=="fishing_mortality"~"F",
                          what.plot=="U"~"Exploitation_rate",
                         what.plot=="catch"~"Catch",
                          TRUE ~ what.plot)
  for(i in 1:length(vpa_sam_list)) {
    res = vpa_sam_list[[i]]
    if(class(res)=="vpa") {
      # browser()
      if (is.null(res$rep)) {
        stop("Rerun vpa() with TMB=TRUE & sdreport=TRUE!")
      }
      cvdata = tibble(stat0 = names(res$rep$value),CV=res$rep$sd/res$rep$value,model=scenario_name[i]) %>%
        filter(stat0 != "U")
      if (!("U" %in% what.plot)) cvdata = cvdata %>% filter(stat0 != "scale_U")
      if (!("fishing_mortality" %in% what.plot)) cvdata = cvdata %>% filter(stat0 != "F_mean")
      cvdata = cvdata %>%
        mutate(Year = rep(as.numeric(colnames(res$naa)),nrow(cvdata)/ncol(res$naa)))

      cvdata_R = cvdata %>% filter(stat0=="N") %>%
        arrange(Year) %>%
        mutate(Age = rep(as.numeric(rownames(res$naa)),ncol(res$naa))) %>%
        filter(Age == min(Age)) %>%
        select(-Age) %>%
        mutate(stat = "Recruitment")

      # cvdata$stat0 %>% unique()
      cvdata = cvdata %>% filter(stat0 %in% c("SSB","B_total","F_mean","scale_U","catch")) %>%
        mutate(stat = case_when(stat0=="SSB" ~ "SSB",
                                stat0=="F_mean" ~ "F",
                                stat0=="scale_U" ~ "Exploitation_rate",
                                stat0 =="catch" ~ "Catch",
                                TRUE ~ "Biomass")) %>%
        full_join(cvdata_R) %>%
        mutate(stat_f = factor(stat,levels=stat_order)) %>%
        select(-stat0,stat)
      CVdata_all = bind_rows(CVdata_all,cvdata)
    }

    if (class(res)=="sam") {
      # exploitation_rateのCVを(U*(1-U))に変更（2024/09/10）
      if (is.null(res$rep$unbiased)){
        cvdata2 = tibble(stat0 = names(res$rep$value),CV=res$rep$sd/res$rep$value,model=scenario_name[i]) %>%
          filter(stat0 %in% c("exp_logN","ssb","B_total","F_mean","scale_U","Catch_biomass"))
      }else{
        cvdata2 = tibble(stat0 = names(res$rep$value),CV=res$rep$sd/res$rep$unbiased$value,model=scenario_name[i]) %>%
          filter(stat0 %in% c("exp_logN","ssb","B_total","F_mean","scale_U","Catch_biomass"))
      }
      cvdata2 = cvdata2 %>% mutate(stat0 = ifelse(stat0 == "scale_U","Exploitation_rate",stat0))

      cvdata2_R = cvdata2 %>% filter(stat0=="exp_logN") %>%
        mutate(Age = as.numeric(rep(rownames(res$naa),ncol(res$naa)))) %>%
        filter(Age == 0) %>%
        mutate(Year = as.numeric(colnames(res$naa))) %>%
        select(-Age) %>% mutate(stat = "Recruitment")

      cvdata2 = cvdata2 %>% filter(stat0 %in% c("ssb","B_total","F_mean","Exploitation_rate","Catch_biomass")) %>%
        # mutate(Year = rep(as.numeric(colnames(res$naa)),length(unique(.$stat0)))) %>%
        mutate(Year = rep(as.numeric(colnames(res$naa)),length(unique(.$stat0)))) %>%
        mutate(stat = case_when(stat0=="ssb" ~ "SSB",
                                stat0=="F_mean" ~ "F",
                                stat0=="B_total" ~ "Biomass",
                                stat0=="Catch_biomass" ~ "Catch",
                                TRUE~stat0)) %>%
        full_join(cvdata2_R)
      cvdata2 = cvdata2 %>% filter(stat %in% stat_order) %>%
        mutate(stat_f = factor(stat,levels=stat_order)) %>%
        select(-stat0,stat)
      # cvdata2$stat %>% unique
      CVdata_all = bind_rows(CVdata_all,cvdata2)
    }
  }

  data2 = data2 %>% rename(Year = year) %>%
    dplyr::select(-stat)

  CVdata_all = CVdata_all %>% dplyr::select(-stat)
  # CVdata_all

  # Exploitation rateについてはU/1-UのCVから計算に変更（2024/09/10）
  data3 = full_join(data2,CVdata_all) %>%
    mutate(stat_f = factor(stat_f,levels=stat_order)) %>%
    mutate(Cz = ifelse(stat_f != "Exploitation_rate",exp(qnorm(CI+(1-CI)/2)*sqrt(log(1+CV^2))),exp(qnorm(CI+(1-CI)/2)*CV))) %>%
    mutate(lower = ifelse(stat_f != "Exploitation_rate", value/Cz, value/(value+(1-value)*Cz)),
           upper = ifelse(stat_f != "Exploitation_rate", value*Cz, value/(value+(1-value)/Cz))) %>%
    arrange(model,stat_f,Year) %>%
    mutate(Model = factor(model,levels=scenario_name))

  if (!is.null(years)) data3 = data3 %>% dplyr::filter(Year %in% years)

  if (CI==0) {
    g1 = ggplot(data=data3,aes(x=Year,y=value))
  }else{
    g1 = ggplot(data=data3,aes(x=Year,y=value))+
      geom_ribbon(aes(ymax=upper,ymin=lower,fill=Model),alpha=alpha)+
      scale_fill_brewer(palette="Set1",name=legend_name)
  }

  g1 = g1 +
    geom_path(aes(colour=Model,linetype=Model),linewidth=size)+
    facet_wrap(vars(stat_f),scales="free_y",ncol=ncol)+
    theme_SH()+theme_bw(base_size=base_size)+theme(legend.position=legend_position)+
    xlab("Year") + ylab("")+
    # ylim(0,NA)
    scale_colour_brewer(palette="Set1",name=legend_name)+
    scale_linetype_discrete(name=legend_name)+
    guides(colour=guide_legend(nrow=legend_nrow),
           fill=guide_legend(nrow=legend_nrow),
           linetype=guide_legend(nrow=legend_nrow))+
    scale_x_continuous(breaks=scales::pretty_breaks())
  if (isTRUE(log_scale)) {
    g1 = g1 + scale_y_log10()
  } else {
    g1 = g1 + ylim(0,NA)
  }

  g1
}

#' レトロスペクティブ解析の結果をプロットする
#' @param res sam object
#' @param retro_res \code{retro_sam}の結果
#' @param start_year プロットを開始する年
#' @import dplyr
#' @export

retro_plot = function(res,retro_res,start_year=NULL,scale=1000, forecast=FALSE,
                      base_size=14, plot_mohn=TRUE, mohn_position = "upperleft") {

  mohn_res = retro_res$mohn
  if (isTRUE(forecast)) {
    mohn_res = retro_res$mohn_forecast
    if (class(res)=="vpa") warning("'forecast=TRUE' is not possible for VPA")
  }

  res_tibble = convert_sam_tibble(res) %>% mutate(id = 0)
  maxyear = res_tibble$year %>% max
  # if(isTRUE(forecast)) maxyear <- maxyear + 1
  # if (isTRUE(res$input$last.catch.zero)) {
  #   res_tibble = res_tibble %>% dplyr::filter(year<maxyear)
  #   maxyear <- maxyear - 1
  # }

  for (i in 1:length(retro_res$Res)) {
    if (isTRUE(res$input$last.catch.zero)) {
      res_tibble <- full_join(res_tibble, convert_sam_tibble(retro_res$Res[[i]]) %>%
                                mutate(id = i) %>% dplyr::filter(year<as.numeric(forecast)+maxyear-id))
    } else {
      res_tibble <- full_join(res_tibble, convert_sam_tibble(retro_res$Res[[i]]) %>% mutate(id = i))
    }
  }

  res_tibble -> data
  # browser()
  if (is.null(start_year)) start_year <- max(as.numeric(colnames(res$faa)))-15
  data2 = data %>% dplyr::filter(stat %in% c("SSB","biomass","Recruitment","fishing_mortality")) %>%
    group_by(id,year,stat) %>%
    summarise(value=mean(value)) %>%
    ungroup() %>%
    mutate(value = if_else(stat == "fishing_mortality",value,value/scale)) %>%
    mutate(terminal_year = max(data$year)-id+as.numeric(forecast)-as.numeric(res$input$last.catch.zero))
  # if (!(isTRUE(forecast) & class(res)[1]=="sam")) {
  #   data2 = data2 %>% mutate(terminal_year = terminal_year-1)  }
  data2 = data2 %>% dplyr::filter(year <= terminal_year) %>%
    mutate(stat_f = factor(stat,levels=c("biomass","SSB","Recruitment","fishing_mortality"),labels=c("Biomass","SSB","Recruitment","F"))) %>%
    filter(year >= start_year) %>%
    filter(year <= as.numeric(terminal_year)) %>%
    mutate(terminal_year = factor(terminal_year, levels=unique(terminal_year)))
  data_full = data2 %>% filter(id==0)
  data_others = data2 %>% filter(id>0)
  data_term = data2 %>%
    mutate(term_year = as.numeric(as.character(terminal_year))) %>%
    filter(year==term_year)
  g1 = ggplot(data=data_others,aes(x=year,y=value))+
    geom_path(data=data_full,colour="black",linewidth=1)+
    geom_path(aes(colour=terminal_year),linewidth=1)+
    facet_wrap(vars(stat_f),scales="free_y",ncol=2)+
    geom_point(data=filter(data_term,id==0),colour="black",size=2)+
    geom_point(data=filter(data_term,id>0),aes(colour=terminal_year),size=2)+
    theme_SH()+theme_bw(base_size=base_size)+theme(legend.position="none")+
    xlab("Year") + ylab("")+
    ylim(0,NA) +
    scale_x_continuous(breaks=scales::breaks_pretty())
  if (isTRUE(plot_mohn)) {
    mohn = tibble(rho = mohn_res, stat = names(mohn_res)) %>%
      dplyr::filter(stat != "N") %>%
      mutate(stat_f = factor(stat,levels=c("B","SSB","R","F"),labels=c("Biomass","SSB","Recruitment","F"))) %>%
      mutate(label=sprintf("rho == %.2f",rho)) %>%
      full_join(data2 %>% group_by(stat_f) %>% summarise(ymax = max(value)))
    if (mohn_position=="upperleft") {
      g1 = g1 + geom_text(data=mohn,parse=TRUE,aes(x=start_year,y=ymax,label=label,hjust=0,vjust=1))
    } else {
      if (mohn_position=="upperright") {
        end_year = data2$year %>% max
        g1 = g1 + geom_text(data=mohn,parse=TRUE,aes(x=end_year,y=ymax,label=label,hjust=1,vjust=1))
      } else {
        if (mohn_position=="bottomleft") {
          g1 = g1 + geom_text(data=mohn,parse=TRUE,aes(x=start_year,y=0,label=label,hjust=0,vjust=0))
        } else {
          if (mohn_position=="bottomright") {
            end_year = data2$year %>% max
            g1 = g1 + geom_text(data=mohn,parse=TRUE,aes(x=end_year,y=0,label=label,hjust=1,vjust=0))
          } else {
            stop("Inappropriate 'mohn_position'!")
          }
        }
      }
    }
  }
  g1
}

#' Indexの当てはまりについてプロットする関数
#'
#' @export

index_plot = function(samvpa_list,model_name=NULL, fleet_no = NULL,
                      scales=c("free","free_x","free"),base_size=16) {
  if (class(samvpa_list)[1] == "list") {
    nmodel = length(samvpa_list)
    dat_index = samvpa_list[[1]]$input$dat$index
    for(i in 1:(nmodel-1)) {
      dat_index2 = samvpa_list[[i+1]]$input$dat$index
      if (sum(dat_index2-dat_index,na.rm=T)!=0) stop("The models using different datasets can not be compared!")
    }
  }else{
    nmodel = 1
    dat = samvpa_list$input$dat
    samvpa_list[[1]] <- samvpa_list
    dat_index = samvpa_list[[1]]$input$dat$index
  }
  if(is.null(model_name)) model_name = as.character(1:nmodel)
  index_obs = as_tibble(dat_index) %>%
    mutate(id = 1:n()) %>%
    pivot_longer(cols=-id,names_to="Year",values_to="obs")
  for (i in 1:nmodel) {
    res = samvpa_list[[i]]
    index_obs = index_obs %>% full_join(
      as_tibble(res$pred.index) %>% mutate(id=1:n()) %>%
        pivot_longer(cols=-id,names_to="Year",values_to=model_name[i])
    )
  }

  if (is.null(fleet_no)) fleet_no = 1:nrow(dat_index)

  index_pred = index_obs %>%
    pivot_longer(cols=all_of(model_name),names_to = "Model",values_to="pred") %>%
    na.omit() %>%
    mutate(resid = log(obs/pred)) %>%
    mutate(Year = as.numeric(Year)) %>%
    mutate(Fleet = str_c("Fleet ",as.character(fleet_no[id])))

  index_obs = na.omit(index_obs) %>%
    mutate(Year = as.numeric(Year)) %>%
    mutate(Fleet = str_c("Fleet ",as.character(fleet_no[id])))

  g_index = ggplot(data=NULL,aes(x=Year))+
    geom_point(data=index_obs,aes(y=obs),colour="black",size=1.5)+
    geom_path(data=index_pred,aes(y=pred,colour=Model),linewidth=1)+
    facet_wrap(vars(Fleet),nrow=2,scales=scales[1])+
    scale_colour_brewer(palette="Set1",name="")+
    theme_bw(base_size=base_size)+theme(legend.position="top")+
    ylab("Index value")+ylim(0,NA)
  # g_index

  g_resid = ggplot(data=index_pred,aes(x=Year,y=resid)) +
    # geom_point(data=index_obs,aes(y=obs),size=1.5) +
    geom_point(aes(colour=Model,shape=Model),size=1.5)+
    facet_wrap(vars(Fleet),nrow=2,scales=scales[2])+
    theme_bw(base_size=base_size)+theme(legend.position="top")+
    scale_colour_brewer(palette="Set1",name="")+
    scale_shape_discrete(name="")+
    geom_hline(yintercept=0)+
    stat_smooth(level=0.8,se=FALSE,aes(group=Model,colour=Model))+
    ylab("Residual")

  index_pred2 = index_pred %>%
    mutate(model_no = map_dbl(1:n(), function(i) which(model_name==Model[i])))
  q_tmp = sapply(1:nrow(index_pred2), function(i) samvpa_list[[index_pred2$model_no[i]]]$q[index_pred2$id[i]])
  b_tmp = sapply(1:nrow(index_pred2), function(i) samvpa_list[[index_pred2$model_no[i]]]$b[index_pred2$id[i]])
  index_pred2 = index_pred2 %>% mutate(q=q_tmp,b=b_tmp) %>%
    mutate(abund = (pred/q)^(1/b))

  # colnames(index_pred2)

  index_pred_tmp = index_pred2 %>% group_by(Fleet,Model,model_no,id,q,b) %>%
    summarise(max_abund = max(abund))

  index_pred_curve = map_dfr(1:nrow(index_pred_tmp), function(i) {
   data.frame(Fleet=index_pred_tmp$Fleet[i],Model=index_pred_tmp$Model[i],
              model_no=index_pred_tmp$model_no[i],id=index_pred_tmp$id[i],
              q=index_pred_tmp$q[i],b=index_pred_tmp$b[i],abund=seq(0,index_pred_tmp$max_abund[i],length=201)) %>%
      mutate(pred = q*abund^b)
  })

  g_abund = ggplot(data=NULL,aes(x=abund))+
    geom_path(data=index_pred_curve,aes(y=pred,colour=Model),linewidth=1)+
    geom_point(data=index_pred2,aes(y=obs,colour=Model),size=1.5)+
    facet_wrap(vars(Fleet),nrow=2,scales=scales[3])+
    theme_bw(base_size=base_size)+theme(legend.position="top")+
    scale_colour_brewer(palette="Set1",name="")+
    ylab("Index")+xlab("Abundance")

  # g_abund

  return( list(index=g_index,resid=g_resid,abund=g_abund) )
}

#' Indexの当てはまりについてプロットする関数
#'
#' @export

index_plot2 = function(samres, index_name = NULL,nrow=2,
                      scales=c("free","free_x","free"),base_size=14) {
  # browser()
  dat_index = samres$input$dat$index
  index_obs = as_tibble(dat_index) %>%
    mutate(id = 1:n()) %>%
    pivot_longer(cols=-id,names_to="Year",values_to="obs")
  res = samres
  index_obs = index_obs %>% full_join(
      as_tibble(res$pred.index) %>% mutate(id=1:n()) %>%
        pivot_longer(cols=-id,names_to="Year",values_to="pred")
    )

  if (is.null(index_name)) index_name = as.character(1:nrow(dat_index))

  index_pred = index_obs %>%
    # pivot_longer(cols=all_of(model_name),names_to = "Model",values_to="pred") %>%
    na.omit() %>%
    mutate(resid = log(obs/pred)) %>%
    mutate(Year = as.numeric(Year)) %>%
    mutate(Fleet = index_name[id]) %>%
    mutate(Fleet = factor(Fleet,levels=index_name))

  index_obs = na.omit(index_obs) %>%
    mutate(Year = as.numeric(Year)) %>%
    mutate(Fleet = index_name[id]) %>%
    mutate(Fleet = factor(Fleet,levels=index_name))

  g_index = ggplot(data=NULL,aes(x=Year))+
    geom_point(data=index_obs,aes(y=obs),colour="black",size=1.5)+
    geom_path(data=index_pred,aes(y=pred),colour="blue",linewidth=1)+
    facet_wrap(vars(Fleet),nrow=nrow,scales=scales[1])+
    scale_colour_brewer(palette="Set1",name="")+
    theme_bw(base_size=base_size)+theme(legend.position="top")+
    ylab("Index value")+ylim(0,NA)
  # g_index

  g_resid = ggplot(data=index_pred,aes(x=Year,y=resid)) +
    # geom_point(data=index_obs,aes(y=obs),size=1.5) +
    geom_point(size=1.5)+
    facet_wrap(vars(Fleet),nrow=nrow,scales=scales[2])+
    theme_bw(base_size=base_size)+theme(legend.position="top")+
    scale_colour_brewer(palette="Set1",name="")+
    scale_shape_discrete(name="")+
    geom_hline(yintercept=0)+
    stat_smooth(level=0.8,se=TRUE)+
    ylab("Residual")

  index_pred2 = index_pred
  # %>% mutate(model_no = map_dbl(1:n(), function(i) which(model_name==Model[i])))
  q_tmp = sapply(1:nrow(index_pred2), function(i) samres$q[index_pred2$id[i]])
  b_tmp = sapply(1:nrow(index_pred2), function(i) samres$b[index_pred2$id[i]])
  index_pred2 = index_pred2 %>% mutate(q=q_tmp,b=b_tmp) %>%
    mutate(abund = (pred/q)^(1/b))

  # colnames(index_pred2)

  index_pred_tmp = index_pred2 %>% group_by(Fleet,id,q,b) %>%
    summarise(max_abund = max(abund))

  index_pred_curve = map_dfr(1:nrow(index_pred_tmp), function(i) {
    data.frame(Fleet=index_pred_tmp$Fleet[i],id=index_pred_tmp$id[i],
               q=index_pred_tmp$q[i],b=index_pred_tmp$b[i],abund=seq(0,index_pred_tmp$max_abund[i],length=201)) %>%
      mutate(pred = q*abund^b)
  })

  g_abund = ggplot(data=NULL,aes(x=abund))+
    geom_path(data=index_pred_curve,aes(y=pred),colour="blue",linewidth=1)+
    geom_point(data=index_pred2,aes(y=obs),size=1.5)+
    facet_wrap(vars(Fleet),nrow=nrow,scales=scales[3])+
    theme_bw(base_size=base_size)+theme(legend.position="top")+
    scale_colour_brewer(palette="Set1",name="")+
    ylab("Index")+xlab("Abundance")

  # g_abund

  return( list(index=g_index,resid=g_resid,abund=g_abund) )
}


#' Catch at ageの当てはまりについてプロットする関数
#'
#' @export

caa_plot = function(samres,
         scales=c("free","free_x"),base_size=16) {
  dat = samres$input$dat
  caa_obs = as_tibble(dat$caa) %>% mutate(Age=0:(n()-1)+samres$input$rec.age) %>%
    pivot_longer(cols=-Age,names_to="Year",values_to="obs") %>%
    mutate(Year = as.numeric(Year))

  if(isTRUE(samres$input$last.catch.zero)) caa_obs = caa_obs %>% filter(Year < max(Year))

  caa_pred = as_tibble(samres$caa) %>% mutate(Age=0:(n()-1)+samres$input$rec.age) %>%
    pivot_longer(cols=-Age,names_to="Year",values_to="pred") %>%
    mutate(Year = as.numeric(Year))

  maxage = caa_obs$Age %>% max

  caa_dat = full_join(caa_obs,caa_pred) %>%
    mutate(resid = log(obs/pred)) %>%
    mutate(Age = if_else(Age==maxage,str_c("Age ",Age,"+"),str_c("Age ",Age)))

  g_caa = ggplot(data=caa_dat,aes(x=Year))+
    geom_point(aes(y=obs),size=1.5) +
    geom_path(aes(y=pred),linewidth=0.8)+
    facet_wrap(vars(Age),nrow=2,scales=scales[1]) +
    theme_bw(base_size=base_size)+ylab("Catch at age")

  # g_caa
  #
  g_caa_resid = ggplot(data=caa_dat,aes(x=Year,y=resid))+
    geom_point(size=1.5) +
    facet_wrap(vars(Age),nrow=2,scales=scales[2]) +
    geom_hline(yintercept=0)+
    stat_smooth(aes(group=Age),size=0.8)+ylab("Residual")+
    theme_bw(base_size=base_size)

  # g_caa_resid

  return( list(caa = g_caa,resid=g_caa_resid) )

}

#' 再生産関係についてプロットする関数
#'
#' @param samres SAMの結果オブジェクト
#' @export

plot_SR_simple = function(samres,length=1000,base_size=14,...) {
  get_SR = get_predSR(samres,length=length,...)
  g_SR = ggplot(data=NULL,aes(x=SSB,y=R))+
    geom_point(data=get_SR$est,colour="darkgray",size=2,alpha=0.5)+
    geom_path(data=get_SR$pred,linewidth=1)+
    theme_bw(base_size=base_size)
  return( g_SR )
}


#' ブートストラップについてプロットする関数

#' @param samres \code{sam()}の結果オブジェクト
#' @param boores \code{boo_sam}の結果オブジェクト
#' @param draw_deltaCI デルタ法の信頼区間を描くか（デフォルト:FALSE）
#' @param scenario_name 推定値の結果、ブートストラップの結果の凡例に使う名前
#' @inheritParams plot_samvpa
#' @export

plot_boosam = function(samres,
                       boores,
                       CI = 0.95,
                       what.plot = c("biomass", "SSB", "Recruitment", "U"),
                       draw_deltaCI = FALSE,
                       scenario_name = c("Estiamte", "Bootstrap"),
                       alpha=0.4,
                       size=1,
                       base_size=14,
                       log_scale=FALSE,
                       legend_name="Scenario",
                       legend_nrow=1,
                       legend_position="top",
                       years = NULL,
                       ncol=2
                       ) {

  est_res = plot_samvpa(samres,CI=CI, what.plot = what.plot)

  est_data = est_res$data %>%
    mutate(scenario = scenario_name[1])

  if(!isTRUE(draw_deltaCI)) {
    est_data$lower <- est_data$upper <- NA
  }


  boo_res = plot_samvpa(boores, what.plot = what.plot, CI=0)
  boo_data = boo_res$data

  boo_data2 = boo_data %>% group_by(Year,stat2,stat_f) %>%
    summarise(median = median(value),
              mean = mean(value),
              SD = sd(value),
              CV = SD/mean) %>% suppressMessages() %>%
    ungroup() %>%
    rename(value = median)

    # Exploitation rateについてはU/1-UのCVから計算に変更（2024/09/10）
    boo_data3 = boo_data2 %>%
      mutate(Cz = ifelse(stat_f != "Exploitation_rate",exp(qnorm(CI+(1-CI)/2)*sqrt(log(1+CV^2))),exp(qnorm(CI+(1-CI)/2)*CV))) %>%
      mutate(lower = ifelse(stat_f != "Exploitation_rate", value/Cz, value/(value+(1-value)*Cz)),
             upper = ifelse(stat_f != "Exploitation_rate", value*Cz, value/(value+(1-value)/Cz))) %>%
      mutate(scenario = scenario_name[2])

  data3 = bind_rows(boo_data3,est_data) %>%
    mutate(Model = fct_inorder(scenario))

  if (CI==0) {
    g1 = ggplot(data=data3,aes(x=Year,y=value))
  }else{
    g1 = ggplot(data=data3,aes(x=Year,y=value))+
      geom_ribbon(aes(ymax=upper,ymin=lower,fill=Model),alpha=alpha)+
      scale_fill_brewer(palette="Set1",name=legend_name)
  }

  g1 = g1 +
    geom_path(aes(colour=Model,linetype=Model),linewidth=size)+
    facet_wrap(vars(stat_f),scales="free_y",ncol=ncol)+
    frasyr::theme_SH()+theme_bw(base_size=base_size)+theme(legend.position=legend_position)+
    xlab("Year") + ylab("")+
    # ylim(0,NA)
    scale_colour_brewer(palette="Set1",name=legend_name)+
    scale_linetype_discrete(name=legend_name)+
    guides(colour=guide_legend(title=NULL, nrow=legend_nrow),
           fill=guide_legend(title=NULL, nrow=legend_nrow),
           linetype=guide_legend(title=NULL, nrow=legend_nrow))+
    scale_x_continuous(breaks=scales::pretty_breaks())
  if (isTRUE(log_scale)) {
    g1 = g1 + scale_y_log10()
  } else {
    g1 = g1 + ylim(0,NA)
  }

  g1

}



#' OSA residualをプロットする関数
#'
#' @param osares \code{do_osa_resid}の結果オブジェクト
#'
#' @encoding UTF-8
#'
#' @export

plot_osa_resid <- function(osares) {

  ## caa
  osa_resid_caa = osares %>%
    filter(fleet==1) %>% #Fleet=1がCatch at age, それ以外がIndexを表す
    filter(!is.nan(residual))

  g_caa = osa_resid_caa %>%
    ggplot(aes(x=year,y=age,colour=residual,size=abs(residual))) +
    geom_point() +
    scale_colour_gradient2(high="red",low="blue",mid="gray")+
    xlab("Year")+
    ylab("Age")+
    scale_y_continuous(breaks=0:100)+
    theme_bw()

  ## index
  osa_resid_index = osares %>%
    filter(fleet>1) %>%
    filter(!is.nan(residual))

  g_index = osa_resid_index %>%
    ggplot(aes(x=year,y=fleet-1,colour=residual,size=abs(residual))) +
    geom_point() +
    scale_colour_gradient2(high="red",low="blue",mid="gray")+
    xlab("Fishing year")+
    ylab("Index ID")+
    scale_y_continuous(breaks=1:100)+
    theme_bw()

  ## qq plot
  p <- osares %>%
    filter(!is.nan(residual)) %>%
    ggplot(aes(sample = residual)) +
    stat_qq(distribution = qnorm) +  # 標準正規分布の分位数を使用
    stat_qq_line(distribution = qnorm) +  # 標準正規分布に基づく直線
    labs(x = "Theoretical Quantiles",
         y = "Sample Quantiles") +
    theme_SH()

  return(list(caa = g_caa, index = g_index, qq = p))
}

#' Popsimの結果ををプロットする関数
#'
#' @inheritParams sumup_popsim
#'
#' @encoding UTF-8
#'
#' @export

plot_popsim = function(
    res_true,
    fit2PS,
    CI = 0.95,
    what.plot = c("biomass","SSB","Recruitment","U"),
    Age =NULL,
    scenario_name = c("True", "Simulation"),
    alpha=0.4,
    size=1,
    base_size=14,
    log_scale=FALSE,
    legend_name="Scenario",
    legend_nrow=1,
    legend_position="top",
    years = NULL,
    ncol=2,
    sim_colour = "blue") {

  tmp = sumup_popsim(res_true,fit2PS,CI=CI)
  tbl_wide = tmp$summary %>%
    filter(stat %in% what.plot) %>%
    filter(!is.null(age) | age %in% Age) %>%
    mutate(stat_f = factor(stat,levels=what.plot)) %>%
    arrange(stat_f) %>%
    mutate(stat2 = ifelse(is.na(age), stat,
                          mutate(str_c(stat,as.character(age))))) %>%
    mutate(stat_f2 = fct_inorder(stat2)) %>%
    arrange(stat_f2,year) %>%
    mutate(year = as.integer(year))

   ggplot(tbl_wide,aes(x=year,group=stat_f2)) +
    geom_ribbon(aes(ymin=lower,ymax=upper),alpha=alpha,fill=sim_colour) +
    facet_wrap(vars(stat_f2),scales="free_y",ncol=ncol) +
    geom_path(aes(y=Median),colour=sim_colour,linewidth=0.7) +
    lemon::geom_pointline(aes(y=value_true), distance = 0,
                          linetype="dotted", linewidth=0.5)+
    ylim(0,NA) + ylab("Value") + xlab("Year")+
    frasyr::theme_SH()+
    scale_x_continuous(breaks=scales::pretty_breaks())
}


#' レトロの結果を使って資源量指標値に対するhindcast cross validationをプロットする関数
#'
#' @param show_mase MASEの結果を載せるかどうか
#' @param use_index 特定のIndexを使う場合、\code{use_index = c(1,3)}のように指定する
#' @param years プロットする期間を指定する場合、\code{years = 2015:2024}のように指定する
#' @inheritParams calc_mase
#'
#'
#' @encoding UTF-8
#'
#' @export

plot_hindcastCV = function(samres,
                           retrores,
                           h=1,
                           log = FALSE,
                           index_name = NULL,
                           show_mase = TRUE,
                           mase_position = "upperright",
                           use_index = NULL,
                           years = NULL
) {
  res_mase = calc_mase(samres = samres,
                       retrores = retrores,
                       h = h,
                       log = log,
                       index_name = index_name)

  if(is.null(use_index)) {
    use_index = 1:nrow(samres$input$dat$index)
  }

  dat_removed = res_mase$removed %>%
    filter(idx %in% use_index)
  dat_full = res_mase$full %>% filter(idx %in% use_index)

  if(!is.null(years)) {
    dat_removed = filter(dat_removed, year %in% years)
    dat_full = filter(dat_full, year %in% years)
  }
  dat_mase = res_mase$mase %>% filter(idx %in% use_index)
  dat_cv = res_mase$cv %>% filter(idx %in% use_index)

  dat_both = bind_rows(dat_full,dat_removed) %>%
    group_by(idx, index) %>%
    summarise(ymax = max(obs, pred_full, pred_cond, na.rm=T)) %>%
    ungroup()

  if (!mase_position %in% c("upperright", "upperleft", "bottomright", "bottomleft")) {
    stop("Invalid value for 'mase_position'. Must be one of 'upperright', 'upperleft', 'bottomright', or 'bottomleft'.")
  }

  dat_mase2 = left_join(dat_mase,dat_both) %>%
    mutate(year = max(dat_full$year), #upperright
           y = ymax,
           label = sprintf("MASE == %.2f",MASE)) %>%
    mutate(hjust=1,vjust=1)

  if(mase_position == "upperleft") {
    dat_mase2 = dat_mase2 %>%
      mutate(year = min(dat_full$year),
             hjust = 0)
  }
  if(str_detect(mase_position,"bottom")) {
    dat_mase2 = dat_mase2 %>%
      mutate(y = 0, vjust = 0)
    if(str_detect(mase_position,"left")) {
      dat_mase2 = dat_mase2 %>%
        mutate(vjust = 0, hjust = 0,
               year = min(dat_full$year))
    }
  }

  gg <- ggplot(data = dat_removed, aes(x=year)) +
    geom_path(linewidth=0.8,aes(y=pred_cond,colour=as.factor(retro_id)))+
    geom_point(data = dat_cv, size=2,
               aes(x = year_target, y=pred_cond,colour=as.factor(retro_id)))+
    geom_path(data=dat_full,aes(y=pred_full),colour="black",linewidth=0.8)+
    geom_point(data=dat_full,aes(y=obs),colour="black",size=2)+
    frasyr::theme_SH()+
    scale_x_continuous(breaks=scales::pretty_breaks()) +
    ylab("Value") + xlab("Year") +
    ylim(0,NA)

  if (length(use_index) > 1) {
    gg <- gg + facet_wrap(vars(index), scales="free_y")
  }

  if (isTRUE(show_mase)) {
    gg <- gg +
      geom_text(data=dat_mase2,parse=TRUE,size=4,
                aes(y=y,label=label,hjust=hjust,vjust=vjust))

  }
  gg
}



#' Plot age-aggregated biomass factors
#'
#' Draws a stacked bar chart of the age-aggregated contributions returned by
#' [decompose_biomass_factors()]. Positive and negative contributions are
#' stacked on opposite sides of zero.
#'
#' @param x Output from [decompose_biomass_factors()].
#' @param type Output scale. `"percent"` (default) uses
#'   `percent_aggregated`; `"absolute"` uses `age_aggregated`.
#' @param scale Positive divisor applied to values when `type = "absolute"`.
#'
#' @return A `ggplot` object. The `effect` variable in the plot data is a
#'   factor whose levels follow the row order of the aggregated matrix.
#'
#' @export
plot_biomass_factors <- function(x,
                                 type = c("percent", "absolute"),
                                 scale = 1) {
  type <- match.arg(type)
  if (length(scale) != 1L || !is.finite(scale) || scale <= 0) {
    stop("scale must be one positive finite number.", call. = FALSE)
  }

  component <- if (type == "percent") "percent_aggregated" else "age_aggregated"
  values <- x[[component]]

  if (is.null(values) || !is.matrix(values)) {
    stop("x must be an output from decompose_biomass_factors().", call. = FALSE)
  }
  if (is.null(rownames(values)) || is.null(colnames(values))) {
    stop(component, " must have row and column names.", call. = FALSE)
  }

  effect_levels <- rownames(values)
  years <- suppressWarnings(as.numeric(colnames(values)))
  if (any(!is.finite(years))) {
    stop("Column names of aggregated results must be numeric years.", call. = FALSE)
  }

  plot_data <- data.frame(
    effect = factor(rep(effect_levels, times = ncol(values)),
                    levels = effect_levels),
    year = rep(years, each = nrow(values)),
    value = as.vector(values),
    stringsAsFactors = FALSE
  )

  if (!is.null(x$annual_change) && length(x$annual_change) == ncol(values)) {
    valid_years <- years[!is.na(x$annual_change)]
    plot_data <- plot_data[plot_data$year %in% valid_years, , drop = FALSE]
  }
  plot_data <- plot_data[is.finite(plot_data$value), , drop = FALSE]
  if (type == "absolute") plot_data$value <- plot_data$value / scale

  y_label <- if (type == "percent") {
    quantity <- if (identical(x$target, "ssb")) {
      "SSB"
    } else {
      "biomass"
    }
    paste0("Contribution to ", quantity, " change (%)")
  } else {
    quantity <- if (identical(x$target, "ssb")) {
      "SSB"
    } else {
      "biomass"
    }
    paste0("Contribution to ", quantity, " change")
  }

  year_breaks <- seq(
    ceiling(min(plot_data$year) / 5) * 5,
    floor(max(plot_data$year) / 5) * 5,
    by = 5
  )

  ggplot2::ggplot(
    plot_data,
    ggplot2::aes(x = year, y = value, fill = effect)
  ) +
    ggplot2::geom_col() +
    ggplot2::geom_hline(yintercept = 0, linewidth = 0.3) +
    ggplot2::scale_x_continuous(breaks = year_breaks) +
    ggplot2::labs(x = "Year", y = y_label, fill = "Factor") +
    ggplot2::theme_bw()
}
#' Convert matrix or data.frame to numeric matrix
#'
#' @keywords internal
.as_numeric_matrix <- function(x, name = deparse(substitute(x))) {
  if (is.data.frame(x)) {
    x <- as.data.frame(x, check.names = FALSE)

    x <- lapply(x, function(z) {
      if (is.factor(z)) z <- as.character(z)
      suppressWarnings(as.numeric(z))
    })

    x <- as.data.frame(x, check.names = FALSE)
    x <- as.matrix(x)
  }

  if (!is.matrix(x)) {
    stop(name, " must be a matrix or data.frame.", call. = FALSE)
  }

  storage.mode(x) <- "numeric"

  if (!is.numeric(x)) {
    stop(name, " could not be converted to a numeric matrix.", call. = FALSE)
  }

  if (anyNA(x)) {
    warning(name, " contains NA after numeric conversion.", call. = FALSE)
  }

  x
}
