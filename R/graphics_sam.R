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

#' SAMで推定された再生産関係の予測値
#'
#' @param samres samの結果オブジェクト
#' @export
#'
get_predSR <- function(samres,max.ssb.pred=1.3){
  SR = samres$SR
  a = samres$rec.par["a"]
  b = samres$rec.par["b"]
  gamma = samres$input$gamma

  if (SR=="HS") SRF <- function(x,a,b) ifelse(x>b,b*a,x*a)
  if (SR=="BH") SRF <- function(x,a,b) a*x/(1+b*x)
  if (SR=="RI") SRF <- function(x,a,b) a*x*exp(-b*x)
  if (SR=="Mesnil") SRF <- function(x,a,b) 0.5*a*(x+sqrt(b^2+gamma^2/4)-sqrt((x-b)^2+gamma^2/4))

  assertthat::assert_that(SR %in% c("HS","BH","RI","Mesnil"))

  data_SR = get.SRdata(samres) %>% as.data.frame()

  ssb_c = seq(from=0,to=max(data_SR$SSB)*max.ssb.pred,length=100)
  R_c = purrr::map_dbl(ssb_c,SRF,a=a,b=b)
  pred_data = data.frame(SSB=ssb_c,R=R_c)

  list(pred=pred_data,est=data_SR)
  # g1 = ggplot(data=NULL,aes(x=SSB,y=R))+
  #   geom_path(data=pred_data)+
  #   geom_point(data=data_SR)+
  #   frasyr::theme_SH()
  # g1

}

#' SAM or VPAの結果を描くグラフ
#'
#' @param samres samの結果オブジェクト
#' @export
#' @encoding UTF-8

plot_samvpa <- function(vpa_sam_list,CI=0.95,scenario_name=NULL,
                     alpha=0.4,size=1,base_size=16,log_scale=FALSE,
                     legend_name="Scenario",legend_nrow=1, legend_position="top",
                     what.plot = c("biomass","SSB","Recruitment","U"),years = NULL
){

  g0 = frasyr::plot_vpa(vpa_sam_list)

  data = g0$data
  # data$stat %>% unique()
  # data$id %>% unique()
  data2 = data %>% dplyr::filter(stat %in% what.plot) %>%
    group_by(id,year,stat) %>%
    summarise(value=mean(value)) %>%
    ungroup() %>%
    mutate(value = if_else(stat == "fishing_mortality"|stat == "U",value,value/1000))
  # data2 %>% filter(stat == "U")
  # class(data2$stat)
  data2 = data2 %>%
    mutate(stat = as.character(stat)) %>%
    mutate(stat2 = case_when(stat=="biomass" ~ "Biomass",
                             stat=="fishing_mortality"~"F",
                             stat=="U"~"Exploitation_rate",
                             TRUE ~ stat))
  what.plot_f = case_when(what.plot=="biomass" ~ "Biomass",
                          what.plot=="fishing_mortality"~"F",
                          what.plot=="U"~"Exploitation_rate",
                          TRUE ~ what.plot)
  # data2$stat2 %>% unique
  data2 = data2 %>%
    # mutate(model = if_else(id=="1","VPA","SAM")) %>%
    mutate(stat_f = factor(stat2,levels=what.plot_f,labels=what.plot_f))
  # data2$stat_f %>% unique()
  # nrow(data2)
  if (!is.null(scenario_name)) {
    data2 = data2 %>% mutate(model = scenario_name[sapply(1:nrow(data2), function(i) which(data2$id[i]==unique(data2$id)))])
  } else{
    scenario_name = unique(as.character(data2$id))
    # data2 = data2 %>% mutate(model = id)
    data2 = data2 %>% mutate(model = scenario_name[sapply(1:nrow(data2), function(i) which(data2$id[i]==unique(data2$id)))])
  }

  CVdata_all = tibble()
  stat_order = case_when(what.plot=="biomass" ~ "Biomass",
                          what.plot=="fishing_mortality"~"F",
                          what.plot=="U"~"Exploitation_rate",
                          TRUE ~ what.plot)
  for(i in 1:length(vpa_sam_list)) {
    res = vpa_sam_list[[i]]
    if(class(res)=="vpa") {
      if (is.null(res$rep)) {
        stop("Rerun vpa() with TMB=TRUE & sdreport=TRUE!")
      }
      cvdata = tibble(stat0 = names(res$rep$value),CV=res$rep$sd/res$rep$value,model=scenario_name[i])
      if (!("U" %in% what.plot)) cvdata = cvdata %>% filter(stat0 != "U")
      if (!("fishing_mortality" %in% what.plot)) cvdata = cvdata %>% filter(stat0 != "F_mean")
      cvdata = cvdata %>%
        mutate(Year = rep(as.numeric(colnames(res$naa)),nrow(cvdata)/ncol(res$naa)))

      cvdata_R = cvdata %>% filter(stat0=="N") %>%
        arrange(Year) %>%
        mutate(Age = rep(as.numeric(rownames(res$naa)),ncol(res$naa))) %>%
        filter(Age == 0) %>%
        select(-Age) %>%
        mutate(stat = "Recruitment")

      # cvdata$stat0 %>% unique()
      cvdata = cvdata %>% filter(stat0 %in% c("SSB","B_total","F_mean","U")) %>%
        mutate(stat = case_when(stat0=="SSB" ~ "SSB",stat0=="F_mean" ~ "F", stat0=="U" ~ "Exploitation_rate",TRUE ~ "Biomass")) %>%
        full_join(cvdata_R) %>%
        mutate(stat_f = factor(stat,levels=stat_order)) %>%
        select(-stat0,stat)
      CVdata_all = bind_rows(CVdata_all,cvdata)
    }

    if (class(res)=="sam") {
      if (is.null(res$rep$unbiased)){
        cvdata2 = tibble(stat0 = names(res$rep$value),CV=res$rep$sd/res$rep$value,model=scenario_name[i]) %>%
          filter(stat0 %in% c("exp_logN","ssb","B_total","F_mean","Exploitation_rate"))
      }else{
        cvdata2 = tibble(stat0 = names(res$rep$value),CV=res$rep$sd/res$rep$unbiased$value,model=scenario_name[i]) %>%
          filter(stat0 %in% c("exp_logN","ssb","B_total","F_mean","Exploitation_rate"))
      }

      cvdata2_R = cvdata2 %>% filter(stat0=="exp_logN") %>%
        mutate(Age = as.numeric(rep(rownames(res$naa),ncol(res$naa)))) %>%
        filter(Age == 0) %>%
        mutate(Year = as.numeric(colnames(res$naa))) %>%
        select(-Age) %>% mutate(stat = "Recruitment")

      cvdata2 = cvdata2 %>% filter(stat0 %in% c("ssb","B_total","F_mean","Exploitation_rate")) %>%
        mutate(Year = rep(as.numeric(colnames(res$naa)),4)) %>%
        mutate(stat = case_when(stat0=="ssb" ~ "SSB",stat0=="F_mean" ~ "F", stat0=="B_total" ~ "Biomass", TRUE~stat0)) %>%
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

  data3 = full_join(data2,CVdata_all) %>%
    mutate(stat_f = factor(stat_f,levels=stat_order)) %>%
    mutate(Cz = exp(qnorm(CI+(1-CI)/2)*sqrt(log(1+CV^2)))) %>%
    mutate(lower = value/Cz, upper = value*Cz) %>%
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
    geom_path(aes(colour=Model,linetype=Model),size=size)+
    facet_wrap(vars(stat_f),scales="free_y",ncol=2)+
    theme_SH()+theme_bw(base_size=base_size)+theme(legend.position=legend_position)+
    xlab("Year") + ylab("")+
    # ylim(0,NA)
    scale_colour_brewer(palette="Set1",name=legend_name)+
    scale_linetype_discrete(name=legend_name)+
    guides(colour=guide_legend(nrow=legend_nrow),
           fill=guide_legend(nrow=legend_nrow),
           linetype=guide_legend(nrow=legend_nrow))
  if (isTRUE(log_scale)) {
    g1 = g1 + scale_y_log10()
  } else {
    g1 = g1 + ylim(0,NA)
  }

  g1
}

#' レトロスペクティブ解析の結果をプロットする
#' @param res フルデータを使った解析結果 (vpa or sam)
#' @param retro_res
#' @import dplyr
#' @export

retro_plot = function(res,retro_res,start_year=NULL,scale=1000, forecast=FALSE,
                      base_size=16, plot_mohn=TRUE, mohn_position = "upperleft") {

  mohn_res = retro_res$mohn
  if (isTRUE(forecast)) {
    mohn_res = retro_res$mohn_forecast
    if (class(res)=="vpa") warning("'forecast=TRUE' is not possible for VPA")
  }

  res_tibble = convert_sam_tibble(res) %>% mutate(id = 0)
  maxyear = res_tibble$year %>% max
  if (class(res)=="vpa" & isTRUE(res$input$last.catch.zero)) res_tibble = res_tibble %>% dplyr::filter(year<maxyear)

  for (i in 1:length(retro_res$Res)) {
    if (class(res)=="vpa" & isTRUE(res$input$last.catch.zero)) {
      res_tibble <- full_join(res_tibble, convert_sam_tibble(retro_res$Res[[i]]) %>%
                                mutate(id = i) %>% dplyr::filter(year<maxyear-i))
    } else {
      res_tibble <- full_join(res_tibble, convert_sam_tibble(retro_res$Res[[i]]) %>% mutate(id = i))
    }
  }

  res_tibble -> data
  if (is.null(start_year)) start_year <- max(as.numeric(colnames(res$faa)))-15
  data2 = data %>% dplyr::filter(stat %in% c("SSB","biomass","Recruitment","fishing_mortality")) %>%
    group_by(id,year,stat) %>%
    summarise(value=mean(value)) %>%
    ungroup() %>%
    mutate(value = if_else(stat == "fishing_mortality",value,value/scale)) %>%
    mutate(terminal_year = max(data$year)-id)
  if (isTRUE(forecast) & class(res)=="sam") {
    data2 = data2 %>% mutate(terminal_year = terminal_year+1)  }
  data2 = data2 %>% dplyr::filter(year <= terminal_year) %>%
    mutate(stat_f = factor(stat,levels=c("biomass","SSB","Recruitment","fishing_mortality"),labels=c("Biomass","SSB","Recruitment","F"))) %>%
    filter(year >= start_year) %>%
    filter(year <= as.numeric(terminal_year)) %>%
    mutate(terminal_year = factor(terminal_year, levels=unique(terminal_year)))
  g1 = ggplot(data=data2,aes(x=year,y=value))+
    geom_path(aes(colour=terminal_year,linetype=terminal_year),size=1)+
    facet_wrap(vars(stat_f),scales="free_y",ncol=2)+
    theme_SH()+theme_bw(base_size=base_size)+theme(legend.position="none")+
    xlab("Year") + ylab("")+
    ylim(0,NA)
  if (isTRUE(plot_mohn)) {
    mohn = tibble(rho = mohn_res, stat = names(mohn_res)) %>%
      dplyr::filter(stat != "N") %>%
      mutate(stat_f = factor(stat,levels=c("B","SSB","R","F"),labels=c("Biomass","SSB","Recruitment","F"))) %>%
      mutate(label=sprintf("rho == %.2f",rho)) %>%
      full_join(data2 %>% group_by(stat_f) %>% summarise(ymax = max(value)))
    if (mohn_position=="upperleft") {
      g1 = g1 + geom_text(data=mohn,parse=TRUE,aes(x=start_year,y=ymax,label=label,hjust=0,vjust=1))
    } else {
      end_year = data2$year %>% max
      g1 = g1 + geom_text(data=mohn,parse=TRUE,aes(x=end_year,y=0,label=label,hjust=1,vjust=0))
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
    geom_path(data=index_pred,aes(y=pred,colour=Model),size=1)+
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
    geom_path(data=index_pred_curve,aes(y=pred,colour=Model),size=1)+
    geom_point(data=index_pred2,aes(y=obs,colour=Model),size=1.5)+
    facet_wrap(vars(Fleet),nrow=2,scales=scales[3])+
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

  caa_pred = as_tibble(samres$caa) %>% mutate(Age=0:(n()-1)+samres$input$rec.age) %>%
    pivot_longer(cols=-Age,names_to="Year",values_to="pred") %>%
    mutate(Year = as.numeric(Year))

  maxage = caa_obs$Age %>% max

  caa_dat = full_join(caa_obs,caa_pred) %>%
    mutate(resid = log(obs/pred)) %>%
    mutate(Age = if_else(Age==maxage,str_c("Age ",Age,"+"),str_c("Age ",Age)))

  g_caa = ggplot(data=caa_dat,aes(x=Year))+
    geom_point(aes(y=obs),size=1.5) +
    geom_path(aes(y=pred),size=0.8)+
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
